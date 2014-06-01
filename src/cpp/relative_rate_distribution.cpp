/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2010 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	define NO_IMPORT_ARRAY
#endif

#include <numeric>
#include <functional>
#include <cmath>
#include <boost/lambda/lambda.hpp>

#include "multivariate_probability_distribution.hpp"
#include "relative_rate_distribution.hpp"
//#include "ncl/nxsexception.h"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Construct the equivalent of a flat Beta distribution by default.
*/
RelativeRateDistribution::RelativeRateDistribution()
  : DirichletDistribution(), dim(2)
	{
	sum_params = std::accumulate(dirParams.begin(), dirParams.end(), 0.0);
    _subset_proportions = SubsetProportionsShPtr(new SubsetProportions());
    _subset_proportions->setSubsetProportions(std::vector<double>(2, 0.5));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a RelativeRateDistribution based on the parameter values supplied in the vector `params' and the
|   coefficients supplied in the vector 'coeffs'.
*/
RelativeRateDistribution::RelativeRateDistribution(
  const double_vect_t & params, /**< is the vector of parameters */
  const double_vect_t & coeffs) /**< is the vector of coefficients */
  : DirichletDistribution(params), dim(2)
	{
	PHYCAS_ASSERT(params.size() > 1);
	dim = (unsigned)params.size();
	PHYCAS_ASSERT(coeffs.size() == dim);
	sum_params = std::accumulate(dirParams.begin(), dirParams.end(), 0.0);
    _subset_proportions = SubsetProportionsShPtr(new SubsetProportions());
    _subset_proportions->setSubsetProportions(coeffs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `dirParams', `scratchSpace' and `paramDistributions' vectors to the values contained in their
|   counterparts in `other'.
*/
RelativeRateDistribution::RelativeRateDistribution(
  const RelativeRateDistribution & other)	/* the relative rate distribution to clone */
	{
	dirParams.resize(other.dirParams.size());
	std::copy(other.dirParams.begin(), other.dirParams.end(), dirParams.begin());

	scratchSpace.resize(other.scratchSpace.size());
	std::copy(other.scratchSpace.begin(), other.scratchSpace.end(), scratchSpace.begin());

    paramDistributions.clear();
    for (std::vector<GammaDistribution>::const_iterator it = other.paramDistributions.begin(); it != other.paramDistributions.end(); ++it)
	    paramDistributions.push_back(GammaDistribution(*it));

	dim = other.dim;
	sum_params = other.sum_params;
    _subset_proportions = other._subset_proportions;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the string "RelativeRateDistribution".
*/
std::string RelativeRateDistribution::GetDistributionName() const
	{
	return "RelativeRateDistribution";
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a string similar to "RelativeRateDistribution(1.0,1.0,1.0,0.33333,0.33333,0.33333)", where the first 3
|   arguments are the underlying parameters while the last 3 are the probabilities associated with each relative rate.
*/
std::string RelativeRateDistribution::GetDistributionDescription() const
	{
	PHYCAS_ASSERT(dim > 1);
	std::string s;
    double_vect_t::const_iterator it;
	s << "RelativeRateDistribution(";
	for (it = dirParams.begin(); it != dirParams.end(); ++it)
		{
		if (it == dirParams.begin())
			s << boost::str(boost::format("%#.5f") % (*it));
		else
			s << boost::str(boost::format(",%#.5f") % (*it));
		}
    std::vector<double> const & p = _subset_proportions->getSubsetProportions();
	for (it = p.begin(); it != p.end(); ++it)
		{
		if (it == p.begin())
			s << boost::str(boost::format("%#.5f") % (*it));
		else
			s << boost::str(boost::format(",%#.5f") % (*it));
		}
	s << ')';
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is needed because in Python the parameters of the distribution are supplied to the constructor via
|	a tuple, which requires extra sets of parentheses. For example, in Phycas, one could sample from a 5-parameter
|	relative rate distribution as follows:
|
|	>>> RelativeRateDistribution((0.8,0.9,1.0,1.1,1.2),(0.2,0.2,0.2,0.2,0.2)).sample()
|
|	Note the parentheses around the parameters and coefficients needed in the Python representation.
*/
std::string RelativeRateDistribution::GetDescriptionForPython() const
	{
	PHYCAS_ASSERT(dim > 1);
	std::string s;
    double_vect_t::const_iterator it;
	s << "RelativeRateDistribution((";
	for (it = dirParams.begin(); it != dirParams.end(); ++it)
		{
		if (it == dirParams.begin())
			s << boost::str(boost::format("%#.5f") % (*it));
		else
			s << boost::str(boost::format(",%#.5f") % (*it));
		}
	s << "),(";
    std::vector<double> const & p = _subset_proportions->getSubsetProportions();
	for (it = p.begin(); it != p.end(); ++it)
		{
		if (it == p.begin())
			s << boost::str(boost::format("%#.5f") % (*it));
		else
			s << boost::str(boost::format(",%#.5f") % (*it));
		}
	s << "))";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles representing the mean of this RelativeRateDistribution. The mean is `dim' times the
|	mean vector for a Dirichlet random variable with the same parameters.
*/
double_vect_t RelativeRateDistribution::GetMean() const
	{
    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t retvect(dim, 0.0);
	for (unsigned i = 0; i < dim; ++i)
		{
		retvect[i] = dirParams[i]/(p[i]*sum_params);
		}
	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles representing the variance of this RelativeRateDistribution. The variance is `dim'
|	squared times the variance vector for a Dirichlet random variable with the same parameters.
*/
double_vect_t RelativeRateDistribution::GetVar() const
	{
    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t retvect(dim, 0.0);
	double denom = sum_params*sum_params*(sum_params + 1.0);
	for (unsigned i = 0; i < dim; ++i)
		{
		retvect[i] = dirParams[i]*(sum_params - dirParams[i])/(denom*p[i]*p[i]);
		}
	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles representing the square root of this RelativeRateDistribution. The square root is `dim'
|	times the square root vector for a Dirichlet random variable with the same parameters.
*/
double_vect_t RelativeRateDistribution::GetStdDev() const
	{
    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t retvect(dim, 0.0);
	double denom = sum_params*sqrt(sum_params + 1.0);
	for (unsigned i = 0; i < dim; ++i)
		{
		double numer = dirParams[i]*(sum_params - dirParams[i]);
		retvect[i] = sqrt(numer)/(denom*p[i]);
		}
	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Approximates the cumulative distribution function evaluated at the supplied point `x'. The precision of the
|	approximation is controlled by `nsamples'. The approximation is done using a brute force approach: `nsamples'
|	samples are drawn from this RelativeRateDistribution, and the proportion of those samples that are inside the region
|	defined by `x' is returned as the approximated CDF. The supplied point should be a vector of length k, where k is
|	one fewer than the number of parameters of the RelativeRateDistribution. If `x' has length greater than k, all
|	elements of `x' after the (k-1)th. will be ignored.
*/
double RelativeRateDistribution::ApproxCDF(
  const double_vect_t & x,		/**< is the point to be evaluated */
  unsigned nsamples)			/**< is the number of samples to use in approximating the CDF */
  const
	{
	PHYCAS_ASSERT(nsamples > 0);
	PHYCAS_ASSERT(x.size() >= dim - 1);
    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t tmp(x.size(), 0.0);
	for (unsigned i = 0; i < dim; ++i)
		{
		tmp[i] = p[i]*x[i];
		}
	return DirichletDistribution::ApproxCDF(tmp, nsamples);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Samples and returns a `dim'-dimensional point from this distribution.
*/
double_vect_t RelativeRateDistribution::Sample() const
	{
    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t x(dim, 0.0);

	double sum = 0.0;
	for (unsigned i = 0; i < dim; ++i)
		{
		const double y = paramDistributions[i].Sample();
		scratchSpace[i] = y;
		sum += y;
		}

	for (unsigned i = 0; i < dim; ++i)
		{
		x[i] = scratchSpace[i]/(p[i]*sum);
		}
	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the probability density function evaluated at the supplied point `x'. To derive
|	the relative rate distribution density, imagine a random variable Y ~ Beta(a,b) and assume X=Y/p. X has a relative
|	rate distribution with density
|>
|	f_X(x) = f_Y(pX) dy/dx
|	dy/dx = p
|>
|	Thus, the density of X will be p times the Beta density evaluated at pX. This generalizes to the Dirichlet: the
|	density of a value x of a random variable X having an n-parameter relative rate distribution will be p_1 p_2 ... p_n
|	times the corresponding Dirichlet pdf evaluated at the point (p_1 x_1, p_2 x_2, ..., p_n x_n), where p_i is the
|	weight associated with subset i and the sum of the p_i values equals 1.
*/
double RelativeRateDistribution::GetLnPDF(
  const double_vect_t & x) const 	/**< is the point at which to evaluate the density */
	{
	unsigned sz = (unsigned)x.size();

	if (sz == dim || sz == dim - 1)
		{
		double_vect_t tmp(dim, 0.0);

        std::vector<double> const & p = _subset_proportions->getSubsetProportions();

		std::transform(x.begin(), x.end(), p.begin(), tmp.begin(), boost::lambda::_1*boost::lambda::_2);

#if 0
        //POL CHECK SUBSET RELRATES
        std::cerr << "\nIn RelativeRateDistribution::GetLnPDF:" << std::endl;
        std::cerr << "i" << "\t" << "p[i]" << "\t" << "x[i]" << "\t" << "p[i]*x[i]" << "\t" << "tmp[i]" << std::endl;
        double doof_sum = 0.0;
        for (unsigned i = 0; i < sz; ++i)
            {
            std::cerr << i << "\t" << p[i] << "\t" << x[i] << "\t" << (p[i]*x[i]) << "\t" << tmp[i] << std::endl;
            doof_sum += p[i]*x[i];
            }
        std::cerr << "sum p[i]*x[i] = " << doof_sum;
        std::cerr << std::endl;
#endif

		double sum_tmp = std::accumulate(tmp.begin(), tmp.end(), 0.0);
        double log_prod_p = _subset_proportions->getLogProdProportions();
		if (sz == dim - 1)
			{
			// last element of tmp needs to be calculated because supplied vector x was one element short
			tmp[sz] = 1.0 - sum_tmp;
			sum_tmp = 1.0;
			double ln_pdf = log_prod_p + DirichletDistribution::GetLnPDF(tmp);
			return ln_pdf;
			}
		if (fabs(sum_tmp - 1.0) > 0.001)
			throw XProbDist(boost::str(boost::format("Vector supplied to RelativeRateDistribution::GetLnPDF is not valid: sum of transformed elements was %g but was expected to be 1.0 (tolerance +-0.001)") % sum_tmp));
		double ln_pdf = log_prod_p + DirichletDistribution::GetLnPDF(tmp);
		return ln_pdf;
		}
	else
		throw XProbDist(boost::str(boost::format("Vector supplied to RelativeRateDistribution::GetLnPDF has dimension %d but dimension should be either %d or %d") % sz % dim % (dim-1)));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the probability density function evaluated at the supplied point `x'.
*/
double RelativeRateDistribution::GetRelativeLnPDF(
  const double_vect_t & x) const
	{
	return GetLnPDF(x);
	}

#if defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a 1-dimensional vector (first row, then second row, etc.)
*/
double_vect_t RelativeRateDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;

	double c = sum_params;
	double denom = c*c*(c + 1.0);

    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t V;
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			if (i == j)
				{
				double var_i = dirParams[i]*(c - dirParams[i])/(p[i]*p[i]*denom);
				V.push_back(var_i);
				}
			else
				{
				double cov_ij = -dirParams[i]*dirParams[j]/(p[i]*p[j]*denom);
				V.push_back(cov_ij);
				}
			}
		}

	return V;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Sets parameters of the distribution from the vector of means and vector of variances (note: v is not the variance
|	covariance matrix, but just a single-dimensional array of variances - the covariances are not needed). Specifically,
|	if X is a relative rate vector, and Yi = Xi*pi, then E[Yi] = pi*E[Xi] and Var(Yi) = pi^2*Var(Xi). Thus, simple
|	transformations of the means (E[Xi]) and variances (Var(Xi)) of the elements of X (supplied in the form of `m' and
|	`v') can be turned over to the DirichletDistribution::SetMeanAndVariance function to parameterize the underlying
|	Dirichlet base class data members.
*/
void RelativeRateDistribution::SetMeanAndVariance(
  const double_vect_t & m, 	/**< is the vector of means */
  const double_vect_t & v) 	/**< is the vector of variances */
	{
	if (dim != (unsigned)m.size())
		throw XProbDist(boost::str(boost::format("Mean vector supplied to SetMeanAndVariance has dimension %d but RelativeRateDistribution has dimension %d") % (unsigned)m.size() % dim));
	if (dim != (unsigned)v.size())
		throw XProbDist(boost::str(boost::format("Variance vector supplied to SetMeanAndVariance has dimension %d but RelativeRateDistribution has dimension %d") % (unsigned)v.size() % dim));

    double_vect_t const & p = _subset_proportions->getSubsetProportions();
	double_vect_t transformed_means(dim);
	std::transform(m.begin(), m.end(), p.begin(), transformed_means.begin(), boost::lambda::_1*boost::lambda::_2);

	double_vect_t transformed_variances(dim);
	std::transform(v.begin(), v.end(), p.begin(), transformed_variances.begin(), boost::lambda::_1*boost::lambda::_2*boost::lambda::_2);

	DirichletDistribution::SetMeanAndVariance(transformed_means, transformed_variances);
	sum_params = std::accumulate(dirParams.begin(), dirParams.end(), 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of parameters.
*/
unsigned RelativeRateDistribution::GetNParams() const
	{
	return dim;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object, calls the new object's SetLot member function (passing the
|	supplied Lot object `other'), and returns a pointer to it. The caller is expected to manage the new object.
*/
RelativeRateDistribution * RelativeRateDistribution::cloneAndSetLot(Lot * other) const
	{
    RelativeRateDistribution * clone = new RelativeRateDistribution(*this);
	clone->SetLot(other);
	return clone;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object.
*/
RelativeRateDistribution * RelativeRateDistribution::Clone() const
	{
    return new RelativeRateDistribution(*this);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `_subset_proportions' to the supplied `subset_proportions'.
*/
void RelativeRateDistribution::setSubsetProportions(SubsetProportionsShPtr subset_proportions)
	{
    _subset_proportions = subset_proportions;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Copies the supplied coefficiences in `coeff' to the data member `p'.
*/
//void RelativeRateDistribution::SetCoefficients(
//  const double_vect_t & coeff)	/**< is the new vector of coefficients */
//	{
//	// Throw exception if coeff has an incorrect number of elements
//	if (dim != (unsigned)coeff.size())
//		throw XProbDist(boost::str(boost::format("Coefficient vector supplied to RelativeRateDistribution::SetCoefficients has dimension %d but RelativeRateDistribution has dimension %d") % (unsigned)coeff.size() % dim));
//
//	// Throw exception if any elements of coeff are negative or zero
//	for (double_vect_t::const_iterator it = coeff.begin(); it != coeff.end(); ++it)
//		{
//		double x = *it;
//		if (x <= 0.0)
//			throw XProbDist(boost::str(boost::format("Coefficient vector supplied to RelativeRateDistribution::SetCoefficients has at least one element (%g) with value less than or equal to zero") % x));
//		}
//
//	// Don't assume that the elements of coeff add up to 1.0
//	double sum_coeff = std::accumulate(coeff.begin(), coeff.end(), 0.0);
//	std::transform(coeff.begin(), coeff.end(), p.begin(), boost::lambda::_1/sum_coeff);
//
//	// Recalculate log_prod_p
//	log_prod_p = 0.0;
//	for (unsigned i = 0; i < dim - 1; ++i)
//		{
//		log_prod_p += log(p[i]);
//		}
//    }


} // namespace phycas


