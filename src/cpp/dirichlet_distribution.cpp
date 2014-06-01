/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

//#include "ncl/nxsexception.h"
#include "dirichlet_distribution.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `dirParams', `scratchSpace' and `paramDistributions' vectors to the values contained in their
|   counterparts in `other'.
*/
DirichletDistribution::DirichletDistribution(
  const DirichletDistribution & other)	/* the dirichlet distribution to clone */
	{
	dirParams.resize(other.dirParams.size());
	std::copy(other.dirParams.begin(), other.dirParams.end(), dirParams.begin());

	scratchSpace.resize(other.scratchSpace.size());
	std::copy(other.scratchSpace.begin(), other.scratchSpace.end(), scratchSpace.begin());

    paramDistributions.clear();
    for (std::vector<GammaDistribution>::const_iterator it = other.paramDistributions.begin(); it != other.paramDistributions.end(); ++it)
	    paramDistributions.push_back(GammaDistribution(*it));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object, calls the new object's SetLot member function (passing the
|	supplied Lot object `other'), and returns a pointer to it. The caller is expected to manage the new object.
*/
DirichletDistribution * DirichletDistribution::cloneAndSetLot(Lot * other) const
	{
    DirichletDistribution * clone = new DirichletDistribution(dirParams);
	clone->SetLot(other);
	return clone;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object.
*/
DirichletDistribution * DirichletDistribution::Clone() const
	{
    return new DirichletDistribution(dirParams);
    }

const GammaDistribution &DirichletDistribution::GetDistributionOnParameter(
  unsigned i) const
  	{
  	PHYCAS_ASSERT( i < paramDistributions.size());
	return paramDistributions[i];
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	Construct the equivalent of a flat Beta distribution by default.
*/
DirichletDistribution::DirichletDistribution()
	{
    //std::cerr << "DirichletDistribution::DirichletDistribution()" << std::endl;
	dirParams.clear();
	scratchSpace.clear();
	paramDistributions.clear();
	dirParams.push_back(1.0);
	dirParams.push_back(1.0);
	scratchSpace.push_back(0.0);
	scratchSpace.push_back(0.0);
	paramDistributions.push_back(GammaDistribution(1.0, 1.0));
	paramDistributions.push_back(GammaDistribution(1.0, 1.0));
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a DirichletDistribution based on the parameter values supplied in the vector `params'.
*/
DirichletDistribution::DirichletDistribution(const std::vector<double> & params)
	{
    //std::cerr << "DirichletDistribution::DirichletDistribution(params)" << std::endl;
	initialize(params);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a DirichletDistribution based on the parameter values supplied in the vector `params'.
*/
DirichletDistribution::~DirichletDistribution()
    {
    //std::cerr << "DirichletDistribution::~DirichletDistribution()" << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes a DirichletDistribution based on the parameter values supplied in the vector `params'.
*/
void DirichletDistribution::initialize(
  const double_vect_t & params)	/**< is the vector of Dirichlet parameters, which should sum to 1.0 */
	{
	unsigned params_size = (unsigned)params.size();
	PHYCAS_ASSERT(params_size > 1);

	// Reserve necessary space to save copying
	dirParams.reserve(params_size);
	scratchSpace.reserve(params_size);
	paramDistributions.reserve(params_size);

#	if 0   //@POL Mark, this causes an "out of keys" compiler error on VC
#		if defined NCL_NXS_THROW_UNDEFINED
		if (params.size() < 2)
			throw NxsX_UndefinedException("Illegal Dirichlet", __FILE__, __LINE__);
#		endif
#	endif

	std::vector<double>::const_iterator iter;
	for (iter = params.begin(); iter != params.end(); ++iter)
		{
		double c_i = *iter;
		dirParams.push_back(c_i);
		scratchSpace.push_back(0.0);
    	paramDistributions.push_back(GammaDistribution(c_i, 1.0));
		}

	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the MultivariateProbabilityDistribution::Sample member function.
|	The original random number generator (data member `myLot') can be replaced by calling the
|	MultivariateProbabilityDistribution::ResetLot function. Note that this object does not take ownership of the Lot
|	object whose pointer is specified as `other'. It is assumed that `other' is non-NULL.
|	This overridden member function differs from the base class version ProbabilityDistribution::SetLot in that it
|	calls SetLot for all its component distributions so that the same random number generator is used throughout.
*/
void DirichletDistribution::SetLot(
	Lot * other) /**< is a pointer to the random number generator object to be used subsequently by Sample */
	{
	if (other == NULL)
		throw XProbDist("attempt made to install a non-existent pseudorandom number generator");
	lot = other;
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(lot);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number seed for the `myLot' data member. Note that if DirichletDistribution::SetLot has been
|	called, calling DirichletDistribution::SetSeed is pointless because you will not be setting the seed for the
|	correct random number generator!
*/
void DirichletDistribution::SetSeed(
  unsigned rnseed)	/**< is the new seed value */
	{
	myLot.SetSeed(rnseed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes the data member `lot' (which is used as the random number generator by the member function
|	MultivariateProbabilityDistribution::Sample) point to the local data member `myLot'. This function only needs to be
|	called if MultivariateProbabilityDistribution::SetLot has been called previously to replace the random number
|	generator used by MultivariateProbabilityDistribution::Sample. This overridden member function differs from the base
|	class version ProbabilityDistribution::ResetLot in that it calls SetLot for all its component distributions so that
|	the same random number generator is used throughout.
*/
void DirichletDistribution::ResetLot()
	{
	lot = &myLot;
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

bool DirichletDistribution::IsDiscrete() const
	{
	return false;
	}

std::string DirichletDistribution::GetDistributionName() const
	{
	if (GetNParams() == 2)
		return "Beta";
	return "Dirichlet";
	}

std::string DirichletDistribution::GetDistributionDescription() const
	{
	PHYCAS_ASSERT(GetNParams() > 1);
	std::string s;
	s << "Dirichlet(";
	//s << MakeStrPrintF("%.5f", dirParams[0]);
	s << str(boost::format("%#.5f") % dirParams[0]);
	for (unsigned i = 1; i < GetNParams(); ++i)
		{
		//s << MakeStrPrintF(",%.5f", dirParams[i]);
		s << str(boost::format(",%#.5f") % dirParams[i]);
		}
	s << ')';
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is needed because in Python the parameters of the distribution are supplied to the constructor via
|	a tuple, which requires an extra set of parentheses. For example, in Phycas, one could sample from a flat Dirichlet
|	as follows:
|
|	begin phycas;
|		sample dist=dirichlet(1.0,1.0,1.0,1.0);
|	end;
|
|	Note only one set of parentheses surround the vector of Dirichlet parameters. In Python, we have:
|
|	>>> DirichletDist((1.0,1.0,1.0,1.0)).sample()
|
|	Note the extra set of parentheses needed in the Python representation.
*/
std::string DirichletDistribution::GetDescriptionForPython() const
	{
	PHYCAS_ASSERT(GetNParams() > 1);
	std::string s;
	s << "Dirichlet((";
	//s << MakeStrPrintF("%.5f", dirParams[0]);
	s << str(boost::format("%#.5f") % dirParams[0]);
	for (unsigned i = 1; i < GetNParams(); ++i)
		{
		//s << MakeStrPrintF(",%.5f", dirParams[i]);
		s << str(boost::format(", %#.5f") % dirParams[i]);
		}
	s << "))";
	return s;
	}

std::vector<double> DirichletDistribution::GetMean() const
	{
	std::vector<double> retvect;

	double sum = 0.0;
	for (std::vector<double>::const_iterator dIt = dirParams.begin(); dIt != dirParams.end(); ++dIt)
		sum += *dIt;

	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		double x = dirParams[i]/sum;
		retvect.push_back(x);
		}

	return retvect;
	}

std::vector<double> DirichletDistribution::GetVar() const
	{
	std::vector<double> retvect;

	double c = 0.0;
	for (std::vector<double>::const_iterator dIt = dirParams.begin(); dIt != dirParams.end(); ++dIt)
		c += *dIt;

	double denom = c*c*(c + 1.0);

	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		double c_i = dirParams[i];
		double v = c_i*(c - c_i)/denom;
		retvect.push_back(v);
		}

	return retvect;
	}

std::vector<double> DirichletDistribution::GetStdDev() const
	{
	std::vector<double> retvect = GetVar();

	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		retvect[i] = sqrt(retvect[i]);
		}

	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Approximates the cumulative distribution function evaluated at the supplied point `x'. The precision of the
|	approximation is controlled by `nsamples'. The approximation is done using a brute force approach: `nsamples'
|	samples are drawn from this Dirichlet distribution, and the proportion of those samples that are inside the region
|	defined by `x' is returned as the approximated CDF. The supplied point should be a vector of length k, where k is
|	one fewer than the number of parameters of the Dirichlet distribution. If `x' has length greater than k, all
|	elements of `x' after the (k-1)th. will be ignored.
*/
double DirichletDistribution::ApproxCDF(
  const std::vector<double> &x,		/**< is the point to be evaluated */
  unsigned nsamples)	/**< is the number of samples to use in approximating the CDF */
  const
	{
	PHYCAS_ASSERT(nsamples > 0);
	unsigned nparams = GetNParams();
	PHYCAS_ASSERT(x.size() >= nparams - 1);
	unsigned num_inside = 0;

	for (unsigned k = 0; k < nsamples; ++k)
		{
		double sum = 0.0;
		for (unsigned i = 0; i < nparams; ++i)
			{
			const double y = paramDistributions[i].Sample();
			scratchSpace[i] = y;
			sum += y;
			}

		bool inside = true;
		for (unsigned j = 0; j < nparams - 1; ++j)
			{
			const double p = scratchSpace[j]/sum;
			if (p > x[j])
				{
				inside = false;
				break;
				}
			}

		if (inside)
			num_inside++;
		}

	return (double)num_inside/(double)nsamples;
	}

std::vector<double> DirichletDistribution::Sample() const
	{
	std::vector<double> x;

	double sum = 0.0;
	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		const double y = paramDistributions[i].Sample();
		scratchSpace[i] = y;
		sum += y;
		}

	for (unsigned j = 0; j < GetNParams(); ++j)
		x.push_back(scratchSpace[j]/sum);

	return x;
	}

double DirichletDistribution::GetLnPDF(
  const std::vector<double> &x)
  const
	{
	PHYCAS_ASSERT(GetNParams() == (unsigned)x.size());
	double retVal = 0.0;

	double c = 0.0;
	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		double c_i = dirParams[i];
		c += c_i;
		retVal += (c_i - 1.0) * std::log(x[i]);
		retVal -= cdf.LnGamma(c_i);
		}

	retVal += cdf.LnGamma(c);

	return retVal;
	}

double DirichletDistribution::GetRelativeLnPDF(
  const std::vector<double> &x)
  const
	{
	PHYCAS_ASSERT(GetNParams() == (unsigned)x.size());
	double retVal = 0.0;

	for (unsigned i = 0; i < GetNParams(); ++i)
		retVal += (dirParams[i] - 1.0) * std::log(x[i]);

	return retVal;
	}

#if defined(PYTHON_ONLY)
#if defined(USING_NUMARRAY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a numarray object.
*/
boost::python::numeric::array DirichletDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;
	unsigned dim = (unsigned)dirParams.size();

	std::vector<int> dims;
	dims.push_back((int)dim);
	dims.push_back((int)dim);

	double c = 0.0;
	for (i = 0; i < dim; ++i)
		{
		c += dirParams[i];
		}
	double denom = c*c*(c + 1.0);

	std::vector<double> V;
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			if (i == j)
				{
				double var_i = dirParams[i]*(c - dirParams[i])/denom;
				V.push_back(var_i);
				}
			else
				{
				double cov_ij = -dirParams[i]*dirParams[j]/denom;
				V.push_back(cov_ij);
				}
			}
		}

	return num_util::makeNum(&V[0], dims);
	}
#else	// USING_NUMARRAY not #defined
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a 1-dimensional vector (first row, then second row, etc.)
*/
std::vector<double> DirichletDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;
	unsigned dim = (unsigned)dirParams.size();

	double c = 0.0;
	for (i = 0; i < dim; ++i)
		{
		c += dirParams[i];
		}
	double denom = c*c*(c + 1.0);

	std::vector<double> V;
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			if (i == j)
				{
				double var_i = dirParams[i]*(c - dirParams[i])/denom;
				V.push_back(var_i);
				}
			else
				{
				double cov_ij = -dirParams[i]*dirParams[j]/denom;
				V.push_back(cov_ij);
				}
			}
		}

	return V;
	}
#endif
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Sets parameters of the distribution from the vector of means and vector of variances (note: v is not the variance
|	covariance matrix, but just a single-dimensional array of variances - the covariances are not needed).
*/
void DirichletDistribution::SetMeanAndVariance(const std::vector<double> &m, const std::vector<double> &v)
	{
	unsigned m_length = (unsigned)m.size();

	// check to make sure mean and variance vectors are the same length
	PHYCAS_ASSERT(m_length == (unsigned)v.size());

	// check to make sure user isn't trying to change the dimension of the distribution
	PHYCAS_ASSERT(m_length == GetNParams());

	// create iterators
	// get pointers to the underlying C arrays
	double *mean     = (double *)(&m[0]);
	double *variance = (double *)(&v[0]);

	// calculate the sum of variances and the sum of squared means
	unsigned i;
	double sum_of_variances = 0.0;
	double sum_of_squared_means = 0.0;
	for (i = 0; i < m_length; ++i)
		{
		sum_of_squared_means += (mean[i]*mean[i]);
		sum_of_variances += variance[i];
		}

	// calculate c, the sum of parameters
	if (sum_of_variances <= 0.0)
		throw XProbDist("sum of variances supplied to DirichletDist.setMeanAndVariance must be positive");
	double c = (1.0 - sum_of_squared_means - sum_of_variances)/sum_of_variances;

	// set parameters of the distribution
	dirParams.clear();
	paramDistributions.clear();
	for (unsigned i = 0; i < m_length; ++i)
		{
		double c_i = c*mean[i];
		dirParams.push_back(c_i);
		paramDistributions.push_back(GammaDistribution(c_i, 1.0));
		}
	}

#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
/*----------------------------------------------------------------------------------------------------------------------
|	Identical to SetMeanAndVariance, but uses numarray objects rather than std::vector. This is temporary, only serving
|	as a template for how to import a numarray and access its data within a member function.
*/
void DirichletDistribution::AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v)
	{
	unsigned m_length = (unsigned)m.nelements();

	// check to make sure mean and variance vectors are the same length
	PHYCAS_ASSERT(m_length == (unsigned)v.nelements());

	// check to make sure user isn't trying to change the dimension of the distribution
	PHYCAS_ASSERT(m_length == GetNParams());

	// get pointers to the underlying C arrays
	double *mean     = (double *)num_util::data(m);
	double *variance = (double *)num_util::data(v);

	// calculate the sum of variances and the sum of squared means
	unsigned i;
	double sum_of_variances = 0.0;
	double sum_of_squared_means = 0.0;
	for (i = 0; i < m_length; ++i)
		{
		sum_of_squared_means += (mean[i]*mean[i]);
		sum_of_variances += variance[i];
		}

	// calculate c, the sum of parameters
	if (sum_of_variances <= 0.0)
		throw XProbDist("sum of variances supplied to DirichletDist.setMeanAndVariance must be positive");
	double c = (1.0 - sum_of_squared_means - sum_of_variances)/sum_of_variances;

	// set parameters of the distribution
	dirParams.clear();
	paramDistributions.clear();
	for (i = 0; i < m_length; ++i)
		{
		double c_i = c*mean[i];
		dirParams.push_back(c_i);
		paramDistributions.push_back(GammaDistribution(c_i, 1.0));
		}
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of parameters.
*/
unsigned DirichletDistribution::GetNParams() const
	{
	return (unsigned)dirParams.size();
	}

} // namespace phycas
