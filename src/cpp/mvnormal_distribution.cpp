/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2009 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#include "mvnormal_distribution.hpp"
#include "ncl/nxsexception.h"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Construct a standard, uncorrelated bivariate normal distribution by default.
*/
MVNormalDistribution::MVNormalDistribution() 
	{
	mean.resize(2, 0.0);
	varcov.Initialize(2, 2);
	varcov.ptr[0][0] = 1.0;
	varcov.ptr[0][1] = 0.0;
	varcov.ptr[1][0] = 0.0;
	varcov.ptr[1][1] = 1.0;
	InitMVNorm();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a MVNormalDistribution based on the mean vector `meanvec' and variance-covariance matrix `varcovmat'.
*/
MVNormalDistribution::MVNormalDistribution(const std::vector<double> & meanvec, const std::vector<double> & varcovmat) 
	{
	unsigned dim = (unsigned)meanvec.size();
	
	//std::cerr << "dim from meanvec = " << dim << std::endl;

	// Reserve necessary space
	mean.resize(dim, 0.0);
	varcov.Initialize(dim, dim);

	// Copy meanvec to mean and varcovmat to varcov
	for (unsigned i = 0; i < dim; ++i)
		{
		//std::cerr << "meanvec[" << i << "] = " << meanvec[i] << std::endl;
		mean[i] = meanvec[i];
		for (unsigned j = 0; j < dim; ++j)
			{
			varcov.ptr[i][j] = varcovmat[i*dim + j];
			//std::cerr << "varcov[" << i << "][" << j << "] = " << varcov.ptr[i][j] << std::endl;
			}
		}
	//std::cerr << "mean.size()   = " << mean.size() << std::endl;
	InitMVNorm();	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `mean' and `varcov' data members to the values contained in their counterparts in `other'.
*/
MVNormalDistribution::MVNormalDistribution(
  const MVNormalDistribution & other)	/* the multivariate normal distribution to clone */
	{
	unsigned dim = (unsigned)other.mean.size();
	
	// Reserve necessary space
	mean.resize(dim, 0.0);
	varcov.Initialize(dim, dim);

	// Copy other.mean to mean and other.varcov to varcov
	for (unsigned i = 0; i < dim; ++i)
		{
		mean[i] = other.mean[i];
		for (unsigned j = 0; j < dim; ++j)
			{
			varcov.ptr[i][j] = other.varcov.ptr[i][j];
			}
		}
	InitMVNorm();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object, calls the new object's SetLot member function (passing the 
|	supplied Lot object `other'), and returns a pointer to it. The caller is expected to manage the new object. 
*/
MVNormalDistribution * MVNormalDistribution::cloneAndSetLot(Lot * other) const
	{
    MVNormalDistribution * clone = new MVNormalDistribution(*this);
	clone->SetLot(other);
	return clone;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
MVNormalDistribution * MVNormalDistribution::Clone() const
	{
    return new MVNormalDistribution(*this);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Fits this MVNormal distribution to the supplied `data', which should be in the form of a one-dimensional vector of
|	double values. This vector will be interpreted as a two-dimensional data matrix with `nrows' rows (observations)
|	and `ncols' columns (variables). Returns true if successful and false if there were problems (e.g. the 
|	variance-covariance matrix could not be inverted). If this function returns false, this MVNormalDistribution
|	object will be left unmodified.
*/
bool MVNormalDistribution::Fit(
  unsigned nrows, 			/**< is the number of rows (observations) */
  unsigned ncols, 			/**< is the number of columns (variables) */
  std::vector<double> data)	/**< is the flattened (1st. row starts at index 0, 2nd. row starts at index ncols, etc.) two-dimensional array */
	{
	PHYCAS_ASSERT(nrows*ncols == data.size());
	PHYCAS_ASSERT(nrows > 2);	// must be at least 3 observations
	PHYCAS_ASSERT(ncols > 1);   // if not at least 2 variables, should use NormalDistribution instead
	
	unsigned row, col, i, j, k;
	
	double mean_denom = (double)nrows;
	double varcov_denom = (double)(nrows - 1);

	// copy current mean vector in case revert necessary
	unsigned dim = (unsigned)mean.size();
	std::vector<double> m(dim, 0.0);
	std::copy(mean.begin(), mean.end(), m.begin());
		
	// copy current varcov matrix in case revert necessary
	ScopedTwoDMatrix<double> s;
	s.Initialize(dim, dim);
	for (i = 0; i < dim; ++i)
		for (j = 0; j < dim; ++j)
			{
			s.ptr[i][j] = varcov.ptr[i][j];
			}

	// resize mean and varcov and zero out both
	mean.resize(ncols);
	varcov.Initialize(ncols, ncols);
	for (i = 0; i < ncols; ++i)
		{
		mean[i] = 0.0;
		for (j = 0; j < ncols; ++j)
			varcov.ptr[i][j] = 0.0;
		}
	
	// estimate mean from data
	unsigned n = (unsigned)data.size();
	for (i = 0; i < n; ++i)
		{
		col = i % ncols;
		mean[col] += data[i]/mean_denom;
		}
		
	// estimate covariance matrix from data
	for (row = 0; row < nrows; ++row)
		{
		k = row*ncols;
		for (i = 0; i < ncols; ++i)
			{
			double di = data[k + i] - mean[i];
			for (j = 0; j < ncols; ++j)
				{
				double dj = data[k + j] - mean[j];
				varcov.ptr[i][j] += di*dj/varcov_denom;
				}
			}
		}
		
	bool ok = InitMVNorm();

	if (!ok)
		{
		// revert to previous mean and varcov
		mean.resize(dim);
		varcov.Initialize(dim, dim);
		for (i = 0; i < ncols; ++i)
			{
			mean[i] = m[i];
			for (j = 0; j < ncols; ++j)
				varcov.ptr[i][j] = s.ptr[i][j];
			}
		}
		
	return ok;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets parameters of the distribution from the vector of means and vector of variances and covariances. If dim is the
|	length of the vector `m', then the ith. dim elements of `v' represent the ith. row of the variance/covariance
|	matrix.
*/
void MVNormalDistribution::SetMeanAndVariance(const std::vector<double> & m, const std::vector<double> & v)
	{
	unsigned dim = (unsigned)m.size();
	PHYCAS_ASSERT(dim*dim == v.size());
	
	// Reserve necessary space
	mean.resize(dim, 0.0);
	varcov.Initialize(dim, dim);

	// Copy other.mean to mean and other.varcov to varcov
	unsigned k = 0;
	for (unsigned i = 0; i < dim; ++i)
		{
		mean[i] = m[i];
		for (unsigned j = 0; j < dim; ++j)
			{
			varcov.ptr[i][j] = v[k++];
			}
		}
	InitMVNorm();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Set varcov[i][j] = varcov[j][i] = (varcov[i][j] + varcov[j][i])/2.0 to ensure that the variance-covariance matrix
|	is symmetric.
*/
void MVNormalDistribution::MakeVarCovSymmetric() 
	{
	unsigned dim = (unsigned)mean.size();
	for (unsigned i = 0; i < dim - 1; ++i)
		{
		for (unsigned j = i + 1; j < dim; ++j)
			{
			double v = (varcov.ptr[i][j] + varcov.ptr[j][i])/2.0;
			varcov.ptr[i][j] = v;
			varcov.ptr[j][i] = v;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assuming the `mean' and `varcov' data members are populated with values, this function computes the matrix L 
|	(Cholesky decomposition such that varcov = L*Ltranspose), which is used for sampling from the distribution.
|	The determinant of `varcov' can be obtained as the square of the product of diagonal elements of L, the Cholesky
|	decomposition, precalculated in InitMVNorm(). If varcov cannot be inverted, returns false; otherwise returns true.
*/
bool MVNormalDistribution::InitMVNorm() 
	{
	bool ok = true;
	lot = &myLot;
	
	unsigned i, j, k;
	unsigned dim = (unsigned)mean.size();

	//std::cerr << "dim = " << dim << std::endl;

	MakeVarCovSymmetric();
	
	// set up standard normal distribution
	stdnorm.SetMeanAndVariance(0.0, 1.0);
	stdnorm.SetLot(lot);
	
	L.Initialize(dim,dim);
	
	// Using Cholesky algorithm described in http://en.wikipedia.org/wiki/Cholesky_decomposition
	for (j = 0; j < dim; ++j)
		{
		double v = varcov.ptr[j][j];
		L.ptr[j][j] = sqrt(v);
		//std::cerr << boost::str(boost::format("L[%d][%d] = %.5f") % j % j % L.ptr[j][j]) << std::endl;
		for (i = 0; i < j; ++i)
			{
			v = L.ptr[i][j];
			L.ptr[j][j] -= v*v;
			//std::cerr << boost::str(boost::format("L[%d][%d] -= %.5f") % j % j % (v*v)) << std::endl;
			}
		for (i = j + 1; i < dim; ++i)
			{
			L.ptr[j][i] = 0.0;
			L.ptr[i][j] = varcov.ptr[i][j];
			//std::cerr << boost::str(boost::format("L[%d][%d] = %.5f") % i % j % L.ptr[i][j]) << std::endl;
			for (k = 0; k < j; ++k)
				{
				L.ptr[i][j] -= L.ptr[i][k]*L.ptr[j][k];
				//std::cerr << boost::str(boost::format("L[%d][%d] -= (%.5f)*(%.5f)") % i % j % L.ptr[i][k] % L.ptr[j][k]) << std::endl;
				}
			L.ptr[i][j] /= L.ptr[j][j];
			//std::cerr << boost::str(boost::format("L[%d][%d] /= %.5f") % i % j % L.ptr[j][j]) << std::endl;
			}
		}
		
	// compute the inverse of the variance-covariance matrix
	std::vector<double> workvec1(dim, 0.0);
	std::vector<int> workvec2(dim, 0);
	invVarCov.Initialize(dim, dim);
	ScopedTwoDMatrix<double> workmatrix;
	workmatrix.Initialize(dim, dim);
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			double d = varcov.ptr[i][j];
			invVarCov.ptr[i][j] = d;
			workmatrix.ptr[i][j] = d;
			}
		}
	int result = InvertMatrix(workmatrix.ptr, (int)dim, &workvec1[0], &workvec2[0], invVarCov.ptr);
	if (result != 0)
		ok = false;

	// compute the determinant of varcov using the diagonal elements of L
	double tmp = 1.0;
	for (i = 0; i < dim; ++i)
		{
		tmp *= L.ptr[i][i];
		}
	PHYCAS_ASSERT(tmp > 0.0);
	logDetVarCov = 2.0*log(tmp);
	//std::cerr << "logDetVarCov = " << logDetVarCov << std::endl;
	
	double pi = 4.0*std::atan(1.0);
	nLogTwoPi = (double)dim*log(2.0*pi);
	
	return ok;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assuming the `mean' and `varcov' data members are initialized, this function computes the matrices L and D (Cholesky
|	decomposition such that varcov = L*D*Ltranspose), which is used for sampling from the distribution.
*/
std::string MVNormalDistribution::DebugMVNorm() const
	{
	unsigned i, j, dim;
	dim = (unsigned)mean.size();
	
	std::string s = "mean:\n\n";
	for (i = 0; i < dim; ++i)
		s << boost::str(boost::format(" %12.5f") % mean[i]);
	s << "\nvarcov:\n\n";
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			s << boost::str(boost::format(" %12.5f") % varcov.ptr[i][j]);
		s << "\n";
		}
	s << "\nL:\n\n";
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			s << boost::str(boost::format(" %12.5f") % L.ptr[i][j]);
		s << "\n";
		}
	s << "\ninvVarCov:\n\n";
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			s << boost::str(boost::format(" %12.5f") % invVarCov.ptr[i][j]);
		s << "\n";
		}
	return s;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the MultivariateProbabilityDistribution::Sample member function. 
|	The original random number generator (data member `myLot') can be replaced by calling the 
|	MultivariateProbabilityDistribution::ResetLot function. Note that this object does not take ownership of the Lot 
|	object whose pointer is specified as `other'. It is assumed that `other' is non-NULL.
*/
void MVNormalDistribution::SetLot(
	Lot * other) /**< is a pointer to the random number generator object to be used subsequently by Sample */
	{
	if (other == NULL)
		throw XProbDist("attempt made to install a non-existent pseudorandom number generator in a MVNormalDistribution object");
	lot = other;
	stdnorm.SetLot(lot);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number seed for the `myLot' data member. Note that if MVNormalDistribution::SetLot has been 
|	called, calling MVNormalDistribution::SetSeed is pointless because you will not be setting the seed for the 
|	random number generator being used!
*/
void MVNormalDistribution::SetSeed(
  unsigned rnseed)	/**< is the new seed value */
	{
	myLot.SetSeed(rnseed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes the data member `lot' (which is used as the random number generator by the member function 
|	MultivariateProbabilityDistribution::Sample) point to the local data member `myLot'. This function only needs to be
|	called if MultivariateProbabilityDistribution::SetLot has been called previously to replace the random number 
|	generator used by MultivariateProbabilityDistribution::Sample. 
*/
void MVNormalDistribution::ResetLot()
	{
	lot = &myLot;
	stdnorm.SetLot(lot);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false because the multivariate normal is a continuous distribution.
*/
bool MVNormalDistribution::IsDiscrete() const 
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns "Bivariate Normal" if the dimension is 2, or "Multivariate Normal" otherwise.
*/
std::string MVNormalDistribution::GetDistributionName() const 
	{
	if (mean.size() == 2)
		return "Bivariate Normal";
	return "Multivariate Normal";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string that, if printed, shows how to construct an object just like this one. Because of their length,
|	the mean vector and variance-covariance matrix are shown separately above the constructor call.
*/
std::string MVNormalDistribution::GetDistributionDescription() const 
	{
	unsigned dim = (unsigned)mean.size();
	PHYCAS_ASSERT(dim > 1);
	std::string s;
	s << "mu = (\n  ";
	for (unsigned i = 0; i < dim; ++i) 
		{
		s << boost::str(boost::format("%#.5f") % mean[i]);
		if (i < dim - 1)
			s << ",\n  ";
		else
			s << "\n  )\n";
		}
	s << "Sigma = (\n  ";
	for (unsigned i = 0; i < dim; ++i) 
		{
		s << "(";
		for (unsigned j = 0; j < dim; ++j) 
			{
			s << boost::str(boost::format("%#.5f") % varcov.ptr[i][j]);
			if (j < dim - 1)
				s << ", ";
			}
		if (i < dim - 1)
			s << "),\n  ";
		else
			s << ")\n  ";
		}
	s << ")\nMVNormal(mu, Sigma)";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is needed because in Python the parameters of the distribution are supplied to the constructor via
|	a tuple, which requires an extra set of parentheses. In this case, we can just return the result provided by 
|	GetDistributionDescription because, due to their length, the arguments have been split out from the function call
|	and thus there should be no confusion about the extra parentheses.
*/
std::string MVNormalDistribution::GetDescriptionForPython() const 
	{
	return GetDistributionDescription() ;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a copy of the mean vector.
*/
std::vector<double> MVNormalDistribution::GetMean() const 
	{
	std::vector<double> mean_copy(mean.size());
	std::copy(mean.begin(), mean.end(), mean_copy.begin());
	return mean_copy;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a copy of the variance-covariance matrix as a vector (first row, second row, etc).
*/
std::vector<double> MVNormalDistribution::GetVar() const 
	{
	unsigned dim = (unsigned)mean.size();
	std::vector<double> varcov_copy(dim*dim);
	unsigned k = 0;
	for (unsigned i = 0; i < dim; ++i)
		for (unsigned j = 0; j < dim; ++j)
			varcov_copy[k++] = varcov.ptr[i][j];
	return varcov_copy;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a copy of the square-root of the variance-covariance matrix as a vector (first row, second row, etc). The
|	square root function is applied separately to each element of the varcov matrix.
*/
std::vector<double> MVNormalDistribution::GetStdDev() const 
	{
	unsigned dim = (unsigned)mean.size();
	std::vector<double> stdev(dim*dim);
	unsigned k = 0;
	for (unsigned i = 0; i < dim; ++i)
		for (unsigned j = 0; j < dim; ++j)
			stdev[k++] = sqrt(varcov.ptr[i][j]);
	return stdev;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Approximates the cumulative distribution function evaluated at the supplied point `x'. The precision of the 
|	approximation is controlled by `nsamples'. The approximation is done using a brute force approach: `nsamples'
|	samples are drawn from this multivariate normal distribution, and the proportion of those samples that are inside
|	the region defined by `x' is returned as the approximated CDF. The supplied point should be a vector of length equal
|	to the number of parameters of the mean vector.
*/
double MVNormalDistribution::ApproxCDF(
  const std::vector<double> & x,	/**< is the point to be evaluated */
  unsigned nsamples)				/**< is the number of samples to use in approximating the CDF */
  const 
	{
	PHYCAS_ASSERT(nsamples > 0);
	unsigned dim = GetNParams();
	PHYCAS_ASSERT(x.size() == dim);
	unsigned num_inside = 0;
	unsigned i, j;
	std::vector<double> y(dim, 0.0);

	for (unsigned k = 0; k < nsamples; ++k)
		{
		// draw a random sample (this part is identical to body of Sample() function)
		std::copy(mean.begin(), mean.end(), y.begin());
		for (j = 0; j < dim; ++j)
			{
			double z = stdnorm.Sample();
			for (i = j; i < dim; ++i)
				{
				y[i] += z*L.ptr[i][j];
				}
			}
			
		// check if sample less than x
		bool y_less_than_x = true;
		for (i = 0; i < dim; ++i)
			{
			if (y[i] > x[i])
				{
				y_less_than_x = false;
				break;
				}
			}
		
		if (y_less_than_x)
			num_inside++;
		}

	return (double)num_inside/(double)nsamples;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Samples from this multivariate normal distribution by drawing a vector Z of independent, univariate standard normal
|	variates and then returning the vector mean + L*Z, where L is the Cholesky decomposition computed in InitMVNorm().
*/
std::vector<double> MVNormalDistribution::Sample() const 
	{
	unsigned i, j;
	unsigned dim = (unsigned)mean.size();
	
	/*
	
	 /  \     /           \   /  \      /                             \
	| m1 |   | L11  0   0  | | z1 |    | m1 + L11*z1                   |
	| m2 | + | L21 L22  0  | | z2 | =  | m2 + L21*z1 + L22*z2          |
	| m3 |   | L31 L32 L33 | | z3 |    | m3 + L31*z1 + L32*z2 + L33*z3 |
	 \  /     \           /   \  /      \                             /
	 
	*/
	
	std::vector<double> x(dim, 0.0);
	std::copy(mean.begin(), mean.end(), x.begin());
	for (j = 0; j < dim; ++j)
		{
		double z = stdnorm.Sample();
		for (i = j; i < dim; ++i)
			{
			x[i] += z*L.ptr[i][j];
			}
		}

	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the multivariate normal density of the supplied vector `x'. Here is the density function in pseudo-latex:
|	
|	f(x) = (2 \pi)^{-n/2} \det(varcov)^{-0.5} exp{ -0.5 (x - mean)^T varcov^{-1} (x-mean)}
|	
|	The log of f(x) is thus:
|	
|	log[f(x)] = -0.5 \left[ n \log(2 \pi) + \log(\det(varcov)) + (x - mean)^T varcov^{-1} (x-mean) \right]
|	
*/
double MVNormalDistribution::GetLnPDF(
  const std::vector<double> & x)
  const 
	{
	unsigned dim = (unsigned)mean.size();
	PHYCAS_ASSERT(dim == (unsigned)x.size());
	unsigned i, j;
	
	// compute vector y = x - mean
	std::vector<double> y(dim, 0.0);	
	for (i = 0; i < dim; ++i)
		{
		y[i] = x[i] - mean[i];
		}
		
	double log_density = 0.0;
	
	// perform multiplication y^T varcov^{-1} y
	for (i = 0; i < dim; ++i)
		{
		double tmp = 0.0;
		for (j = 0; j < dim; ++j)
			{
			tmp += y[j]*invVarCov.ptr[j][i];
			}
		tmp *= y[i];
		log_density += tmp;
		}
	
	// finish computation of log_density
	log_density += nLogTwoPi + logDetVarCov;
	log_density *= -0.5;
	
	return log_density;
	}

double MVNormalDistribution::GetRelativeLnPDF(
  const std::vector<double> & x)
  const 
	{
	return GetLnPDF(x);
	}

#if defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a 1-dimensional vector (first row, then second row, etc.)
*/
std::vector<double> MVNormalDistribution::GetVarCovarMatrix()
	{
	unsigned i, j, k;
	unsigned dim = (unsigned)mean.size();
	std::vector<double> V(dim*dim, 0.0);
	k = 0;
	for (i = 0; i < dim; ++i)
		for (j = 0; j < dim; ++j)
			V[k++] = varcov.ptr[i][j];
	return V;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of parameters.
*/
unsigned MVNormalDistribution::GetNParams() const 
	{
	return (unsigned)mean.size();
	}

} // namespace phycas
