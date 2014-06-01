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

#include "lognormal.hpp"
#include "ncl/nxsexception.h"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `mean' to 0.0 and `sd' to 1.0. See initialize function for details.
*/
LognormalDistribution::LognormalDistribution()
  	{
	initialize(0.0, 1.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the `logmean' and `logsd' data members to `logm' and `logs', respectively. See initialize function for
|	details.
*/
LognormalDistribution::LognormalDistribution(
  double logm,		/* the log mean parameter */
  double logs)		/* the log standard deviation parameter */
  	{
	initialize(logm, logs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the `logmean' and `logsd' data members to `other.logmean' and `other.logsd', respectively. See the
|	initialize function for details.
*/
LognormalDistribution::LognormalDistribution(
  const LognormalDistribution & other)	/* the lognormal distribution to clone */
  	{
	initialize(other.logmean, other.logsd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
LognormalDistribution::~LognormalDistribution()
	{
	//std::cerr << "Deleting a LognormalDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the `logmean' and `logsd' data members to `logm' and `logs', respectively. Also initializes `pi_const'
|	and `sqrt2_const'. The following is relevant to the calculation of `pi_const':
|>
|	      /|
|	     / |     sin(theta) = y/r
|	  r /  |     sin(90 deg)  = sin(pi/2) = 1.0
|	   /   | y   asin(1.0) = pi/2
|	  /    |     2.0*asin(1.0) = pi
|	 /_____|
|	    x
|>
*/
void LognormalDistribution::initialize(
  double logm,		/**< is the mean of the log of the lognormal random variable */
  double logs)		/**< is the standard deviation of the log of the lognormal random variable */
  	{
	PHYCAS_ASSERT(logs > 0.0);
	logmean = logm;
	logsd = logs;
	pi_const = 2.0*std::asin(1.0);
	sqrt2_const = std::sqrt(2.0);
	ComputeLnConst();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object, calls the new object's SetLot member function (passing the
|	supplied Lot object `other'), and returns a pointer to it. The caller is expected to manage the new object.
*/
LognormalDistribution * LognormalDistribution::cloneAndSetLot(Lot * other) const
	{
    LognormalDistribution * clone = new LognormalDistribution(logmean, logsd);
	clone->SetLot(other);
	return clone;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object.
*/
LognormalDistribution * LognormalDistribution::Clone() const
	{
    return new LognormalDistribution(logmean, logsd);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false because the univariate lognormal is a continuous distribution.
*/
bool LognormalDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Lognormal", which is the name of this distribution.
*/
std::string LognormalDistribution::GetDistributionName() const
	{
	return "Lognormal";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Lognormal(<logmean>, <logstddev>)", with <logmean> and <logstddev> replaced with the current
|	values of the `logmean' and `logsd' data members.
*/
std::string LognormalDistribution::GetDistributionDescription() const
	{
	return boost::str(boost::format("Lognormal(%#.5f, %#.5f)") % logmean % logsd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected value of the lognormal distribution as currently specified. The mean of a lognormal(logmean,
|	logsd) random variable is exp(logmean + logsd^2/2).
*/
double LognormalDistribution::GetMean() const
	{
	return std::exp(logmean + logsd_squared/2.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the lognormal distribution as currently specified. The variance of a lognormal(logmean,
|	logsd) random variable is (exp(logsd^2) - 1)*exp(2*logmean + logsd^2)

*/
double LognormalDistribution::GetVar() const
	{
	double term1 = (std::exp(logsd_squared) - 1.0);
	double term2 = std::exp(2.0*logmean + logsd_squared);
	return term1*term2;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected standard deviation of the lognormal distribution as currently specified, which is simply the
|	square root of the value returned by GetVar().
*/
double LognormalDistribution::GetStdDev() const
	{
	return std::sqrt(GetVar());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative lognormal distribution for x (integral of lognormal density function from negative infinity to x).
*/
double LognormalDistribution::GetCDF(
  double x)	 const	/**< is the value for which the cumulative distribution function is to be evaluated */
	{
	double term = (std::log(x) - logmean)/logsd;
	return cdf.CumNorm(term, 0.0, 1.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from this lognormal distribution.
*/
double LognormalDistribution::Sample() const
	{
	double z = cdf.SampleNorm(lot->Uniform(), 0.0, 1.0);
	double v = logsd*z + logmean;
	return std::exp(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value returned by GetLnPDF function.
*/
double LognormalDistribution::GetRelativeLnPDF(
  double x) const	/**< the value for which the density function is to be evaluated */
	{
	return GetLnPDF(x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the normal distribution is
|>
|	                 1              /  (log(x) - logmean)^2 \
|	f(x) = -------------------- exp| ----------------------  |
|	       x logsd (2 pi)^{0.5}     \       2 logsd^2       /
|>
|	This function returns the natural log of the density function at `x':
|>
|	log[f(x)] = (log(x) - logmean)^2/(2 logsd^2) - log(x) - log(logsd) - 0.5 log(2 pi)
|<
|	The sum of the last two terms is precalculated and available as the variable `ln_const'.
*/
double LognormalDistribution::GetLnPDF(
  double x)   const /**< the value for which the density function is to be evaluated */
	{
	double logx = std::log(x);
	double term1 = std::pow(logx - logmean, 2.0);
	double term2 = 2.0*logsd_squared;
	double term3 = term1/term2;
	double lnpdf = -term3 - logx + ln_const;
	return lnpdf;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the `logmean' and `logsd'
|	parameters.	The relationships below...
|>
|	m = exp(logmean + logsd_squared/2)
|   v = [exp(logsd_squared) - 1]*exp(2*logmean + logsd_squared)
|>
|	allow us to estimate logmean and logsd_squared (and hence logsd) from m and v:
|>
|	logsd_squared = log(v/m^2 + 1)
|	logmean       = log(m) - logsd_squared/2
|>
*/
void LognormalDistribution::SetMeanAndVariance(
  double m, 	/* the desired mean of the lognormal distribution */
  double v)		/* the desired variance of the lognormal distribution */
  	{
	PHYCAS_ASSERT(m > 0.0);
	PHYCAS_ASSERT(v > 0.0);
	double m_squared = pow(m, 2.0);
	logsd_squared = std::log(v/m_squared + 1.0);
	logmean = std::log(m) - logsd_squared/2.0;
	logsd = std::sqrt(logsd_squared);
	ComputeLnConst();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the value of the data member `ln_const' (the log of that part of the normal density function that depends
|	only on the standard deviation parameter, and which is recalculated whenever `sd' changes).
*/
void LognormalDistribution::ComputeLnConst()
  	{
	logsd_squared = std::pow(logsd, 2.0);
	ln_const = -std::log(logsd) - 0.5*std::log(2.0*pi_const);
	}

} // namespace phycas
