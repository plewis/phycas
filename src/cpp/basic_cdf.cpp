/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#include "basic_cdf.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the gamma function for the value `a'.
*/
double CDF::LnGamma(double a) const
	{
	double alpha = a;
	return gamln(&alpha);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns x where the integral of the beta probability density function from 0 to x equals the supplied value `p'
*/
double CDF::BetaQuantile(
  double p,		/**< is the integral of the probability density function from 0 to x (the value x is returned) */
  double alpha,	/**< is the first parameter of the Beta distribution */
  double beta)	const /**< is the second parameter of the Beta distribution */
	{
    int which    = 2;       // tells cdfbet function you want to calculate X,Y given P,Q,A,B
	double P	 = p;       // integral from 0 to X
	double Q	 = 1.0 - p; // integral from X to 1
	double A	 = alpha;   // first shape parameter of the beta distribution
	double B	 = beta;    // second shape parameter of the beta distribution
    double X     = 0.0;     // upper limit of integration (calculated and returned)
    double Y     = 0.0;     // lower limit of integration (calculated but not returned)
    int status   = 0;       // remains 0 if everything went well
    double bound = 0.0;     // unimportant if status is 0

    cdfbet(&which, &P, &Q, &X, &Y, &A, &B, &status, &bound);
	return X;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative beta distribution function for value `x' and parameters `alpha' and `beta'. This is the integral
|	of the beta probability density function from 0 to x.
*/
double CDF::CumBeta(
  double x,		/**< is the integral of the probability density function from 0 to x will be returned */
  double alpha,	/**< is the first parameter of the Beta distribution */
  double beta)	const /**< is the second parameter of the Beta distribution */
	{
	double p	= x;
	double q	= 1.0 - x;
	double a	= alpha;
	double b	= beta;
	double cum	= 0.0;
	double ccum	= 0.0;

	cumbet(&p, &q, &a, &b, &cum, &ccum);
	return cum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative gamma distribution function for the value `x', shape parameter `alpha' and scale parameter
|	`beta'. This is the integral of the gamma probability density function from 0 up to `x'. The gamma density as
|	defined here is
|>
|	       x^{alpha - 1} e^{-x/beta}
|	f(x) = -------------------------
|	       beta^{alpha} Gamma(alpha)
|>
|	Note: the CDF library for which this CDF class is a wrapper defines the scale parameter as the inverse of beta!
*/
double CDF::CumGamma(
  double x,					/**< is the upper limit of the interval starting at 0.0 for which the integral will be computed */
  double alpha,				/**< is the shape parameter of the gamma distribution */
  double beta) const		/**< is the scale parameter of the gamma distribution */
	{
	double X		= x;
	int status		= 0;
	int which		= 1;	// compute p given x, shape and scale
	double P		= 0.0;
	double Q		= 0.0;
	double shape	= alpha;
	double scale	= 1.0/beta; // CDF library uses inverse of normal definition of scale parameter
	double bound	= 0.0;
	cdfgam(&which, &P, &Q, &X, &shape, &scale, &status, &bound);
	PHYCAS_ASSERT(status == 0);
	return P;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `x' such that the integral of the gamma probability density function from 0 up to `x' is equal to
|	the supplied parameter `p' (with shape parameter `alpha' and scale parameter `beta' provided). Useful for simulating
|	draws from a gamma distribution. The gamma density as defined here is
|>
|	       x^{alpha - 1} e^{-x/beta}
|	f(x) = -------------------------
|	       beta^{alpha} Gamma(alpha)
|>
|	Note: the CDF library for which this CDF class is a wrapper defines the scale parameter as the inverse of beta!
*/
double CDF::SampleGamma(
  double p,				/**< is the integral of the gamma probability density function from 0 up to `x' */
  double alpha,			/**< is the shape parameter of the gamma distribution */
  double beta) const	/**< is the scale parameter of the gamma distribution */
	{
	double X = 0.0;
	int status = 0;
	int which = 2;	// compute x given p, shape and scale
	double P = p;
	double Q = 1.0 - P;
	double shape = alpha;
	double scale = 1.0 / beta; // CDF library uses inverse of normal definition of scale parameter
	double bound = 0.0;
	cdfgam(&which, &P, &Q, &X, &shape, &scale, &status, &bound);
	if (status == 10)
		{
		// If alpha is tiny (alpha < 0.0019) and beta is correspondingly huge, then cdfgam may
		// fail with status 10, in which case we return 0.0
		PHYCAS_ASSERT(alpha < 0.0019 && beta > 1.0/0.0019);
		return 0.0;
		}
	return X;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative normal distribution function for the value `x', mean `mu' and standard deviation `sigma'. This
|	is the integral of the normal probability density function from negative infinity up to `x'. The normal density is
|>
|	                 1             /  (x - mu)^2   \
|	f(x) = ------------------- exp| --------------- |
|	        sigma (2 pi)^{0.5}     \   2 sigma^2   /
|>
*/
double CDF::CumNorm(
  double x,					/**< is the upper limit of the interval starting at 0.0 for which the integral will be computed */
  double mu,				/**< is the mean of the normal distribution */
  double sigma) const		/**< is the standard deviation of the distribution */
	{
	PHYCAS_ASSERT(sigma > 0.0);
	double X		= (x - mu)/sigma;
	int status		= 0;
	int which		= 1;	// compute p given x, mean and standard deviation
	double P		= 0.0;
	double Q		= 0.0;
	double bound	= 0.0;
	cdfnor(&which, &P, &Q, &X, &mu, &sigma, &status, &bound);
	PHYCAS_ASSERT(status == 0);
	return P;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `x' such that the integral of the normal probability density function from negative infinity up to
|	`x' is equal to the supplied parameter `p' (with mean `mu' and standard devition `sigma' provided). Useful for
|	simulating draws from a normal distribution. The normal density is
|>
|	                 1             /  (x - mu)^2   \
|	f(x) = ------------------- exp| --------------- |
|	        sigma (2 pi)^{0.5}     \   2 sigma^2   /
|>
*/
double CDF::SampleNorm(
  double p,				/**< is the integral of the normal probability density function from negative infinity up to `x' */
  double mu,			/**< is the mean of the normal distribution */
  double sigma) const	/**< is the standard deviation of the normal distribution */
	{
	double X = 0.0;
	int status = 0;
	int which = 2;	// compute x given p, shape and scale
	double P = p;
	double Q = 1.0 - P;
	double bound = 0.0;
	cdfnor(&which, &P, &Q, &X, &mu, &sigma, &status, &bound);
	PHYCAS_ASSERT(status == 0);
	return X;
	}
