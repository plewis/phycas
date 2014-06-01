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

#if !defined(LOGNORMAL_HPP)
#define LOGNORMAL_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
//#include "ncl/nxsdefs.h"

#include <boost/shared_ptr.hpp>
//#include <boost/format.hpp>

#include "probability_distribution.hpp"

namespace phycas
{

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	The Lognormal distribution with two parameters, the mean and standard deviation of the natural logarithm of the variable.
*/
class LognormalDistribution : public ProbabilityDistribution
	{
	public:
					LognormalDistribution();
					LognormalDistribution(double mean, double stddev);
					LognormalDistribution(const LognormalDistribution & other);
					~LognormalDistribution();

        LognormalDistribution * cloneAndSetLot(Lot * other) const;
        LognormalDistribution * Clone() const;
		bool		IsDiscrete() const;
		std::string	GetDistributionName() const;
		std::string	GetDistributionDescription() const;
		double		GetMean() const;
		double		GetVar() const;
		double		GetStdDev() const;
		double		GetCDF(double x) const;
		double		Sample() const;
		double		GetLnPDF(double x) const;
		double		GetRelativeLnPDF(double x) const;
		void		SetMeanAndVariance(double mean, double var);

	protected:
		void 			initialize(double m, double s);
		void			ComputeLnConst();

	protected:
		double		logmean;			/**< the mean parameter of the lognormal distribution */
		double		logsd;				/**< the standard deviation parameter of the lognormal distribution */
		double		logsd_squared;		/**< the square of the logsd parameter (i.e. the variance of the log of this lognormal random variable) */
		double		ln_const;			/**< the natural logarithm of the constant part of the density function */
		double		pi_const;			/**< precalculated (in constructor) value of pi */
		double		sqrt2_const;		/**< precalculated (in constructor) value of sqrt(2.0) */

	};

} // namespace phycas

#endif
