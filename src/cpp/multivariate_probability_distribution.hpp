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

#if !defined(MULTIVARIATE_PROBABILITY_DISTRIBUTION_HPP)
#define MULTIVARIATE_PROBABILITY_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include "ncl/nxsdefs.h"

#include "states_patterns.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include "probability_distribution.hpp"
#include "basic_cdf.hpp"
#include "basic_lot.hpp"
#include "phycas_string.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/tuple.hpp>
#	include <boost/python/numeric.hpp>
#	include "thirdparty/num_util/num_util.h"
#endif
#include "xprobdist.hpp"
//class XUnderflow{};

namespace phycas
{

class MultivariateProbabilityDistribution
	{
	public:
									MultivariateProbabilityDistribution() {lot = &myLot;}
							virtual ~MultivariateProbabilityDistribution() {}

		virtual void				SetLot(Lot *other)									= 0;
		virtual void				ResetLot()													= 0;
		virtual void				SetSeed(unsigned rnseed)									= 0;

		virtual	bool				IsDiscrete() const											= 0;
		virtual std::string 		GetDistributionName() const									= 0;
		virtual std::string 		GetDistributionDescription() const							= 0;
		virtual std::vector<double>				GetMean() const												= 0;
		virtual std::vector<double> 				GetVar() const												= 0;
		virtual std::vector<double> 				GetStdDev() const											= 0;
		virtual std::vector<double>				Sample() const												= 0;
		virtual double				ApproxCDF(const std::vector<double> &x, unsigned nsamples = 10000) const	= 0;
		virtual double				GetLnPDF(const std::vector<double> &) const								= 0;
		virtual double				GetRelativeLnPDF(const std::vector<double> &) const						= 0;
		virtual void 				SetMeanAndVariance(const std::vector<double> &m, const std::vector<double> &v)		= 0;
#		if defined(PYTHON_ONLY)
#			if defined(USING_NUMARRAY)
				virtual void 			AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v) = 0;
				virtual boost::python::numeric::array	GetVarCovarMatrix()						= 0;
#			else
				//virtual void 			AltSetMeanAndVariance(std::vector<double> m, std::vector<double> v)				= 0;
				virtual std::vector<double>			GetVarCovarMatrix()										= 0;
#			endif
#		endif
		virtual unsigned			GetNParams()	const										= 0;

		CDF					cdf;
		Lot					myLot;
		Lot *				lot;
	};

typedef boost::shared_ptr<MultivariateProbabilityDistribution> MultivarProbDistShPtr;

} // namespace phycas

#endif
