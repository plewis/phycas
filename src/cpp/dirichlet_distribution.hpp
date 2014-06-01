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

#if !defined(DIRICHLET_DISTRIBUTION_HPP)
#define DIRICHLET_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

//#include <cmath>
//#include "ncl/nxsdefs.h"

//#include "phycas/src/states_patterns.hpp"

//#include <boost/shared_ptr.hpp>
//#include <boost/format.hpp>
#include "multivariate_probability_distribution.hpp"
//#include "phycas/src/basic_cdf.hpp"
//#include "phycas/src/basic_lot.hpp"
//#include "phycas/src/phycas_string.hpp"
//#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
//#	include <boost/python/tuple.hpp>
//#	include <boost/python/numeric.hpp>
//#	include "thirdparty/num_util/num_util.h"
//#endif
//#include "phycas/src/xprobdist.hpp"

namespace phycas
{

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	The Dirichlet distribution, with the number of parameters determined by the number of params supplied to the constructor.
*/
class DirichletDistribution : public MultivariateProbabilityDistribution
	{
	public:
                                                    DirichletDistribution();
                                                    DirichletDistribution(const std::vector<double> & params);
                                                    DirichletDistribution(const DirichletDistribution & other);
		virtual                                     ~DirichletDistribution();

        DirichletDistribution *                     cloneAndSetLot(Lot * other) const;
        DirichletDistribution *                     Clone() const;
		virtual void								SetLot(Lot * other);
		virtual void								ResetLot();
		virtual void								SetSeed(unsigned rnseed);

		virtual bool								IsDiscrete() const;
		virtual std::string 						GetDistributionName() const;
		virtual std::string 						GetDistributionDescription() const;
		std::string                                 GetDescriptionForPython() const;
		virtual std::vector<double>					GetMean() const;
		virtual std::vector<double> 				GetVar() const;
		virtual std::vector<double> 				GetStdDev() const;
		virtual std::vector<double>					Sample() const;
		virtual double								ApproxCDF(const std::vector<double> &x, unsigned nsamples = 10000) const;
		virtual double								GetLnPDF(const std::vector<double> &x) const;
		virtual double								GetRelativeLnPDF(const std::vector<double> &x) const;
		virtual void 								SetMeanAndVariance(const std::vector<double> &m, const std::vector<double> &v);
#		if defined(PYTHON_ONLY)
#			if defined(USING_NUMARRAY)
				void                                AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v);
				boost::python::numeric::array       GetVarCovarMatrix();
#			else
				void                                AltSetMeanAndVariance(std::vector<double> m, std::vector<double> v);
				std::vector<double>					GetVarCovarMatrix();
#			endif
#		endif

		virtual unsigned							GetNParams() const;
		const GammaDistribution &                   GetDistributionOnParameter(unsigned i) const;


    protected:

		void                                        initialize(const double_vect_t & params);

    protected:

		std::vector<double>                         dirParams;
		std::vector<GammaDistribution>              paramDistributions;
		mutable std::vector<double>                 scratchSpace;
	};

} // namespace phycas

#endif
