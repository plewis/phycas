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

#if !defined(RELATIVE_RATE_DISTRIBUTION_HPP)
#define RELATIVE_RATE_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include "ncl/nxsdefs.h"

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>

#include "subset_proportions.hpp"
#include "dirichlet_distribution.hpp"

namespace phycas
{

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	This is the distribution of a relative rates vector X = (x1, x2, ..., xn) when the vector Y = (p1 x1, p2 x2, ..., pn xn) ~ Dirichlet(c1, c2, ..., cn). It is
|	useful as a prior for relative rates of subsets in a partition model, where the coefficients p1, p2, ..., pn are the relative sizes of the subsets.
*/
class RelativeRateDistribution : public DirichletDistribution
	{
	public:
											RelativeRateDistribution();
											RelativeRateDistribution(const std::vector<double> & params, const double_vect_t & coeffs);
                        					RelativeRateDistribution(const RelativeRateDistribution & other);
		virtual								~RelativeRateDistribution() {}

        RelativeRateDistribution * 			cloneAndSetLot(Lot * other) const;
        RelativeRateDistribution * 			Clone() const;

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
		// virtual void 								SetCoefficients(const std::vector<double> & coeff);
		virtual unsigned							GetNParams() const;
#		if defined(PYTHON_ONLY)
		//void								AltSetMeanAndVariance(std::vector<double> m, std::vector<double> v);
		double_vect_t						GetVarCovarMatrix();
#		endif

        void                                setSubsetProportions(SubsetProportionsShPtr subset_proportions);

    private:

		unsigned							dim;                    /**< The dimension, which equals the number of parameters (used to initialze the coefficients) */
        SubsetProportionsShPtr              _subset_proportions;    /**< The coefficients used to weight the relative rates */
		double								sum_params;             /**< The sum of the dirichlet parameters stored in the data member `dirParams', which is provided by the base class */
	};

typedef boost::shared_ptr<RelativeRateDistribution> RelativeRateDistributionShPtr;

} // namespace phycas

#endif
