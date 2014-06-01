/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2012 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#if !defined(TREE_LENGTH_DISTRIBUTION_HPP)
#define TREE_LENGTH_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include "ncl/nxsdefs.h"

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include "basic_cdf.hpp"
#include "basic_lot.hpp"
#include "basic_tree.hpp"
#include "probability_distribution.hpp"
#include "multivariate_probability_distribution.hpp"

namespace phycas
{

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	This is the joint distribution of tree length and edge lengths. The marginal distribution of tree length is Gamma(alphaT, betaT), while the marginal
|   of edge lengths conditional on tree length is Dirichlet(a_1, a_2, ..., a_m, b_1, b_2, ..., b_n) where a_1=a_2=...=a_m=alpha, b_1=b_2=...=b_n=c*alpha,
|   m is the number of external edges, and n is the number of internal edges. This prior was described in the following paper:
|
|   Rannala, Bruce, Tianqi Zhu, and Ziheng Yang. 2012. Tail paradox, partial identifiability, and influential priors in Bayesian branch length inference.
|   Mol. Biol. Evol. 29(1):325-335.
|
|   Note 1: this class follows the conventions of the paper in that the gamma distribution mean equals shape/scale, whereas everywhere else a gamma distribution
|   appears in Phycas it has the property that mean equals shape*scale.
|
|   Note 2: this class is unusual in that it represents a probability distribution but is not derived from either ProbabilityDistribution or
|   MultivariateProbabilityDistribution. This is because it is a compound distribution combining a univariate (gamma) distribution with a multivariate
|   (Dirichlet) distribution, and requires a tree as input to its GetLnPDF function.
*/
class TreeLengthDistribution   // TREE_LENGTH_DISTRIBUTION
	{
	public:
											TreeLengthDistribution();
											TreeLengthDistribution(double alphaT, double betaT, double alpha, double c);
                        					TreeLengthDistribution(const TreeLengthDistribution & other);
		virtual								~TreeLengthDistribution() {}

		virtual void						SetLot(Lot * other);
		virtual void						ResetLot();
		virtual void						SetSeed(unsigned rnseed);

        TreeLengthDistribution * 			cloneAndSetLot(Lot * other) const;
        TreeLengthDistribution * 			Clone() const;

		virtual std::string                 GetDistributionName() const;
		virtual std::string 				GetDistributionDescription() const;
		virtual std::vector<double>			Sample(unsigned num_external, unsigned num_internal);
		virtual double						GetLnPDF(TreeShPtr t) const;
		virtual double						GetRelativeLnPDF(TreeShPtr t) const;

        double                              getShape() const {return _alphaT;}
        double                              getScale() const {return _betaT;}
        double                              getExtEdgelenParam() const {return _alpha;}
        double                              getIntExtEdgelenRatio() const {return _c;}

    protected:

        void                                SetupSamplingDistributions(unsigned num_external, unsigned num_internal);

    private:

		double								_alphaT;                /**< The shape parameter of the gamma tree length distribution */
		double								_betaT;                 /**< The scale paramter of the gamma tree length distribution */
		double								_alpha;                 /**< The parameter governing the Dirichlet distribution for external edge lengths (internal edge lengths have parameter c*alpha) */
		double								_c;                     /**< The ratio of the mean internal/external edge length (normally less than 1.0) */

        CDF                                 _cdf;                   /**< Used in GetLnPDF to compute log of gamma function */
		Lot                                 _myLot;                 /**< Own random number generator */
		Lot *                               _lot;                   /**< Points to either _myLot or an external random number generator object */

        ProbDistShPtr                       _tldist;                /**< Points to a gamma distribution of tree lengths used when sampling from this distribution */
        MultivarProbDistShPtr               _eldist;                /**< Points to a Dirichlet distribution of edge length proportions used when sampling from this distribution */
        unsigned                            _num_internal_edges;    /**< The number of internal edges specified the last time Sample was called */
        unsigned                            _num_external_edges;    /**< The number of external edges specified the last time Sample was called */
	};

typedef boost::shared_ptr<TreeLengthDistribution> TreeLengthDistributionShPtr;

} // namespace phycas

#endif

