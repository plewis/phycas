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

#if ! defined(PRIOR_DISTRIBUTION_HPP)
#define PRIOR_DISTRIBUTION_HPP

#include "probability_distribution.hpp"
#include "multivariate_probability_distribution.hpp"

namespace phycas
{

class Tree;
typedef boost::shared_ptr<Tree> TreeShPtr;

class TopoProbCalculator;
typedef boost::shared_ptr<TopoProbCalculator> TopoProbCalculatorShPtr;

class TreeLengthDistribution;
typedef boost::shared_ptr<TreeLengthDistribution> TreeLengthDistributionShPtr;

class JointPriorManager;

/*----------------------------------------------------------------------------------------------------------------------
|   Stores a univariate or multivariate probability distribution along with other information needed to comply with a
|   request for the density at the current value of the parameter (or parameter vector). This is a helper class for use
|   with JointPriorManager.
*/
class PriorDistribution
	{
    friend class JointPriorManager;

	public:
                                    PriorDistribution();
                                    ~PriorDistribution();

        void                        initUnivariate(std::string name, ProbDistShPtr dist, double value);
        void                        initMultivariate(std::string name, MultivarProbDistShPtr dist, const double_vect_t & value);

        void                        initExternalEdgeLenDistribution(std::string name, ProbDistShPtr dist, TreeShPtr tree);
        void                        initInternalEdgeLenDistribution(std::string name, ProbDistShPtr dist, TreeShPtr tree);

        void                        initTreeLenDistribution(std::string name, TreeLengthDistributionShPtr dist, TreeShPtr tree);
        void                        initTopologyDistribution(std::string name, TopoProbCalculatorShPtr topo_prior_calculator, TreeShPtr tree);

        void                        debugShowCurrValue(std::ostream & out) const;

    private:

        enum                        edgelen_t {edgelen_none, edgelen_all, edgelen_internal, edgelen_external};
        edgelen_t                   _edgelen_type;

        bool                        _neverBeforeInitialized();
        double                      _recalculateLogDensity();

        std::string                 _name;

        ProbDistShPtr               _univar_dist;
        MultivarProbDistShPtr       _multivar_dist;
        TopoProbCalculatorShPtr     _topology_dist;
        TreeLengthDistributionShPtr _treelen_dist;

        double                      _univar_value;
        double_vect_t               _multivar_value;
        TreeShPtr                   _tree;

        double                      _current_density;
    };

typedef boost::shared_ptr<PriorDistribution>    PriorDistrShPtr;

} // namespace phycas

#endif
