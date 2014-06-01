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

#if ! defined(JOINT_PRIOR_MANAGER_HPP)
#define JOINT_PRIOR_MANAGER_HPP

#include "prior_distribution.hpp"

namespace phycas
{

class Tree;
typedef boost::shared_ptr<Tree> TreeShPtr;

typedef std::vector<PriorDistrShPtr> prior_distr_vect_t;
typedef std::map< std::string, PriorDistrShPtr > prior_distr_map_t;

/*----------------------------------------------------------------------------------------------------------------------
|	Maintains a list of prior distribution objects, reference distribution objects, and current parameter values for all
|   model parameters. Can be asked to compute the joint prior or the joint reference density, and all parameters are
|   expected to update their current values here so that the joint prior or reference density reported is always
|   accurate. Also manages sampling from the prior or reference distribution when necessary.
|
|   Here is a list of standard names for parameters:
|       Note: subset index substituted for %d when it is a prefix, with first subset = 1
|       Note: node number substituted for %d when it is a suffix
|
|       Univariate:
|           Transition/transversion rate ratio:                     %d_kappa
|           Non-synonymous/synonymous rate ratio:                   %d_omega
|           Nucleotide frequency parameters:                        %d_freqA, %d_freqC, %d_freqG, %d_freqT
|           Codon frequency parameters:                             %d_freqAAA, %d_freqAAC, %d_freqAAG, %d_freqAAT, ..., $d_freqTTT
|           GTR exchangeability parameters:                         %d_rAC, %d_rAG, %d_rAT, %d_rCG, %d_rCT, %d_rGT
|           Internal edge length parameters in fixed tree:          intedge_%d
|           External edge length parameters in fixed tree:          extedge_%d
|           Shape parameter of gamma distributed relative rates:    %d_gamma_shape
|           Proportion of invariable sites:                         %d_pinvar
|
|       Multivariate:
|           State frequencies:                                      %d_state_freqs
|           GTR exchangeabilities:                                  %d_relrates
|           Subset relative rates:                                  subset_relrates
|
|       Topology:
|           Tree topology:                                          tree_topology
|
|       Edge lengths:
|           External edge lengths:                                  external_edgelen
|           Internal edge lengths:                                  internal_edgelen
|
|       Edge length hyperparameters:
|           External edge length hyperparameter:                    external_hyper
|           Internal edge length hyperparameter:                    internal_hyper
|           Universal edge length hyperparameter:                   edgelen_hyper
|
|       Tree length:                                                tree_length
*/
class JointPriorManager
	{
	public:
								JointPriorManager();
								~JointPriorManager();

        void                    skipNextDebugCheck();

        void                    recalcLogJointPrior();
        double                  getLogJointPrior() const;
        double                  getLogTopologyPrior() const;

        bool                    debugCheckLogJointPrior(std::string called_from) const;
        std::string             debugPriorBreakdown(const char * msg = 0, const char * before_text = "\n", const char * after_text = "\n") const;

        void                    addUnivariateDistribution(std::string keystr, ProbDistShPtr dist, double initial_value);
        void                    addMultivariateDistribution(std::string keystr, MultivarProbDistShPtr dist, double_vect_t const & initial_value);

        TopoProbCalculatorShPtr getTopoProbCalculator();
        void                    addTopologyDistribution(std::string keystr, TopoProbCalculatorShPtr topo_prior_calculator, TreeShPtr tree);

        void                    addTreeLengthDistribution(std::string keystr, TreeLengthDistributionShPtr dist, TreeShPtr tree);

        void                    addInternalEdgelenDistribution(std::string keystr, ProbDistShPtr dist, TreeShPtr tree);
        void                    addExternalEdgelenDistribution(std::string keystr, ProbDistShPtr dist, TreeShPtr tree);

        void                    addEdgelenHyperprior(std::string keystr, ProbDistShPtr dist, TreeShPtr t, double initial_value);

        double                  univariateModified(std::string name, double v);
        double                  univariateModifiedNoDebugCheck(std::string name, double v);

        double                  multivariateModified(std::string name, const double_vect_t & v);

        double                  allEdgeLensModified(TreeShPtr t);
        double                  externalAndInternalEdgeLensModified(std::string name1, std::string name2, TreeShPtr t);
        double                  externalEdgeLensModified(std::string name, TreeShPtr t);
        double                  internalEdgeLensModified(std::string name, TreeShPtr t);

        double                  treeLengthModified(std::string name, TreeShPtr t);

        double                  edgeLenHyperparamModified(std::string name, TreeShPtr t, double v);

        double                  topologyModified(std::string name, TreeShPtr t);

        bool                    isTreeLengthPrior() const;

        bool                    _checkDistrMapKey(std::string name) const;
        std::string             _debugShowDistrMapKeys();

    private:

        void                    _replaceInternalEdgelenPriorMean(double new_mean, TreeShPtr tree);
        void                    _replaceExternalEdgelenPriorMean(double new_mean, TreeShPtr tree);
        void                    _replaceExternalAndInternalEdgelenPriorMean(double new_external_mean, double new_internal_mean, TreeShPtr tree);

        std::string             _internal_hyper_key;    /**< Key used by edgeLenHyperparamModified to determine which type of hyperparameter has been modified; set in constructor */
        std::string             _external_hyper_key;    /**< Key used by edgeLenHyperparamModified to determine which type of hyperparameter has been modified; set in constructor */
        std::string             _universal_hyper_key;   /**< Key used by edgeLenHyperparamModified to determine which type of hyperparameter has been modified; set in constructor */
        std::string             _internal_edgelens_key; /**< Key used by _replaceInternalEdgelenPriorMean to access the internal edge length prior; set in constructor */
        std::string             _external_edgelens_key; /**< Key used by _replaceExternalEdgelenPriorMean to access the external edge length prior; set in constructor */

        double                  _log_joint_prior;       /**< Current value of log(joint prior) */
        prior_distr_vect_t      _distr_vect;            /**< Vector of all prior distributions managed; allows fast traversal in order added */
        prior_distr_map_t       _distr_map;             /**< Map of all prior distributions managed; allows fast lookups for individual parameters */

        mutable std::vector<std::string> _debug_tmp;    /**< Temporary workspace used to enumerate parameters that are out of sync in debugCheckLogJointPrior function */
        mutable bool            _skip_next_debug_check; /**< if set to true, next call to debugCheckLogJointPrior will be a no-op */
	};

typedef boost::shared_ptr<JointPriorManager>    JointPriorManagerShPtr;
typedef boost::weak_ptr<JointPriorManager>		JointPriorManagerWkPtr;

} // namespace phycas

#endif
