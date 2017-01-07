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

#if ! defined(PARTITION_MODEL_HPP)
#define PARTITION_MODEL_HPP

#include "subset_proportions.hpp"
#include "model.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Container for information about a data partition: stores number of subsets, a vector of the number of rate
|   categories for each subset, the number of states for each subset, and a vector of shared pointers to the model
|   applied to each subset.
*/
class PartitionModel
{
	friend class TreeLikelihood;
	friend class TipData;
	friend class InternalData;

public:

        								PartitionModel();
    virtual								~PartitionModel();

	// Utilities
    void                                createJointPriorManager();
    unsigned                        	getNumSubsets() const;
    unsigned                        	getTotalNumPatterns() const;
    unsigned                        	getTotalNumSites() const;
    unsigned                            getNumFreeParameters() const;
    std::vector<std::string>            getPWKParameterNames() const;
    std::vector<std::string>            getFreeParameterNames() const;
    std::vector<std::string>            getAllParameterNames() const;
    std::vector<double>                 getUntransformedParameters() const;
    std::vector<double>                 getTransformedParameters() const;
    void                                setTransformedParameters(const std::vector<double> v, TreeShPtr tree);
    double                              getLogDetJacobian() const;
    double                              calcLogDetJacobian() const;
    void                                appendUntransformedSubsetRelativeRates(std::vector<double> & values) const;
    void                                appendTransformedSubsetRelativeRates(std::vector<double> & values) const;

    void                                debugShowTransformedRelativeRates(const std::string msg) const;

	// Accessors
	ModelShPtr							getModel(unsigned i) {return subset_model[i];}
	const ModelVect					&	getModelsVect() const;
	double								getSubsetRelRate(unsigned i) const;
	const std::vector<double>		&	getSubsetRelRatesVect() const;
	MultivarProbDistShPtr				getSubsetRelRatePrior() const;
	const std::vector<unsigned>		&	getNumPatternsVect() const;
	unsigned							getNumPatterns(unsigned i) const {return subset_num_patterns[i];}
	const std::vector<unsigned>		&	getNumStatesVect() const;
	unsigned							getNumStates(unsigned i) const {return subset_num_states[i];}
	const std::vector<unsigned>		&	getNumRatesVect() const;
	unsigned							getNumRates(unsigned i) const {return subset_num_rates[i];}
	unsigned	 						getNumSites(unsigned i) const;
	const std::vector<unsigned>		&	getNumSitesVect() const;
    const std::vector<unsigned>     &   getSiteAssignments() const;
    SubsetProportionsShPtr              getSubsetProportions();
    JointPriorManagerShPtr              getJointPriorManager();

	// Setters
	void								addModel(ModelShPtr m);
	void								setModelsVect(const std::vector<ModelShPtr> & models);
	void								setSubsetRelRatesVect(const std::vector<double> & rrates);
	void								setSubsetRelRatePrior(MultivarProbDistShPtr rrate_prior);
	void								setNumPatternsVect(const std::vector<unsigned> & npatterns);
	void								setNumSitesVect(const std::vector<unsigned> & nsites);
	void								setNumStatesVect(const std::vector<unsigned> & nstates);
	void								setNumRatesVect(const std::vector<unsigned> & nrates);
    void                                setSiteAssignments(const std::vector<unsigned> & v);

private:

	MultivarProbDistShPtr				subset_relrate_prior;	/**< The prior on subset relative rates */
    std::vector<double>           		subset_relrates;		/**< The relative rates for each partition subset */
    std::vector<unsigned>           	subset_num_patterns;	/**< The number of site patterns for each partition subset */
    std::vector<unsigned>           	subset_num_sites;		/**< The number of sites in each partition subset */
    std::vector<unsigned>           	subset_num_states;		/**< The number of states for each partition subset */
    std::vector<unsigned>           	subset_num_rates;		/**< The number of relative rate categories for each partition subset */
    std::vector<unsigned>           	site_assignments;		/**< The index of the model assigned to each site */
    ModelVect							subset_model;			/**< The substitution model for each partition subset */
    mutable std::vector<std::string>    free_param_names;       /**< After a call to getFreeParameterNames(), holds the names of free parameters in all subset models */
    mutable std::vector<std::string>    param_names;            /**< After a call to getAllParameterNames(), holds the names of all parameters in all subset models */
    mutable std::vector<double>         free_param_values;      /**< After a call to getTransformedParameters(), holds the transformed values of all free parameters in all subset models */
    mutable std::vector<double>         param_values;           /**< After a call to getUntransformedParameters(), holds the untransformed values of all parameters in all subset models */
    SubsetProportionsShPtr              _subset_proportions;    /**< Holds the vector of subset proportions (_subset_proportions[i] is the proportion of sites in subset i) */
    JointPriorManagerShPtr              _joint_prior_manager;   /**< Maintains the value of the prior density for the current state of the chain */
};  // class PartitionModel

typedef boost::shared_ptr<PartitionModel> PartitionModelShPtr;

} // namespace phycas

#endif
