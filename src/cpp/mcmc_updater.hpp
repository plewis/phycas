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

#if ! defined(MCMC_UPDATER_HPP)
#define MCMC_UPDATER_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "tree_manip.hpp"				// for TreeManip
#include "slice_sampler.hpp"
#include "probability_distribution.hpp"
#include "multivariate_probability_distribution.hpp"
#include "relative_rate_distribution.hpp"
#include "lognormal.hpp"
#include "states_patterns.hpp"		// for double_vect
#include "topo_prior_calculator.hpp"
#include "tree_length_distribution.hpp"	// for TreeLengthDistributionShPtr

namespace phycas
{

class Tree;
typedef boost::shared_ptr<Tree>					TreeShPtr;

class TreeLikelihood;
typedef boost::shared_ptr<TreeLikelihood>		TreeLikeShPtr;

class Model;
typedef boost::shared_ptr<Model>				ModelShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>		ChainManagerWkPtr;

class MCMCUpdater;
typedef boost::shared_ptr<MCMCUpdater>			MCMCUpdaterShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Class encapsulating both MCMC moves (i.e. Metropolis-Hastings updates) and model parameters (slice sampling
|	updates). Provides member functions for specifying a name to be used for this updater, as well as objects needed to
|	compute the poseterior density: a tree, a substitution model, a pseudo-random number generator, and a likelihood
|	calculator. These should all be set before adding the MCMCUpdater to the MCMCChainManager. Inside
|	MCMCChainManager::finalize(), the MCMCChainManager calls the setChainManager() itself for all MCMCUpdater objects
|	that have been added to it. The key virtual member function of this interface is update(). MCMC move objects
|	override update() to perform a Metropolis-Hastings proposal and either accept or reject it. Parameter objects
|	override update() to perform a single slice sampling update of the model parameter they are managing. Below is a
|	list of some of the differences between moves and parameters:
|
|	Moves usually have these properties:
|	1) do NOT use the createSliceSampler() or getSliceSampler() member functions and thus leave their `slice_sampler'
|		data member empty
|	2) do NOT use the `curr_value' data member or the getCurrValue() member function
|	4) do NOT use the operator() member function
|	5) DO use the getLnHastingsRatio(), getLnJacobian(), proposeNewState(), accept() and revert() member functions
|	6) some use the `tree_manipulator' data member to effect topological rearrangements
|	7) the `has_slice_sampler' data member is usually false (exception is TreeScalerMove)
|	8) the `is_move' data member is true
|
|	Parameters usually have these properties:
|	1) do NOT use the createSliceSampler() or getSliceSampler() member functions and thus leave their `slice_sampler'
|		data member empty
|	2) DO use the `curr_value' data member or the getCurrValue() member function
|	4) DO usethe operator() member function
|	5) do NOT use the getLnHastingsRatio(), getLnJacobian(), proposeNewState(), accept() and revert() member functions
|	6) never use the `tree_manipulator' data member
|	7) the `has_slice_sampler' data member is usually true (exception is EdgeLenMasterParam)
|	8) the `is_move' data member is false
|
|	The `weight' data member determines how many times in a single update cycle the update() function of this
|	particular MCMCUpdater will be called. Typically, parameters have weight 1 whereas moves (e.g. LargetSimonMove)
|	have a higher weight so that their proposals will be attempted several times for every parameter update.
*/
class MCMCUpdater : public AdHocDensity, public boost::enable_shared_from_this<MCMCUpdater>
	{
	public:
                                    MCMCUpdater();
                                    virtual	~MCMCUpdater();

		// Predicates
		bool                        isFixed() const;
		virtual bool                isPriorSteward() const;
		bool                        useWorkingPrior() const;
		bool                        isParameter() const;
		bool                        isMasterParameter() const;
		bool                        isHyperParameter() const;
		bool                        hasSliceSampler() const;
		bool                        isMove() const;
		bool                        computesUnivariatePrior() const;
		bool                        computesMultivariatePrior() const;
		virtual bool                computesTopologyPrior() const;
        virtual bool                computesTreeLengthPrior() const;

		// Accessors
		const std::string &         getName() const;
		unsigned                    getWeight() const;
		double                      getLnLike() const;
		double                      getLnPrior() const;
		virtual std::string         getPriorDescr() const;
		virtual std::string         getWorkingPriorDescr() const;
		virtual double              sampleWorkingPrior() const;
		double_vect_t               sampleMultivariateWorkingPrior() const;
		LotShPtr                    getLot();
        std::string                 getDebugInfo() const;
        void                        setPower(double p);
        double                      getPower() const;
        void                        setStandardHeating();
        void                        setLikelihoodHeating();
        bool                        isStandardHeating() const;
        bool                        isLikelihoodHeating() const;
        bool                        isNoHeating() const;
        double                      setCurrLnPrior(double x);

        double                      getTuningParameter() const;
        void                        setTuningParameter();
        void                        autotune();

        TreeLikeShPtr               getTreeLikelihood();
        TreeShPtr                   getTree();

		// Accessors used only by parameters
		SliceSamplerShPtr           getSliceSampler();

		// Accessors used only by moves
		ProbDistShPtr               getPriorDist();

		// Modifiers
		virtual void                setName(const std::string & s);
		virtual void                setWeight(unsigned w);
		virtual void                setBoldness(double b);
        virtual void                setPosteriorTuningParam(double x);
        virtual void                setPriorTuningParam(double x);
        void                        setMaxUnits(unsigned max_units);
		virtual void                setLot(LotShPtr p);
		virtual void                setTree(TreeShPtr p);
		virtual void                setTreeLikelihood(TreeLikeShPtr p);
		virtual void                setModel(ModelShPtr p);
		virtual void                setChainManager(ChainManagerWkPtr p);
        void                        setSaveDebugInfo(bool save_info);
		virtual void                setPrior(ProbDistShPtr p);
		virtual void                setMultivarPrior(MultivarProbDistShPtr p);
		virtual void                setWorkingPrior(ProbDistShPtr p);
		virtual void                setMultivariateWorkingPrior(MultivarProbDistShPtr p);
		virtual void                setReferenceDistribution(FocalTreeTopoProbCalculatorShPtr p);
		void                        setUseWorkingPrior(bool b);
		virtual double              recalcWorkingPrior() const;

		// Modifiers used only by parameters
		virtual void                setStartingValue(double x);

		// Univariate case
        virtual void                sendCurrValueToModel(double v);
        virtual double              getCurrValueFromModel() const;

		// Multivariate case
        virtual void                sendCurrValuesToModel(const double_vect_t & v);
        virtual void                getCurrValuesFromModel(double_vect_t & v) const;
        virtual double_vect_t       listCurrValuesFromModel();

		// Utilities
		void                        releaseSharedPointers();
		virtual bool                update();
		virtual double              recalcLike();
		void                        recalcRelativeRates();

		// Utilites for reporting diagnostics
		double                      getNumAttempts() const;
		double                      getNumAccepts() const;
		void                        resetDiagnostics();

		// Utilities used only by parameters
		void                        fixParameter();
		void                        freeParameter();
		void                        createSliceSampler();
		virtual double              operator()(double);

		// Utilities related to working priors used in steppingstone sampling
		bool                        isWorkingPrior() const;
		double                      calcLnWorkingPriorPDF(double x) const;
		void                        fitBetaWorkingPrior();
		void                        fitGammaWorkingPrior();
		void                        fitLognormalWorkingPrior();

		// Note: some member functions could be made pure virtuals were it not for a bug in the
		// boost::lambda library that causes compiles to fail if attempting to use boost::lambda::bind
		// with a shared_ptr to an abstract base class. Even if the shared_ptr points to a derived class,
		// the boost::lambda library attempts (without success) to create an object of the abstract base
		// class. We have two choices: give up using boost::lambda::bind (painful) or avoid abstract base
		// classes that are going to be used with shared_ptr (also painful, but less so). Hopefully, the
		// bug will eventually get fixed and we can have our cake and eat it too. For a discussion of
		// this bug by the author of the boost::lambda library (Jaakko Jarvi), see
		// http://lists.boost.org/Archives/boost/2003/10/55139.php

	protected:

		// Virtual functions needed only by moves
		//
		virtual double              getLnHastingsRatio() const;
		virtual double              getLnJacobian() const;
		virtual void                proposeNewState();
		virtual void                accept();
		virtual void                revert();

	protected:

		std::string                 name;					/**< Storage for the name of this parameter to be used in reporting to the user */
		unsigned                    weight;					/**< The number of times this updater's update() function should be called in each update cycle */
		TreeShPtr                   tree;					/**< The tree on which the likelihood will be calculated */
		TreeManip                   tree_manipulator;		/**< The object that facilitates topological rearrangements */
		ModelShPtr                  model;					/**< The substitution model to be used in computing the likelihood */
		TreeLikeShPtr               likelihood;				/**< The object that knows how to calculate the likelihood */
		LotShPtr                    rng;					/**< The pseudorandom number generator object used in updating parameter value */

        // one (and only one) of these will be set
		ProbDistShPtr               prior;					/**< The probability distribution serving as the prior for a univariate parameter */
        MultivarProbDistShPtr       mv_prior;               /**< The probability distribution serving as the prior for a multivariate parameter */

		bool                        use_ref_dist;           /**< If true, the (presumably already parameterized) working prior will be used in computing the power posterior in steppingstone sampling */
		ProbDistShPtr               ref_dist;               /**< The probability distribution serving as the working prior for a univariate parameter */
        MultivarProbDistShPtr       mv_ref_dist;            /**< The probability distribution serving as the working prior for a multivariate parameter */
		double_vect_t               fitting_sample;			/**< Storage for sample used in fitting a univariate working prior during steppingstone sampling */
		double_vect_vect_t          mv_fitting_sample;		/**< Storage for sample used in fitting a multivariate working prior during steppingstone sampling */
		SliceSamplerShPtr           slice_sampler;			/**< The slice sampler used by parameters for updating (not used by moves) */
		ChainManagerWkPtr           chain_mgr;				/**< The object that knows how to compute the joint log prior density */
		double                      nattempts;				/**< The number of update attempts made since last call to resetDiagnostics (used only by Metropolis-Hastings moves, slice samplers maintain their own diagnostics) */
		double                      naccepts;				/**< The number of times a proposed move was accepted since last call to resetDiagnostics (used only by Metropolis-Hastings moves, slice samplers maintain their own diagnostics) */
		double                      curr_value;				/**< The current value of this parameter */
		double                      curr_ln_prior;			/**< The value of log prior after most recent call to this object's update() */
		double                      curr_ln_like;			/**< The value of log likelihood after most recent call to this object's update() */
		double                      ln_zero;				/**< The value to return for the log posterior density if the parameter value is out of range */
		bool                        has_slice_sampler;		/**< If true, then createSliceSampler() will be called inside setChainManager, ensuring that all parameters that need a slice sampler have one before they begin updating */
		bool                        is_move;				/**< True if this updater is a move (i.e. encapsulates a Metropolis-Hastings proposal); if false, then this updater is considered a parameter */
		bool                        is_master_param;		/**< True if this updater is a master parameter (does not update any model parameters but can compute the joint prior density for several model parameters) */
		bool                        is_hyper_param;			/**< True if this updater represents a hyperparameter (a model parameter that is part of the prior specification but not the likelihood function) */
		bool                        is_fixed;				/**< If true, update returns immediately so parameter is never updated */
		unsigned                    slice_max_units;		/**< Maximum number of units used by `slice_sampler' */
        std::string                 debug_info;				/**< Information about the last update, only created if save_debug_info is true */
        double                      heating_power;          /**< The power to which the posterior (in standard heating) or just the likelihood (in likelihood heating) is raised. To not heat, specify 1.0. */
        bool                        is_standard_heating;    /**< If true, standard heating is used (posterior is raised to the power `heating_power'); otherwise, likelihood heating is used (just the likelihood is raised to the power `heating_power') */
		bool                        save_debug_info;		/**< If true, information about the last update will be stored in debug_info string */

        FocalTreeTopoProbCalculatorShPtr topo_prob_calc;
	};

typedef std::vector<MCMCUpdaterShPtr>		MCMCUpdaterVect;
typedef MCMCUpdaterVect::iterator			MCMCUpdaterIter;
typedef MCMCUpdaterVect::const_iterator		MCMCUpdaterConstIter;

} // namespace phycas

#endif
