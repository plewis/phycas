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

//#include "phycas/force_include.h"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "mcmc_updater.hpp"
#include "mcmc_param.hpp"
#include "basic_tree.hpp"
#include "mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for the base class of all move and parameter classes used for MCMC. MCMCUpdater-derived objects manage
|	either Metropolis-Hastings proposals (moves) or model quantities updated via slice sampling (parameters) in a
|	Bayesian analysis.
*/
MCMCUpdater::MCMCUpdater()
  :
  use_ref_dist(false),
  nattempts(0.0),
  naccepts(0.0),
  curr_value(0.1),
  curr_ln_prior(0.0),
  curr_ln_like(0.0),
  has_slice_sampler(false),
  is_move(false),
  is_master_param(false),
  is_hyper_param(false),
  is_fixed(false),
  slice_max_units(UINT_MAX),
  heating_power(1.0),
  is_standard_heating(true),
  save_debug_info(false)
	{
	//ln_zero = std::log(std::numeric_limits<double>::denorm_min()); // this doesn't work, lnL can get much lower than the log of dnorm_min!
	ln_zero = -DBL_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Destructor.
*/
MCMCUpdater::~MCMCUpdater()
	{
	//std::cerr << "\n\n>>>>MCMCUpdater dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if one of `prior' or `mv_prior' exists. If any of these exists, then this updater is responsible for
|   calculating some component of the joint prior (i.e. it is a prior steward).
*/
bool MCMCUpdater::isPriorSteward() const
	{
	return (prior || mv_prior);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `nattempts', which stores the number of Metropolis-Hastings proposals attempted
|   since the last call to resetDiagnostics. Assumes this updater is a Metropolis-Hastings move; an assert will trip in
|   debug builds if the `slice_sampler' pointer is not NULL because slice samplers keep their own diagnostics and thus
|   `nattempts' is not used if a slice sampler is being used to perform the updates.
*/
double MCMCUpdater::getNumAttempts() const
	{
	PHYCAS_ASSERT(!slice_sampler);
	return nattempts;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `naccepts', which stores the number of Metropolis-Hastings proposals accepted
|   since the last call to resetDiagnostics. Assumes this updater is a Metropolis-Hastings move; an assert will trip in
|   debug builds if the `slice_sampler' pointer is not NULL because slice samplers keep their own diagnostics and thus
|   `naccepts' is not used if a slice sampler is being used to perform the updates.
*/
double MCMCUpdater::getNumAccepts() const
	{
	PHYCAS_ASSERT(!slice_sampler);
	return naccepts;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets both `naccepts' and `nattempts' to 0.0. Assumes this updater is a Metropolis-Hastings move; an assert will
|   trip in debug builds if the `slice_sampler' pointer is not NULL because slice samplers keep their own diagnostics
|   and thus `nattempts' and `naccepts' are not used if a slice sampler is being used to perform the updates.
*/
void MCMCUpdater::resetDiagnostics()
	{
	PHYCAS_ASSERT(!slice_sampler);
	nattempts = 0.0;
	naccepts = 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_fixed'. While `is_fixed' is true, then the update function always returns
|	immediately without ever modifying the parameter value.
*/
bool MCMCUpdater::isFixed() const
	{
	return is_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fixes the value of this parameter to the current value. Subsequent calls to the update member function will have no
|	effect until freeParameter is called. Sets the value of the data member `is_fixed' to true.
*/
void MCMCUpdater::fixParameter()
	{
	is_fixed = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_fixed' to the default value (false). Subsequent calls to the update member
|	function will update the value of the parameter.
*/
void MCMCUpdater::freeParameter()
	{
	is_fixed = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs an update. For moves this involves calling proposeNewState, computing the posterior density of the new
|	state, and deciding whether to accept or reject the proposed move. For parameters, this involves choosing the next
|	parameter value by slice sampling. This base class version does nothing. In derived classes, update should return
|   true if it changed anything, false if whatever was proposed was rejected and thus nothing was changed.
*/
bool MCMCUpdater::update()
	{
    if (save_debug_info)
        {
        debug_info = "In MCMCUpdater::update(): this base class method does nothing";
        }
    return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the `prior' data member points to something.
*/
bool MCMCUpdater::computesUnivariatePrior() const
	{
	return bool(prior);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the `mv_prior' data member points to something.
*/
bool MCMCUpdater::computesMultivariatePrior() const
	{
	return bool(mv_prior);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This version always returns false; if a derived class implements a tree length prior, it should override this
|   function to return true.
*/
bool MCMCUpdater::computesTreeLengthPrior() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false. Derived classes should override to return true if they compute the topology prior.
*/
bool MCMCUpdater::computesTopologyPrior() const
    {
    return false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the value of the data member `is_move' is false, and vice versa.
*/
bool MCMCUpdater::isParameter() const
	{
	return !is_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_master_param', which is true if this is a master parameter (a parameter
|	object that computes the joint prior density for several model parameters but does not update any of them. This
|	function can be queried if it is important to know whether the data member 'curr_value' is meaningful (it is not
|	in the case of a master parameter).
*/
bool MCMCUpdater::isMasterParameter() const
	{
	return is_master_param;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_hyper_param', which is true if this represents a hyperparameter (a model
|	parameter that is part of the prior specification but does not appear in the likelihood function).
*/
bool	MCMCUpdater::isHyperParameter() const
	{
	return is_hyper_param;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `has_slice_sampler'.
*/
bool	MCMCUpdater::hasSliceSampler() const
	{
	return has_slice_sampler;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_move'.
*/
bool	MCMCUpdater::isMove() const
	{
	return is_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing, just returns 0.0. Derived classes representing MCMC moves should override
|	this to return the log of the Hastings ratio (only necessary if the proposal is not symmetric).
*/
double MCMCUpdater::getLnHastingsRatio() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing, just returns 0.0. Derived classes representing MCMC moves should override
|	this to return the log of the Jacobian (only necessary if the proposal changes the model dimension).
*/
double MCMCUpdater::getLnJacobian() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing MCMC moves should override this to implement the
|	actual proposed move. This function should not do the accept/revert component, however. That should be done in
|	update() using the functions accept() and revert().
*/
void MCMCUpdater::proposeNewState()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing MCMC moves should override this to do anything
|	necessary to accept the state created by the last call to proposeNewState().
*/
void MCMCUpdater::accept()
	{
	nattempts++;
	naccepts++;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing MCMC moves should override this to do anything
|	necessary to revert back to the state that existed before the last call to proposeNewState().
*/
void MCMCUpdater::revert()
	{
	nattempts++;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Gets the value of the `weight' data member used to determine the number of times this updater's update() function
|	should be called in each update cycle.
*/
unsigned MCMCUpdater::getWeight() const
	{
	return weight;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the `weight' data member used to determine the number of times this updater's update() function
|	should be called in each update cycle.
*/
void MCMCUpdater::setWeight(unsigned w)
	{
	weight = w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   In derived classes this function can be used to perform autotuning of Metropolis proposals to match a target
|   acceptance rate.
*/
void MCMCUpdater::autotune()
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the move-specific tuning parameter.
*/
double MCMCUpdater::getTuningParameter() const
    {
    return 0.0;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Sets the move-specific tuning parameter.
*/
void MCMCUpdater::setTuningParameter()
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the move-specific data member used to determine the size of a step when exploring the posterior
|   distribution. This base class version does nothing; override in derived classes to set the tuning parameter
|   appropriately.
*/
void MCMCUpdater::setPosteriorTuningParam(
  double x) /**< is the new value of the tuning parameter */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the move-specific data member used to determine the size of a step when exploring the prior
|   distribution. This base class version does nothing; override in derived classes to set the tuning parameter
|   appropriately.
*/
void MCMCUpdater::setPriorTuningParam(
  double x) /**< is the new value of the tuning parameter */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the boldness of proposals. The meaning of boldness is interpreted individually by derived classes.
*/
void MCMCUpdater::setBoldness(
  double b) /**< is the new boldness value */
	{
	//std::cerr << "Calling default (do-nothing) MCMCUpdater::setBoldness for updater named " << name << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of data member `slice_max_units', which is the maximum number of units used for slice sampling
|	updates. If a slice sampler has been created, calls the setMaxUnits function of `slice_sampler'.
*/
void MCMCUpdater::setMaxUnits(unsigned max_units)
	{
	slice_max_units = max_units;
	if (slice_sampler)
		{
		slice_sampler->SetMaxUnits(slice_max_units);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes should override this to compute the log posterior density of
|	the supplied value conditional on all other parameter values.
*/
double MCMCUpdater::operator()(double)
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This creates a SliceSampler object and assigns it to the `slice_sampler' data member, which is a shared_ptr that
|	points to nothing initially. The SliceSampler constructor takes a boost::shared_ptr<Lot> (which we have available)
|	FuncToSampleWkPtr (this object). Thus, the SliceSampler cannot be created in the constructor because "this" does
|	not yet fully exist, hence the need for a createSliceSampler() member function, which needs to be called at some
|	point after the MCMCUpdater-derived object is created and of course before its `slice_sampler' data member starts
|	being used.
*/
void MCMCUpdater::createSliceSampler()
	{
	PHYCAS_ASSERT(!slice_sampler);	// don't want to do this more than once
#if defined(WEAK_FUNCTOSAMPLE)
	slice_sampler.reset(new SliceSampler(rng, FuncToSampleWkPtr(shared_from_this()))); // forces inclusion of "phycas/src/slice_sampler.hpp"
#else
	slice_sampler.reset(new SliceSampler(rng, shared_from_this())); // forces inclusion of "phycas/src/slice_sampler.hpp"
#endif
	slice_sampler->SetMaxUnits(slice_max_units);
	slice_sampler->SetXValue(curr_value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `slice_sampler' data member, which may or may not point to a SliceSampler object. Some MCMCUpdater-
|	derived classes do not use slice sampling to update their parameter's value, so these have an empty `slice_sampler'
|	member. Also, MCMCUpdater-derived classes that represent moves rather than parameters will generally not use their
|	`slice_sampler' data members.
*/
SliceSamplerShPtr MCMCUpdater::getSliceSampler()
	{
	return slice_sampler;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   If this function is called with the value true, the string data member `debug_info' will be cleared and filled upon
|   each update to reflect the last value of the parameter or move. This is useful for debugging for observing exactly
|   where the MCMC analysis begins to diverge from its behavior prior to a change in the code. Call the public member
|   function getDebugInfo to retrieve the `debug_info' string. Calling this function with the value false stops the
|   saving of debug information, thus increasing the efficiency of the MCMC analysis.
*/
void MCMCUpdater::setSaveDebugInfo(
  bool save_info)       /**> is either true (to begin saving debug info) or false (to stop saving debug info) */
	{
	save_debug_info = save_info;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a copy of `debug_info', which holds information about the last update if `save_debug_info' is true. Call
|   the public member function setSaveDebugInfo to change the state of `save_debug_info'. Note that this function leaves
|   the `debug_info' string in its current state, so subsequent calls before the next update will return the same
|   information.
*/
std::string MCMCUpdater::getDebugInfo() const
	{
	return debug_info;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `prior' data member, which points to a ProbabilityDistribution object. This accessor will primarily be
|	of use to moves, as they may use the master edge length updater's prior distribution to directly compute the
|	relative log-prior over only the edge lengths they modify during an update.
*/
ProbDistShPtr MCMCUpdater::getPriorDist()
	{
	return prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager::finalize function calls this function for each MCMCUpdater it knows about. This provides a way
|	for each updater to call back the MCMCChainManager when it needs the joint prior over all parameters.
*/
void MCMCUpdater::setChainManager(
 ChainManagerWkPtr p)		/**< is a pointer to the MCMCChainManager containing this updater */
	{
	chain_mgr = p;
	if (has_slice_sampler)
		createSliceSampler();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Useful primarily for debugging, this function uses the supplied value `x' rather than `curr_value' to set the value
|   of `curr_ln_prior' for this updater.
*/
double MCMCUpdater::setCurrLnPrior(
  double x) /**< is the value to use, temporarily, to set the prior */
    {
	try
		{
		curr_ln_prior = prior->GetLnPDF(x);
		}
	catch(XProbDist &)
		{
		PHYCAS_ASSERT(0);
		}
	return curr_ln_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of data member `use_ref_dist'.
*/
bool MCMCUpdater::useWorkingPrior() const
	{
	return use_ref_dist;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the data member `use_ref_dist'.
*/
void MCMCUpdater::setUseWorkingPrior(
  bool b)	/**< is the new value for `use_ref_dist' */
	{
	use_ref_dist = b;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Retrieves current value for the parameter being managed by this updater, then returns log of the working prior
|	probability density at that value. If this updater is not a prior steward, simply returns 0.0.
*/
double MCMCUpdater::recalcWorkingPrior() const
	{
	if (isFixed() || !isPriorSteward())
		return 0.0;

	// std::cerr << boost::str(boost::format("calling MCMCUpdater::recalcWorkingPrior for %s") % getName()) << std::endl;

	double lnwp = 0.0;
	if (ref_dist)
		{
		double value = getCurrValueFromModel();
		try
			{
			lnwp = ref_dist->GetLnPDF(value);
			}
		catch(XProbDist &)
			{
			PHYCAS_ASSERT(0);
			}
		//std::cerr << boost::str(boost::format("%.8f <-- %.8f <-- %s <-- %s") % lnwp % value % ref_dist->GetDistributionDescription() % getName()) << std::endl;//temp
		}
	else if (mv_ref_dist)
		{
		double_vect_t values;
		getCurrValuesFromModel(values);
		try
			{
			lnwp = mv_ref_dist->GetLnPDF(values);
			}
		catch(XProbDist &)
			{
			PHYCAS_ASSERT(0);
			}
		//std::cerr << boost::str(boost::format("%.8f <-- {") % lnwp);//temp
		//std::copy(values.begin(), values.end(), std::ostream_iterator<double>(std::cerr, " "));//temp
		//std::cerr << boost::str(boost::format("} <-- %s <-- %s") % mv_ref_dist->GetDistributionDescription() % getName()) << std::endl;//temp
		}
	else if (likelihood->getTreeLengthRefDist())   // TREE_LENGTH_DISTRIBUTION
		{
		try
			{
			lnwp = likelihood->getTreeLengthRefDist()->GetLnPDF(tree);
			}
		catch(XProbDist &)
			{
			PHYCAS_ASSERT(0);
			}
		}
    else
        PHYCAS_ASSERT(0);   // should not get here

	return lnwp;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log likelihood just after this updater's update() member function was last called.
*/
double MCMCUpdater::getLnLike() const
	{
	return curr_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the `name' data member.
*/
const std::string & MCMCUpdater::getName() const
	{
	return name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the string returned by the GetDistributionDescription() function of the `prior' or
|   `mv_prior' data member, if either of these is defined.
*/
std::string MCMCUpdater::getPriorDescr() const
	{
	if (prior)
		return prior->GetDistributionDescription();
	else if (mv_prior)
		return mv_prior->GetDistributionDescription();
	else
		return "(no prior defined)";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `name' data member of this move to the supplied string.
*/
void MCMCUpdater::setName(const std::string & s)
	{
	name = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `rng' data member, which is a shared pointer to the pseudorandom number generator (Lot object).
*/
LotShPtr MCMCUpdater::getLot()
	{
	return rng;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `rng' data member to the supplied Lot shared pointer.
*/
void MCMCUpdater::setLot(LotShPtr r)
	{
	rng = r;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' data member to the supplied Tree shared pointer. Also calls the setTree() member function of the
|	`tree_manipulator' data member.
*/
void MCMCUpdater::setTree(TreeShPtr p)
	{
	tree = p;
	tree_manipulator.setTree(p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `likelihood' data member to the supplied TreeLikelihood shared pointer.
*/
void MCMCUpdater::setTreeLikelihood(TreeLikeShPtr L)
	{
	likelihood = L;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `model' data member to the supplied Model shared pointer.
*/
void MCMCUpdater::setModel(ModelShPtr m)
	{
	model = m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `prior' data member to the supplied ProbabilityDistribution shared pointer.
*/
void MCMCUpdater::setPrior(ProbDistShPtr p)
	{
	prior = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `mv_prior' data member to the supplied MultivariateProbabilityDistribution shared pointer.
*/
void MCMCUpdater::setMultivarPrior(MultivarProbDistShPtr p)
	{
	mv_prior = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `ref_dist' data member to the supplied ProbabilityDistribution shared pointer, but only if `use_ref_dist'
|   is true. If  `use_ref_dist' is false, resets `ref_dist' shared pointer.
*/
void MCMCUpdater::setWorkingPrior(ProbDistShPtr p)
	{
    if (use_ref_dist)
        ref_dist = p;
    else
        ref_dist.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `mv_ref_dist' data member to the supplied MultivariateProbabilityDistribution shared pointer, but only if
|   `use_ref_dist' is true. If  `use_ref_dist' is false, resets `mv_ref_dist' shared pointer.
*/
void MCMCUpdater::setMultivariateWorkingPrior(MultivarProbDistShPtr p)
	{
    if (use_ref_dist)
        mv_ref_dist = p;
    else
        mv_ref_dist.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `topo_prob_calc' data member to the supplied FocalTreeTopoProbCalculatorShPtr shared pointer, but only if
|   `use_ref_dist' is true. If  `use_ref_dist' is false, resets `topo_prob_calc' shared pointer.
*/
void MCMCUpdater::setReferenceDistribution(FocalTreeTopoProbCalculatorShPtr p)
	{
    if (use_ref_dist)
        topo_prob_calc = p;
    else
        topo_prob_calc.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the current value to the supplied value `x'. The value is stored in the data member `curr_value' and passed to
|   the slice sampler if and when it is created in the member function MCMCUpdater::createSliceSampler().
*/
void MCMCUpdater::setStartingValue(double x)
	{
    curr_value = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Resets all shared pointers. It is only necessary to reset `model' and `likelihood' in order to avoid problems with
|	shared pointer cycles preventing deletion of objects (in particular, Model->MCMCUpdater->TreeLikelihood->Model);
|	however, we might as well delete all of them. This function is called by MCMCChainManager::releaseUpdaters().
*/
void MCMCUpdater::releaseSharedPointers()
	{
	// model and likelihood must be reset to avoid memory leaks due to shared pointer cycles
	model.reset();
	//likelihood.reset();

	// these are optional
	tree.reset();
	rng.reset();
	prior.reset();
	mv_prior.reset();
	slice_sampler.reset();
	}


/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing model parameters should override this
|	method to set the value of the corresponding parameter in `model' to the supplied value `v'.
*/
void MCMCUpdater::sendCurrValueToModel(double v)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version returns a bogus value. Derived classes representing model parameters should override this
|	method to return the current value of the corresponding parameter in `model'.
*/
double MCMCUpdater::getCurrValueFromModel() const
	{
	return DBL_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing model parameters should override this
|	method to set the values of the corresponding multivariate parameter in `model' to the supplied vector `v'.
*/
void MCMCUpdater::sendCurrValuesToModel(const double_vect_t & v)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply clears the supplied vector `v'. Derived classes representing model parameters should
|	override this method to return the current vector of parameter values from `model'.
*/
void MCMCUpdater::getCurrValuesFromModel(double_vect_t & v) const
	{
	v.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns an empty vector. Derived classes representing model parameters should
|	override this method to return the current vector of parameter values from `model'.
*/
double_vect_t MCMCUpdater::listCurrValuesFromModel()
	{
	double_vect_t v;
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls `likelihood'->calcLnL(`tree') to recalculate `curr_ln_like'. This function is important because the
|	tempting getLnLike() member function only returns the value of `curr_ln_like' (it does not recalculate anything).
*/
double MCMCUpdater::recalcLike()
	{
	likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	curr_ln_like = likelihood->calcLnL(tree);

	//std::cerr << "--> in MCMCUpdater::recalcLike(), name = " << getName();	//POL temp
	//std::cerr << ", curr_ln_like = " << curr_ln_like;	//POL temp
	//std::cerr << std::endl;	//POL temp

	return curr_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log prior just after this updater's update() member function was last called.
*/
double MCMCUpdater::getLnPrior() const
	{
	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `tree' data member.
*/
TreeShPtr MCMCUpdater::getTree()
	{
	return tree;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `likelihood' data member.
*/
TreeLikeShPtr MCMCUpdater::getTreeLikelihood()
	{
	return likelihood;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the power used in heating (variable `heating_power') to the specified value `p'.
*/
void MCMCUpdater::setPower(
  double p)
	{
	heating_power = p;
    if (slice_sampler)
        {
        if (heating_power > 0.0)
            slice_sampler->UseDoublingMethod(false);
        else
            slice_sampler->UseDoublingMethod(true);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the power used in heating (variable `heating_power').
*/
double MCMCUpdater::getPower() const
	{
	return heating_power;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if either the data member `ref_dist' or the data member `mv_ref_dist' points to an object.
*/
bool MCMCUpdater::isWorkingPrior() const
	{
	return (ref_dist || mv_ref_dist);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_standard_heating' to true.
*/
void MCMCUpdater::setStandardHeating()
	{
	is_standard_heating = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `is_standard_heating' to false.
*/
void MCMCUpdater::setLikelihoodHeating()
	{
	is_standard_heating = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the data member `is_standard_heating' is true.
*/
bool MCMCUpdater::isStandardHeating() const
	{
	return is_standard_heating;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the data member `is_standard_heating' is false.
*/
bool MCMCUpdater::isLikelihoodHeating() const
	{
	return !is_standard_heating;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if the data member `heating_power' is set to exactly 1.0.
*/
bool MCMCUpdater::isNoHeating() const
	{
	return (heating_power == 1.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log of the working prior density at the supplied value `x'. Assumes `ref_dist' (or
|	`mv_ref_dist') exists, whichever is appropriate.
*/
double MCMCUpdater::calcLnWorkingPriorPDF(
  double x) const	/**< is the value at which the log density is to be computed */
	{
	PHYCAS_ASSERT(prior && ref_dist);	// if this function is called, updater should be a prior steward and have a ref_dist data member that actually points to something
	return ref_dist->GetLnPDF(x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `ref_dist' data member.
*/
double MCMCUpdater::sampleWorkingPrior() const
	{
	PHYCAS_ASSERT(ref_dist);
	return ref_dist->Sample();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `mv_ref_dist' data member.
*/
double_vect_t MCMCUpdater::sampleMultivariateWorkingPrior() const
	{
	PHYCAS_ASSERT(mv_ref_dist);
	return mv_ref_dist->Sample();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the string returned by the GetDistributionDescription() function of the
|	`ref_dist' or `mv_ref_dist' data member, whichever is appropriate.
*/
std::string MCMCUpdater::getWorkingPriorDescr() const
	{
	if (prior && ref_dist)
		return ref_dist->GetDistributionDescription();
	else if (mv_prior && mv_ref_dist)
		return mv_ref_dist->GetDistributionDescription();
	else
		return "(no working prior defined)";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize a new BetaDistribution, which is then stored in `ref_dist'.
|	Assumes `fitting_sample' has more than 1 element.
*/
void MCMCUpdater::fitBetaWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(fitting_sample.size() > 1);
		double n = (double)fitting_sample.size();
		double sum = 0.0;
		double sum_of_squares = 0.0;
		for (double_vect_t::iterator i = fitting_sample.begin(); i != fitting_sample.end(); ++i)
			{
			double v = (*i);
			sum += v;
			sum_of_squares += v*v;
			}
		double mean = sum/n;
		double variance = (sum_of_squares - n*mean*mean)/(n - 1.0);

		// Let a, b be the parameters of a Beta(a,b) and let phi = a + b
		// Note that:
		//     mean = a/phi
		//   1-mean = b/phi
		// variance = a*b/[phi^2*(phi + 1)]
		// Letting z = mean*(1-mean)/variance,
		// phi can be estimated as z - 1
		// Now, a = mean*phi and b = (1-mean)*phi
		double phi = mean*(1.0-mean)/variance - 1.0;
		double a = phi*mean;
		double b = phi*(1.0 - mean);
		ref_dist = ProbDistShPtr(new BetaDistribution(a, b));
		// std::cerr << boost::str(boost::format("@@@@@@@@@ working prior is Beta(%g, %g) for updater %s: mean = %g, variance = %g") % a % b % getName() % mean % variance) << std::endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize a new GammaDistribution, which is then stored in `ref_dist'.
|	Assumes `fitting_sample' has more than 1 element.
*/
void MCMCUpdater::fitGammaWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(fitting_sample.size() > 1);
		double n = (double)fitting_sample.size();
		double sum = 0.0;
		double sum_of_squares = 0.0;
		for (double_vect_t::iterator i = fitting_sample.begin(); i != fitting_sample.end(); ++i)
			{
			double v = (*i);
			sum += v;
			sum_of_squares += v*v;
			}
		double mean = sum/n;	// shape*scale
		double variance = (sum_of_squares - n*mean*mean)/(n - 1.0);	// shape*scale^2
		double scale = variance/mean;
		PHYCAS_ASSERT(scale > 0.0);
		double shape = mean/scale;
		// std::cerr << boost::str(boost::format("@@@@@@@@@ working prior is Gamma(%g,%g) for updater %s: mean = %g, variance = %g") % shape % scale % getName() % mean % variance) << std::endl;
		ref_dist = ProbDistShPtr(new GammaDistribution(shape, scale));
        ref_dist->SetLot(rng.get());
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize a new LognormalDistribution, which is then stored in `ref_dist'.
|	Assumes `fitting_sample' has more than 1 element.
*/
void MCMCUpdater::fitLognormalWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(fitting_sample.size() > 1);
		double n = (double)fitting_sample.size();
		double sum = 0.0;
		double sum_of_squares = 0.0;
		for (double_vect_t::iterator i = fitting_sample.begin(); i != fitting_sample.end(); ++i)
			{
			double v = (*i);
			double logv = log(v);
			sum += logv;
			sum_of_squares += logv*logv;
			}
		double logmean = sum/n;
		double logvar = (sum_of_squares - n*logmean*logmean)/(n - 1.0);
		double logsd = sqrt(logvar);
		ref_dist = ProbDistShPtr(new LognormalDistribution(logmean, logsd));
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls recalcRelativeRates() function of the `likelihood' data member.
*/
void MCMCUpdater::recalcRelativeRates()
	{
    likelihood->recalcRelativeRates();
    }

}	// namespace phycas
