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
#include <boost/cast.hpp>       //temporary!
#include <boost/math/special_functions/fpclassify.hpp>
#include "probability_distribution.hpp"
#include "relative_rate_distribution.hpp"
#include "model.hpp"
#include "basic_tree_node.hpp"
#include "tree_likelihood.hpp"
#include "xlikelihood.hpp"
#include "mcmc_chain_manager.hpp"
#include "dirichlet_move.hpp"
#include "basic_tree.hpp"
#include "tree_manip.hpp"
#include "gtr.hpp"
#include "hky.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets tuning parameter `psi' and `max_psi' to 300.0 and calls reset().
*/
DirichletMove::DirichletMove() : MCMCUpdater()
	{
	is_move = true;
	dim = 0;
	psi = 300.0;
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Resets the 'dir_forward' and 'dir_reverse' shared pointers.
*/
void DirichletMove::reset()
	{
    dir_forward.reset();
    dir_reverse.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'psi', which is the tuning parameter for this move.
*/
void DirichletMove::setTuningParameter(
  double x) /* is the new value for psi */
	{
	psi = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'dim', which is the number of parameters updated jointly by this move.
*/
void DirichletMove::setDimension(
  unsigned d) /* is the new value for psi */
	{
	dim = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'dim', which is the number of parameters updated jointly by this move.
*/
unsigned DirichletMove::getDimension() const
	{
	return dim;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'psi', which is the tuning parameter for this move.
*/
double DirichletMove::getTuningParameter() const
	{
	return psi;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Chooses a new vector of state frequencies using a sharp Dirichlet distribution centered at the original frequencies.
*/
void DirichletMove::proposeNewState()
	{
    // copy the current parameters from the model to the data member orig_params
    getParams();

    // create vector of Dirichlet parameters for selecting new parameter values (by multiplying
    // each current parameter value by the value `psi')
    c_forward.resize(orig_params.size());
	std::transform(orig_params.begin(), orig_params.end(), c_forward.begin(), 1.0 + boost::lambda::_1*psi);

	//std::cerr << "c_forward:" << std::endl;
	//std::copy(c_forward.begin(), c_forward.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;

	// create Dirichlet distribution and sample from it to get proposed frequencies
    dir_forward = DirichletShPtr(new DirichletDistribution(c_forward));
    dir_forward->SetLot(getLot().get());
    new_params = dir_forward->Sample();

	//std::cerr << "new_params:" << std::endl;
	//std::copy(new_params.begin(), new_params.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;

    // create vector of Dirichlet parameters for selecting old frequencies (needed for Hasting ratio calculation)
    c_reverse.resize(new_params.size());
	std::transform(new_params.begin(), new_params.end(), c_reverse.begin(), 1.0 + boost::lambda::_1*psi);
    dir_reverse = DirichletShPtr(new DirichletDistribution(c_reverse));
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is:
|>
|	Pr(new state -> old state)   Dir(c'[i]*psi) density evaluated at c
|	-------------------------- = -------------------------------------
|	Pr(old state -> new state)   Dir(c[i]*psi) density evaluated at c'
|>
|   where c[i] is the ith current parameter and c'[i] is the ith proposed new parameter.
*/
double DirichletMove::getLnHastingsRatio() const
	{
	double log_prob_reverse_move = dir_reverse->GetLnPDF(orig_params);
	double log_prob_forward_move = dir_forward->GetLnPDF(new_params);
    double log_hastings_ratio = log_prob_reverse_move - log_prob_forward_move;
	return log_hastings_ratio;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double DirichletMove::getLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted. Simply calls the reset() function.
*/
void DirichletMove::accept()
	{
	MCMCUpdater::accept();
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Calls the reset() function after reinstating the original state frequencies
|   and ensuring that all conditional likelihood arrays will be recalculated when the likelihood is next calculated.
*/
void DirichletMove::revert()
	{
	MCMCUpdater::revert();
    setParams(orig_params);

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();
    jpm->multivariateModified(name, orig_params);
    curr_ln_prior = jpm->getLogJointPrior();

    // invalidate all CLAs
    likelihood->useAsLikelihoodRoot(NULL);
    likelihood->storeAllCLAs(tree);

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or
|	revert(), whichever is appropriate.
*/
bool DirichletMove::update()
	{
	if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

	proposeNewState();

	double prev_ln_like = p->getLastLnLike();
    double prev_ln_prior = jpm->getLogJointPrior();
	PHYCAS_ASSERT(!use_ref_dist || mv_ref_dist);
	double prev_ln_ref_dist = (use_ref_dist ? mv_ref_dist->GetLnPDF(orig_params) : 0.0);

    // replace current parameter values with new ones
    setParams(new_params);
    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs

    jpm->multivariateModified(name, new_params);
    curr_ln_prior = jpm->getLogJointPrior();

	double curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
	double curr_ln_ref_dist = (use_ref_dist ? mv_ref_dist->GetLnPDF(new_params) : 0.0);

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;

	if (is_standard_heating)
		{
		prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
		curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
		if (use_ref_dist)
			{
			prev_posterior += (1.0 - heating_power)*prev_ln_ref_dist;
			curr_posterior += (1.0 - heating_power)*curr_ln_ref_dist;
			}
		}
	else
		{
		prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
		curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
		}

	double ln_hastings			= getLnHastingsRatio();
	double ln_accept_ratio		= curr_posterior - prev_posterior + ln_hastings;

    double lnu = std::log(rng->Uniform());

    bool accepted = false;

	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("ACCEPT, prev_ln_ref_dist = %.5f, curr_ln_ref_dist = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\norig_params: ") % prev_ln_ref_dist % curr_ln_ref_dist % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nnew_params:  ";
			for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nc_forward:  ";
			for (std::vector<double>::const_iterator it = c_forward.begin(); it != c_forward.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\ndir_forward:  ";
			debug_info += dir_forward->GetDistributionDescription();
			}
		p->setLastLnLike(curr_ln_like);

		accept();
		accepted = true;
		}
	else
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("REJECT, prev_ln_ref_dist = %.5f, curr_ln_ref_dist = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\norig_params: ") % prev_ln_ref_dist % curr_ln_ref_dist % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nnew_params:  ";
			for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nc_forward:  ";
			for (std::vector<double>::const_iterator it = c_forward.begin(); it != c_forward.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\ndir_forward:  ";
			debug_info += dir_forward->GetDistributionDescription();
			}
		curr_ln_like	= p->getLastLnLike();

		revert();
		accepted = false;
		}

    double inverse_psi = 1.0/psi;
    inverse_psi = p->adaptUpdater(inverse_psi, nattempts, accepted);
    psi = 1.0/inverse_psi;

    return accepted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns an empty vector. Override this function in derived classes to return a vector
|	of parameter values.
*/
double_vect_t DirichletMove::listCurrValuesFromModel()
	{
	double_vect_t v;
	return v;
	}

