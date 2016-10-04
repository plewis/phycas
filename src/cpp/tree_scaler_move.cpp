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

#include "joint_prior_manager.hpp"
#include "probability_distribution.hpp"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "xlikelihood.hpp"
#include "mcmc_chain_manager.hpp"
#include "tree_scaler_move.hpp"
#include "basic_tree.hpp"

#include "boost/format.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor.
*/
TreeScalerMove::TreeScalerMove() : MCMCUpdater()
	{
	is_move = true;
    n = 0;
    m = 0.0;
    mstar = 0.0;
    forward_scaler = 1.0;
    reverse_scaler = 1.0;
    lambda = 0.5;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter used for exploring the posterior
|   distribution in this move.
*/
void TreeScalerMove::setTuningParameter(
  double x) /* is the new value for `lambda' */
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is (`mstar'/`m')^n, where `mstar' is
|	the new tree length and `m' is the tree length before the move is proposed. The value n is the number of edges in
|   the tree.
*/
double TreeScalerMove::getLnHastingsRatio() const
	{
    double ln_hastings = (double)n*(std::log(mstar) - std::log(m));
    return ln_hastings;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Proposes a new value by which to scale the tree length.
*/
void TreeScalerMove::proposeNewState()
	{
    n = tree->GetNNodes() - 1;  // root node's edge does not count
    m = tree->EdgeLenSum();
    mstar = m*exp(lambda*(rng->Uniform() - 0.5));
    forward_scaler = mstar/m;
    reverse_scaler = m/mstar;
    tree->ScaleAllEdgeLens(forward_scaler);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool TreeScalerMove::update()
	{
    if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

    bool using_tree_length_prior = jpm->isTreeLengthPrior();

    // If first edge length master param returns false for useWorkingPrior(), then we are using
    // Mark Holder's variable tree topology reference distribution (see _provideRefDistToUpdaters in MCMCImpl.py)
    MCMCUpdaterVect::const_iterator it = p->getEdgeLenParams().begin();
    bool using_vartopol_prior = !(*it)->useWorkingPrior();

	double prev_ln_prior = jpm->getLogJointPrior();

	double prev_ln_like = p->getLastLnLike();
	double prev_ln_ref_dist	= (use_ref_dist ? recalcWorkingPriorForMove(using_tree_length_prior, using_vartopol_prior) : 0.0);

    proposeNewState();

    if (jpm->isTreeLengthPrior())
        jpm->treeLengthModified("tree_length", tree);
    else
        jpm->allEdgeLensModified(tree);
	curr_ln_prior = jpm->getLogJointPrior();

    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	double curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
	double curr_ln_ref_dist	= (use_ref_dist ? recalcWorkingPriorForMove(using_tree_length_prior, using_vartopol_prior) : 0.0);

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
			debug_info = boost::str(boost::format("ACCEPT, forward_scaler = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f") % forward_scaler % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			}
		p->setLastLnLike(curr_ln_like);
		accept();
		accepted = true;
		}
	else
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("REJECT, forward_scaler = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f") % forward_scaler % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			}
		curr_ln_like = p->getLastLnLike();
		revert();
		accepted = false;
		}

    lambda = p->adaptUpdater(lambda, nattempts, accepted);

    return accepted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in proposeNewState.
*/
void TreeScalerMove::revert()
	{
	MCMCUpdater::revert();
    tree->ScaleAllEdgeLens(reverse_scaler);

    likelihood->useAsLikelihoodRoot(NULL);
    likelihood->storeAllCLAs(tree);	// force CLAs to be recalculated

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();
    jpm->allEdgeLensModified(tree);
    curr_ln_prior = jpm->getLogJointPrior();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void TreeScalerMove::accept()
	{
	MCMCUpdater::accept();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function overrides the base class version to always returns true because this derived class implements a tree
|   length prior.
*/
bool TreeScalerMove::computesTreeLengthPrior() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint log working prior over all edges in the associated tree.
*/
double TreeScalerMove::recalcWorkingPriorForMove(
  bool using_tree_length_prior,         /*< true if using the Ranalla-Yang tree length prior */
  bool using_vartopol_prior) const     /*< true if using the Holder et al. variable topology reference distribution */
	{
    double ln_ref_dist = 0.0;
    if (using_tree_length_prior)
        {
        ln_ref_dist = likelihood->getTreeLengthRefDist()->GetLnPDF(tree);
        }
    else if (using_vartopol_prior)
        {
        // Computes the log of the probability of the tree under Mark Holder's variable tree topology reference distribution
        PHYCAS_ASSERT(topo_prob_calc);
        std::pair<double, double> treeprobs = topo_prob_calc->CalcTopologyLnProb(*tree, true);
        const double ln_ref_topo = treeprobs.first;
        const double ln_ref_edges = treeprobs.second;
        ln_ref_dist = ln_ref_topo + ln_ref_edges;
        }
    else
        {
        // Loop through all EdgeLenMasterParam objects and call the recalcWorkingPrior function of each.
        // Each EdgeLenMasterParam object knows how to compute the working prior for the edge lengths it controls.
        ChainManagerShPtr p = chain_mgr.lock();
        const MCMCUpdaterVect & edge_length_params = p->getEdgeLenParams();
        for (MCMCUpdaterVect::const_iterator it = edge_length_params.begin(); it != edge_length_params.end(); ++it)
            {
            if (!(*it)->isFixed())
                ln_ref_dist += (*it)->recalcWorkingPrior();
            }
        }

	return ln_ref_dist;
	}

}	// namespace phycas
