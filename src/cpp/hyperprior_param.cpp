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

#include "mcmc_param.hpp"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor of HyperPriorParam simply calls the base class (MCMCUpdater) constructor. Like other parameters,
|   sets `is_move' to false and `has_slice_sampler' to true. Because this is a hyperparameter, sets `is_hyper_param' to
|   true. Also sets `is_master_param' to false and `curr_value' to 0.1.
*/
HyperPriorParam::HyperPriorParam()
  : MCMCUpdater()
	{
    //@POL this default constructor should never be used (but needs to be defined to avoid errors related to
    // python/C++ container conversions)
    std::cerr << "*** fatal error: default constructor used for HyperPriorParam ***" << std::endl;
    std::exit(0);
	edge_type = all;
	curr_value = 0.1;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor of HyperPriorParam simply calls the base class (MCMCUpdater) constructor. Like other parameters,
|   sets `is_move' to false and `has_slice_sampler' to true. Because this is a hyperparameter, sets `is_hyper_param' to
|   true. Also sets `is_master_param' to false, `curr_value' to 0.1 and `edgelen_master_param' to the supplied edge
|   length master parameter shared pointer `p'.
*/
HyperPriorParam::HyperPriorParam(
  HyperPriorParam::EdgeTypeEnum t)		/**< determines which type of edges (internal, external, all) are managed by this hyperparameter */
  : MCMCUpdater(), edge_type(t)
	{
	curr_value = 0.1;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
HyperPriorParam::~HyperPriorParam()
	{
	//std::cerr << "\n>>>>> HyperPriorParam dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool HyperPriorParam::update()
	{
	if (is_fixed)
		return false;
	if (slice_sampler)
		{
		slice_sampler->Sample();
		}

	//ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        if (slice_sampler)
            debug_info = str(boost::format("HyperPriorParam %f") % (slice_sampler->GetLastSampledXValue()));
        else
            debug_info = "HyperPriorParam (no slice sampler)";
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set value of hyperparameter in `model' to `v'.
*/
void HyperPriorParam::sendCurrValueToModel(double v)
	{
	if (edge_type == all || edge_type == external)
		model->setExternalEdgelenHyperparam(v);
	else
		model->setInternalEdgelenHyperparam(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return value of hyperparameter from model.
*/
double HyperPriorParam::getCurrValueFromModel() const
	{
	if (edge_type == all || edge_type == external)
		return model->getExternalEdgelenHyperparam();
	else
		return model->getInternalEdgelenHyperparam();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	HyperPriorParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of mu, the mean of the edge length prior distribution. If the supplied
|	value `mu' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior
|	equal to negative infinity).
*/
double HyperPriorParam::operator()(
  double mu) 	/**< is a new parameter value for the prior affecting all edge lengths */
	{
	//double edgeLensLnPrior = 0.0;
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (mu > 0.0)
		{
		sendCurrValueToModel(mu);

		// no need to invalidate all CLAs because only the prior is changing
		//@POL the likelihood is not even needed here because mu does not appear in the likelihood function anywhere
		curr_ln_like = likelihood->calcLnL(tree);

		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
        JointPriorManagerShPtr jpm = p->getJointPriorManager();
        jpm->edgeLenHyperparamModified(name, tree, mu);
        curr_ln_prior = jpm->getLogJointPrior();

        // Store current log-likelihood value (if next updater is a move, it will then not have to calculate it)
		p->setLastLnLike(curr_ln_like);

		if (is_standard_heating)
			if (use_ref_dist)
				{
				PHYCAS_ASSERT(ref_dist);
				double curr_ln_ref_dist = ref_dist->GetLnPDF(mu);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_ref_dist;
				}
			else
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

}
