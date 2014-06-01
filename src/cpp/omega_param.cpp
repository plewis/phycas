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
#include "codon_model.hpp"
#include "tree_likelihood.hpp"
#include "mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its Codon pointer to NULL. Also sets the
|	`curr_value' data member to 1.0 and refreshes `curr_ln_prior' accordingly.
*/
OmegaParam::OmegaParam()
  : MCMCUpdater(), codon(NULL)
	{
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool OmegaParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return false;
	slice_sampler->Sample();

	ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        debug_info = str(boost::format("OmegaParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This override of the base class virtual function uses dynamic_cast to set the `codon' data member from the supplied
|	model shared pointer. Having a Codon pointer allows this parameter to call Codon-specific member functions that are
|	not present in the base class Model (such as setOmega).
*/
void OmegaParam::setModel(ModelShPtr m)
	{
	MCMCUpdater::setModel(m);
	Model * p = m.get();
	codon = dynamic_cast<Codon *>(p);	// forces inclusion of "phycas/src/likelihood_models.hpp"

	//POL tried unsuccessfully to get this to compile as an inlined function, but VC gave me this
	// error (which makes sense):
	//
	// ...mcmc_param.inl(43) : error C2680: 'phycas::Codon *' :
	//  invalid target type for dynamic_cast
	//  'Codon' : class must be defined before using in a dynamic_cast
	//
	// Decided to add this comment because otherwise I will forget this and be tempted to move the
	// function back to mcmc_param.inl
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the value of `omega' in `model' to `v'.
*/
void OmegaParam::sendCurrValueToModel(double v)
	{
	codon->setOmega(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the value of `omega' in `model'.
*/
double OmegaParam::getCurrValueFromModel() const
	{
    return codon->getOmega();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	OmegaParam is a functor whose operator() returns a value proportional to the full-conditional posterior probability
|	density for a particular value of omega, the nonsynonymous/synonymous rate ratio. If the supplied omega value `w' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double OmegaParam::operator()(
  double w)	/**< is a new value for the parameter omega */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (w > 0.0)
		{
		PHYCAS_ASSERT(codon);
		sendCurrValueToModel(w);

		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        JointPriorManagerShPtr jpm = p->getJointPriorManager();
        jpm->univariateModified(name, w);
        curr_ln_prior = jpm->getLogJointPrior();

		if (is_standard_heating)
			if (use_ref_dist)
				{
				PHYCAS_ASSERT(ref_dist);
				double curr_ln_ref_dist = ref_dist->GetLnPDF(w);
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
