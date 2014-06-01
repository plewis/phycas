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

#include "model.hpp"
#include "tree_likelihood.hpp"
#include "mcmc_param.hpp"
#include "mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The StateFreqParam constructor requires the caller to specify a value for `which'. `which' is 0, 1, 2, or 3 for
|	nucleotide models and  determines the particular state frequency (e.g. A, C, G, or T, respectively) this object
|	represents. It calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 1.0 and
|	refreshes `curr_ln_prior' accordingly.
*/
StateFreqParam::StateFreqParam(
  unsigned w)		/**< The 0-based index of the base frequency being managed by this object (0=A, 1=C, 2=G and 3=T) */
  : MCMCUpdater(),  which(w)
	{
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
StateFreqParam::~StateFreqParam()
	{
	//std::cerr << "\n>>>>> StateFreqParam dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool StateFreqParam::update()
	{
	if (is_fixed)
		return false;

	slice_sampler->Sample();

	ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        debug_info = str(boost::format("StateFreqParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the appropriate unnormalized state frequency in `model' to `v'.
*/
void StateFreqParam::sendCurrValueToModel(double v)
	{
	model->setStateFreqUnnorm(which, v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the appropriate unnormalized state frequency in `model'.
*/
double StateFreqParam::getCurrValueFromModel() const
	{
    return model->getStateFreqUnnorm(which);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	StateFreqParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of f, the quantity governing a particular relative base frequency (but
|	f is not a relative base frequency itself, as it is not constrained to lie between 0.0 and 1.0, nor is it
|	constrained by any of the analogous quantities governing the other three base frequencies. If the supplied base
|	frequency parameter value `f' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to
|	a log posterior equal to negative infinity).
*/
double StateFreqParam::operator()(
  double f) /**< is a new value for the base frequency parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

    if (f > 0.0)
		{
		PHYCAS_ASSERT((model->isCodonModel() && which < 61) || which < 4);

		sendCurrValueToModel(f);
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        JointPriorManagerShPtr jpm = p->getJointPriorManager();
        jpm->univariateModified(name, f);
        curr_ln_prior = jpm->getLogJointPrior();

		if (is_standard_heating)
			if (use_ref_dist)
				{
				PHYCAS_ASSERT(ref_dist);
				double curr_ln_ref_dist = ref_dist->GetLnPDF(f);
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

