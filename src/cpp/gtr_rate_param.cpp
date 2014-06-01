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
#include "gtr.hpp"
#include "tree_likelihood.hpp"
#include "mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its GTR pointer to NULL. Also sets the
|	`curr_value' data member to 1.0 and refreshes `curr_ln_prior' accordingly. Assumes `w' is greater than or equal to
|	zero and less than 6.
*/
GTRRateParam::GTRRateParam(
  unsigned w)		/**< The 0-based index of the relative rate being managed by this object (0=AC, 1=AG, 2=AT, 3=CG, 4=CT and 5=GT) */
  : MCMCUpdater(), gtr(NULL), which(w)
	{
	PHYCAS_ASSERT(w >= 0);
	PHYCAS_ASSERT(w < 6);
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool GTRRateParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return false;
	slice_sampler->Sample();

	ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        debug_info = str(boost::format("GTRRateParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This override of the base class virtual function uses dynamic_cast to set the `gtr' data member from the supplied
|	model shared pointer. Having a GTR pointer allows this parameter to call GTR-specific member functions not present
|	in the base class Model (such as setRelRates).
*/
void GTRRateParam::setModel(ModelShPtr m)
	{
	MCMCUpdater::setModel(m);
	Model * p = m.get();
	gtr = dynamic_cast<GTR *>(p);	// forces inclusion of "likelihood_models.hpp"

	// If tempted to move this to mcmc_param.inl, see comment in KappaParam::setModel function
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the appropriate relative rate in `model' to the value `v'.
*/
void GTRRateParam::sendCurrValueToModel(double v)
	{
	gtr->setRelRateUnnorm(which, v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the appropriate relative rate value in `model'.
*/
double GTRRateParam::getCurrValueFromModel() const
	{
    return gtr->getRelRateUnnorm(which);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	GTRRateParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of one of the six relative rates. If the supplied rate value `r' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double GTRRateParam::operator()(
  double r)	/**< is a new value for the relative rate parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

    if (r > 0.0)
		{
		PHYCAS_ASSERT(which < 6);
		PHYCAS_ASSERT(gtr);
		sendCurrValueToModel(r);
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        JointPriorManagerShPtr jpm = p->getJointPriorManager();
        jpm->univariateModified(name, r);
        curr_ln_prior = jpm->getLogJointPrior();

		if (is_standard_heating)
			if (use_ref_dist)
				{
				PHYCAS_ASSERT(ref_dist);
				double curr_ln_ref_dist = ref_dist->GetLnPDF(r);
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
