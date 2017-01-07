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
#include "boost/cast.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 0.5 (or 2.0
|	if `invert' is true.
*/
DiscreteGammaShapeParam::DiscreteGammaShapeParam()
  : MCMCUpdater()
	{
    curr_value = 0.5;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
DiscreteGammaShapeParam::~DiscreteGammaShapeParam()
	{
 	//std::cerr << "\n>>>>> DiscreteGammaShapeParam dying...likelihood use count = " << likelihood.use_count() << std::endl;
 	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool DiscreteGammaShapeParam::update()
	{
	if (is_fixed)
		{
		return false;
		}

	slice_sampler->Sample();

    //POLTMP
    //std::vector<double> v = slice_sampler->DebugSample();
    //if (fabs(slice_sampler->GetLastSampledXValue() - 0.00034951) < 1.e-8)
    //    {
    //    // 0: sampled x
    //    // 1: x-coord of vertical slice
    //    // 2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
    //    // 3: x-coord of left edge of horizontal slice
    //    // 4: x-coord of right edge of horizontal slice
    //    // 5: y-coord of horizontal slice
    //    // 6: horizontal slice interval width
    //    // 7+: x-coord of failed sampling attempts (y-coord all equal to element 5)
    //    std::cerr << "********** curr_value = " << slice_sampler->GetLastSampledXValue() << std::endl;
    //    std::cerr << "********** " << v[0] << " <-- sampled x" << std::endl;
    //    std::cerr << "********** " << v[1] << " <-- x-coord of vertical slice" << std::endl;
    //    std::cerr << "********** " << v[2] << " <-- y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)" << std::endl;
    //    std::cerr << "********** " << v[3] << " <-- x-coord of left edge of horizontal slice" << std::endl;
    //    std::cerr << "********** " << v[4] << " <-- x-coord of right edge of horizontal slice" << std::endl;
    //    std::cerr << "********** " << v[5] << " <-- y-coord of horizontal slice" << std::endl;
    //    std::cerr << "********** " << v[6] << " <-- horizontal slice interval width" << std::endl;
    //    for (unsigned i = 7; i < v.size(); ++i)
    //        std::cerr << "********** " << v[i] << " <-- x-coord of failed sampling attempts (y-coord all equal to element 5)" << std::endl;
    //    std::exit(0);
    //    }

	ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        debug_info = str(boost::format("DiscreteGammaShapeParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the value of `gamma_shape' in `model' to `v'.
*/
void DiscreteGammaShapeParam::sendCurrValueToModel(double v)
	{
    model->setShape(v);	// change the gamma shape parameter in the model
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the value of `gamma_shape' in `model'.
*/
double DiscreteGammaShapeParam::getCurrValueFromModel() const
	{
	double a = model->getShape();
    return a;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	DiscreteGammaShapeParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of the gamma shape parameter. If the supplied gamma shape value `a' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double DiscreteGammaShapeParam::operator()(
  double a)	/**< is a new value for the gamma shape parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (a > 0.0)
		{
		sendCurrValueToModel(a);
		likelihood->recalcRelativeRates();	// must do this whenever model's shape parameter changes
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);
        JointPriorManagerShPtr jpm = p->getJointPriorManager();
        jpm->univariateModified(name, a);
        curr_ln_prior = jpm->getLogJointPrior();

		if (is_standard_heating)
			if (use_ref_dist)
				{
				PHYCAS_ASSERT(ref_dist);
				double curr_ln_ref_dist = ref_dist->GetLnPDF(a);
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
