/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2011 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#include <cmath>
#include <iostream>
#include "model.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 2 and 'root_state' data member to 1.
*/
Irreversible::Irreversible()
  : Model(2), scaling_factor(1.0), scaling_factor_fixed(true), root_present(true)
	{
    state_freq_fixed = true;
    setLossOnly();  // default is a parallel loss model

    state_repr.reserve(2);
	state_repr.push_back("0");
	state_repr.push_back("1");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor_fixed' data member to false.
*/
void Irreversible::freeScalingFactor()
    {
    scaling_factor_fixed = false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor_fixed' data member to true.
*/
void Irreversible::fixScalingFactor()
    {
    scaling_factor_fixed = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
 |	Returns current value of the `scaling_factor' data member.
 */
double Irreversible::getScalingFactor()
    {
    return scaling_factor;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor' data member to the supplied scaling factor `sf'. Assumes `sf' is greater than or equal to zero.
*/
void Irreversible::setScalingFactor(
  double sf)    /**< is the new value of scaling_factor */
    {
    PHYCAS_ASSERT(sf >= 0.0);
    scaling_factor = sf;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `scaling_factor_prior'.
*/
ProbDistShPtr Irreversible::getScalingFactorPrior()
    {
	return scaling_factor_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Irreversible::setScalingFactorPrior(ProbDistShPtr d)
    {
	scaling_factor_prior = d;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `root_present' data member to false and state frequencies to (0.0,1.0).
*/
void Irreversible::setGainOnly()
    {
    // gain only means root state is fixed at 0
    root_present = false;
    std::vector<double> f;
    f.push_back(1.0);   // frequency of state 0
    f.push_back(0.0);   // frequency of state 1
    setStateFreqsUnnorm(f);
    }

/*----------------------------------------------------------------------------------------------------------------------
 |	Constructor sets `num_states' data member to 2 and 'root_state' data member to 1.
 */
void Irreversible::setLossOnly()
{
    // loss only means root state is fixed at 1
    root_present = true;
    std::vector<double> f;
    f.push_back(0.0);   // frequency of state 0
    f.push_back(1.0);   // frequency of state 1
    setStateFreqsUnnorm(f);
}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "IrreversibleLoss", "IrreversibleGain,
|   "IrreversibleLoss+I", "IrreversibleGain+I", "IrreversibleLoss+G", "IrreversibleGain+G", "IrreversibleLoss+G+I", or
|   "IrreversibleGain+G+I".
*/
std::string Irreversible::getModelName() const
	{
	std::string s = "Irreversible";
    if (root_present)
        s += "Loss";
    else
        s += "Gain";
	if (is_pinvar_model)
		s += "+I";
	if (num_gamma_rates > 1)
		s += "+G";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for
|	compatability with MrBayes). The Irreversible model provides additional column labels for the gamma shape parameter (if
|   the number of rates is greater than 1) and the pinvar parameter (if an invariable sites model is being used).
*/
std::string Irreversible::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	return Model::paramHeader();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The only parameters in the Irreversible model are edge lengths, which are not reported by this function, but we must
|	provide an override of the pure virtual base class version.
*/
std::string Irreversible::paramReport(
  unsigned ndecimals,						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
	{
	return Model::paramReport(ndecimals, include_edgelen_hyperparams);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Note that for this model, edgeLength = lambda*t,
|   where lambda is the rate and t the time. Overrides the pure virtual function
|   inherited from the base class Model. For the IrreversibleLoss model, the transition probabilities are:
|>
|	P00 = 1
|	P01 = 0
|	P10 = 1 - exp{-scaling_factor*edgeLength}
|	P11 = exp{-scaling_factor*edgeLength}
|>
|	where the scaling_factor parameter should be a value different than 1 only if edge lengths are generated
|   from characters other than this character. For the IrreversibleGain model, the transition probabilities are:
|>
|	P00 = exp{-scaling_factor*edgeLength}
|	P01 = 1 - exp{-scaling_factor*edgeLength}
|	P10 = 0
|	P11 = 1
|>
*/
void Irreversible::calcPMat(double * * pMat, double edgeLength) const
	{
	// The next two lines fix the "Jockusch" bug (see bug 8 in the BUGS file for details)
    if (edgeLength < 1.e-8)
        edgeLength = 1.e-8; //TreeNode::edgeLenEpsilon;
	const double exp_term = exp(-scaling_factor*edgeLength);
    const double prob_change = 1.0 - exp_term;
	const double prob_nochange = exp_term;
	PHYCAS_ASSERT(prob_change >= 0.0);
    if (root_present)
        {
        // irreversible loss model
        pMat[0][0] = 1.0;
        pMat[0][1] = 0.0;
        pMat[1][0] = prob_change;
        pMat[1][1] = prob_nochange;
        }
    else
        {
        // irreversible gain model
        pMat[0][0] = prob_nochange;
        pMat[0][1] = prob_change;
        pMat[1][0] = 0.0;
        pMat[1][1] = 1.0;
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Irreversible::calcUniformizationLambda() const
	{
    return 1.1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Irreversible::calcLMat(double * * lMat) const
	{
    // The body of this function is totally bogus right now and this model should not be used
    // with uniformization until something better is put in place
    const double lambda = 1.1;
	lMat[0][0] = 0.0;
	lMat[0][1] = 0.0;
	lMat[1][0] = 0.0;
	lMat[1][1] = 0.0;
    return lambda;
	}

double Irreversible::calcUMat(double * * uMat) const
	{
    // The body of this function is totally bogus right now and this model should not be used
    // with uniformization until something better is put in place
    const double lambda = 1.1;
    uMat[0][0] = 0.0;
    uMat[0][1] = 0.0;
    uMat[1][0] = 0.0;
    uMat[1][1] = 0.0;
    return lambda;
	}

