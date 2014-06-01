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
Binary::Binary()
  : Model(2),scaling_factor(1.0),kappa(1.0)
	{
    state_freq_fixed = false;

    state_freqs[0] = 0.5;
    state_freqs[1] = 0.5;

    state_repr.reserve(2);
	state_repr.push_back("0");
	state_repr.push_back("1");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor_fixed' data member to false.
*/
void Binary::freeScalingFactor()
    {
    scaling_factor_fixed = false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor_fixed' data member to true.
*/
void Binary::fixScalingFactor()
    {
    scaling_factor_fixed = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `scaling_factor' data member.
*/
double Binary::getScalingFactor()
    {
    return scaling_factor;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor' data member to the supplied scaling factor `sf'. Assumes `sf' is greater than or equal to zero.
*/
void Binary::setScalingFactor(
  double sf)    /**< is the new value of scaling_factor */
    {
    PHYCAS_ASSERT(sf >= 0.0);
    scaling_factor = sf;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `scaling_factor_prior'.
*/
ProbDistShPtr Binary::getScalingFactorPrior()
    {
	return scaling_factor_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `scaling_factor_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Binary::setScalingFactorPrior(ProbDistShPtr d)
    {
	scaling_factor_prior = d;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_fixed' data member to false.
*/
void Binary::freeKappa()
    {
    kappa_fixed = false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_fixed' data member to true.
*/
void Binary::fixKappa()
    {
    kappa_fixed = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `kappa' data member.
*/
double Binary::getKappa()
    {
    return kappa;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' data member to the supplied forward/reverse rate ratio `k'. Assumes `k' is greater than or equal to zero.
*/
void Binary::setKappa(
  double k)    /**< is the new value of kappa */
    {
    PHYCAS_ASSERT(k >= 0.0);
    kappa = k;
    state_freqs[0] = 1.0/(1.0 + kappa);
    state_freqs[1] = 1.0 - state_freqs[0];

}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member 'kappa_prior'.
*/
ProbDistShPtr Binary::getKappaPrior()
    {
	return kappa_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Binary::setKappaPrior(ProbDistShPtr d)
    {
	kappa_prior = d;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "Binary", "Binary+I", "Binary+G", or "Binary+I+G".
*/
std::string Binary::getModelName() const
	{
	std::string s = "Binary";
	if (is_pinvar_model)
		s += "+I";
	if (num_gamma_rates > 1)
		s += "+G";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for
|	compatability with MrBayes). The Binary model provides additional column labels for the gamma shape parameter (if
|   the number of rates is greater than 1) and the pinvar parameter (if an invariable sites model is being used).
*/
std::string Binary::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	return Model::paramHeader();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes).
*/
std::string Binary::paramReport(
  unsigned ndecimals,						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
    {
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s = "";
    s += str(boost::format(fmt) % scaling_factor);
    s += str(boost::format(fmt) % kappa);
    s += str(boost::format(fmt) % state_freqs[0]);
	s += str(boost::format(fmt) % state_freqs[1]);
	s += Model::paramReport(ndecimals, include_edgelen_hyperparams);
	return s;
    }

/*----------------------------------------------------------------------------------------------------------------------
 |	Computes the transition probability matrix given an edge length v. Overrides the pure virtual function inherited from
 |	the base class Model. For the Binary model, the transition probabilities are:
 |>
 |	P00 = [1/(1+k)] + [k/(1+k)] exp{-v(1+k)^2/(2k)}
 |	P01 = [k/(1+k)] - [k/(1+k)] exp{-v(1+k)^2/(2k)}
 |	P10 = [1/(1+k)] - [1/(1+k)] exp{-v(1+k)^2/(2k)}
 |	P11 = [k/(1+k)] + [k/(1+k)] exp{-v(1+k)^2/(2k)}
 |>
 */
void Binary::calcPMat(double * * pMat, double edgeLength) const
{
	// The next two lines fix the "Jockusch" bug (see bug 8 in the BUGS file for details)
    if (edgeLength < 1.e-8)
        edgeLength = 1.e-8; //TreeNode::edgeLenEpsilon;
	const double exp_term = exp(-edgeLength*scaling_factor*pow(1.0 + kappa,2.0)/(2.0*kappa));
    const double pi0 = state_freqs[0];
    const double pi1 = state_freqs[1];
	pMat[0][0] = pi0 + pi1*exp_term;
	pMat[0][1] = pi1 - pi1*exp_term;
	pMat[1][0] = pi0 - pi0*exp_term;
	pMat[1][1] = pi1 + pi0*exp_term;
}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Binary::calcUniformizationLambda() const
	{
    return 1.1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Binary::calcLMat(double * * lMat) const
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

double Binary::calcUMat(double * * uMat) const
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

