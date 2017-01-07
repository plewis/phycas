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

//#include <cmath>
//#include <iostream>
#include "gtr.hpp"
//#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
//#	include <boost/python/numeric.hpp>
//#	include "phycas/src/thirdparty/num_util/num_util.h"
//#endif
//using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 4, the base frequencies `piA', `piC', `piG' and `piT' to 0.25, and all
|	six relative rates to 1.0.
*/
GTR::GTR()
  : Model(4), rel_rates_fixed(false)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");

	rel_rates.assign(6, 1.0);

	q_matrix.setRelativeRates(rel_rates);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "GTR", "GTR+G", "GTR+I" or "GTR+G+I".
*/
std::string GTR::getModelName() const
	{
	std::string s = "GTR";
	if (num_gamma_rates > 1)
		s += "+G";
	if (is_pinvar_model)
		s += "+I";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Model::createParameters to create the edge length parameters, the edge length hyperparameter, and any
|	parameters related to rate heterogeneity. This function then adds additional GTR-specific parameters to the
|	supplied `parameters' vector. This incudes the four base frequencies and six relative rate parameters.
*/
void	GTR::createParameters(
  JointPriorManagerShPtr jpm,               /**< is the object that keeps the joint prior density up-to-date */
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters,				/**< is the vector of model-specific parameters to fill */
  int subset_pos) 							/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added; otherwise, the `edgelens' and `edgelen_hyperparams' vectors returned will be empty */
	{
	Model::createParameters(jpm, t, edgelens, edgelen_hyperparams, parameters, subset_pos);

    subset_index = subset_pos < 0 ? 0 : subset_pos;
    std::string nm;

	PHYCAS_ASSERT(rel_rate_params.empty());
    PHYCAS_ASSERT(rel_rate_param_prior || rel_rate_prior);
	relrate_name.clear();
	freq_name.clear();
    //if (subset_pos < 0)
    //    {
    //    relrate_name.push_back("rAC");
    //    relrate_name.push_back("rAG");
    //    relrate_name.push_back("rAT");
    //    relrate_name.push_back("rCG");
    //    relrate_name.push_back("rCT");
    //    relrate_name.push_back("rGT");
    //    freq_name.push_back("freqA");
    //    freq_name.push_back("freqC");
    //    freq_name.push_back("freqG");
    //    freq_name.push_back("freqT");
    //    }
    //else
    //    {
    //    unsigned m = subset_pos + 1;
    //    relrate_name.push_back(boost::str(boost::format("rAC_%d") % m));
    //    relrate_name.push_back(boost::str(boost::format("rAG_%d") % m));
    //    relrate_name.push_back(boost::str(boost::format("rAT_%d") % m));
    //    relrate_name.push_back(boost::str(boost::format("rCG_%d") % m));
    //    relrate_name.push_back(boost::str(boost::format("rCT_%d") % m));
    //    relrate_name.push_back(boost::str(boost::format("rGT_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqA_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqC_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqG_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqT_%d") % m));
    //    }
    unsigned m = subset_index + 1;
    relrate_name.push_back(boost::str(boost::format("%d_rAC") % m));
    relrate_name.push_back(boost::str(boost::format("%d_rAG") % m));
    relrate_name.push_back(boost::str(boost::format("%d_rAT") % m));
    relrate_name.push_back(boost::str(boost::format("%d_rCG") % m));
    relrate_name.push_back(boost::str(boost::format("%d_rCT") % m));
    relrate_name.push_back(boost::str(boost::format("%d_rGT") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqA") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqC") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqG") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqT") % m));

    double_vect_t starting_relrates(6, 1.0/12.0);
    starting_relrates[1] = 4.0/12.0;
    starting_relrates[4] = 4.0/12.0;
    if (rel_rate_param_prior)
        {
        // The fact that rel_rate_param_prior points to something means that the user
        // chose to update the relative rates separately using slice sampling. If the
        // user had chosen to update the relative rates jointly (using the RelRateMove
        // Metropolis-Hastings move), then rel_rate_prior would be set and rel_rate_param_prior
        // would be empty.

        MCMCUpdaterShPtr rAC_param = MCMCUpdaterShPtr(new GTRRateParam(0));
		rAC_param->setName(relrate_name[0]);
        rAC_param->setTree(t);
        rAC_param->setStartingValue(starting_relrates[0]);
        rAC_param->setPrior(rel_rate_param_prior);
        nm = boost::str(boost::format(relrate_name[0]));
        jpm->addUnivariateDistribution(nm, rel_rate_param_prior, starting_relrates[0]);
        if (rel_rates_fixed)
            rAC_param->fixParameter();
        parameters.push_back(rAC_param);
        rel_rate_params.push_back(rAC_param);

        MCMCUpdaterShPtr rAG_param = MCMCUpdaterShPtr(new GTRRateParam(1));
		rAG_param->setName(relrate_name[1]);
        rAG_param->setTree(t);
        rAG_param->setStartingValue(starting_relrates[1]);
        rAG_param->setPrior(rel_rate_param_prior);
        nm = boost::str(boost::format(relrate_name[1]));
        jpm->addUnivariateDistribution(nm, rel_rate_param_prior, starting_relrates[1]);
        if (rel_rates_fixed)
            rAG_param->fixParameter();
        parameters.push_back(rAG_param);
        rel_rate_params.push_back(rAG_param);

        MCMCUpdaterShPtr rAT_param = MCMCUpdaterShPtr(new GTRRateParam(2));
		rAT_param->setName(relrate_name[2]);
        rAT_param->setTree(t);
        rAT_param->setStartingValue(starting_relrates[2]);
        rAT_param->setPrior(rel_rate_param_prior);
        nm = boost::str(boost::format(relrate_name[2]));
        jpm->addUnivariateDistribution(nm, rel_rate_param_prior, starting_relrates[2]);
        if (rel_rates_fixed)
            rAT_param->fixParameter();
        parameters.push_back(rAT_param);
        rel_rate_params.push_back(rAT_param);

        MCMCUpdaterShPtr rCG_param = MCMCUpdaterShPtr(new GTRRateParam(3));
		rCG_param->setName(relrate_name[3]);
        rCG_param->setTree(t);
        rCG_param->setStartingValue(starting_relrates[3]);
        rCG_param->setPrior(rel_rate_param_prior);
        nm = boost::str(boost::format(relrate_name[3]));
        jpm->addUnivariateDistribution(nm, rel_rate_param_prior, starting_relrates[3]);
        if (rel_rates_fixed)
            rCG_param->fixParameter();
        parameters.push_back(rCG_param);
        rel_rate_params.push_back(rCG_param);

        MCMCUpdaterShPtr rCT_param = MCMCUpdaterShPtr(new GTRRateParam(4));
		rCT_param->setName(relrate_name[4]);
        rCT_param->setTree(t);
        rCT_param->setStartingValue(starting_relrates[4]);
        rCT_param->setPrior(rel_rate_param_prior);
        nm = boost::str(boost::format(relrate_name[4]));
        jpm->addUnivariateDistribution(nm, rel_rate_param_prior, starting_relrates[4]);
        if (rel_rates_fixed)
            rCT_param->fixParameter();
        parameters.push_back(rCT_param);
        rel_rate_params.push_back(rCT_param);

        MCMCUpdaterShPtr rGT_param = MCMCUpdaterShPtr(new GTRRateParam(5));
		rGT_param->setName(relrate_name[5]);
        rGT_param->setStartingValue(starting_relrates[5]);
        rGT_param->setTree(t);
        rGT_param->setPrior(rel_rate_param_prior);
        nm = boost::str(boost::format(relrate_name[5]));
        jpm->addUnivariateDistribution(nm, rel_rate_param_prior, starting_relrates[5]);
        if (rel_rates_fixed)
            rGT_param->fixParameter();
        parameters.push_back(rGT_param);
        rel_rate_params.push_back(rGT_param);
        }
    else
        {
        nm = boost::str(boost::format("%d_relrates") % (subset_index + 1));
        jpm->addMultivariateDistribution(nm, rel_rate_prior, starting_relrates);
        }

	PHYCAS_ASSERT(freq_params.empty());
    PHYCAS_ASSERT(freq_param_prior || freq_prior);
    double_vect_t starting_freqs(4, 0.25);
    if (freq_param_prior)
        {
        // Only add frequency parameters if freqs will be updated separately
        // The other option is to update the frequencies jointly using the
        // StateFreqMove Metropolis-Hastings move (in which case freq_prior
        // will be set and freq_param_prior will be empty)

	    MCMCUpdaterShPtr freqA_param = MCMCUpdaterShPtr(new StateFreqParam(0));
		freqA_param->setName(freq_name[0]);
	    freqA_param->setTree(t);
	    freqA_param->setStartingValue(starting_freqs[0]);
	    freqA_param->setPrior(freq_param_prior);
        nm = boost::str(boost::format(freq_name[0]));
        jpm->addUnivariateDistribution(nm, freq_param_prior, starting_freqs[0]);
	    if (state_freq_fixed)
		    freqA_param->fixParameter();
	    parameters.push_back(freqA_param);
	    freq_params.push_back(freqA_param);

	    MCMCUpdaterShPtr freqC_param = MCMCUpdaterShPtr(new StateFreqParam(1));
		freqC_param->setName(freq_name[1]);
	    freqC_param->setTree(t);
	    freqC_param->setStartingValue(starting_freqs[1]);
	    freqC_param->setPrior(freq_param_prior);
        nm = boost::str(boost::format(freq_name[1]));
        jpm->addUnivariateDistribution(nm, freq_param_prior, starting_freqs[1]);
	    if (state_freq_fixed)
		    freqC_param->fixParameter();
	    parameters.push_back(freqC_param);
	    freq_params.push_back(freqC_param);

	    MCMCUpdaterShPtr freqG_param = MCMCUpdaterShPtr(new StateFreqParam(2));
		freqG_param->setName(freq_name[2]);
	    freqG_param->setTree(t);
	    freqG_param->setStartingValue(starting_freqs[2]);
	    freqG_param->setPrior(freq_param_prior);
        nm = boost::str(boost::format(freq_name[2]));
        jpm->addUnivariateDistribution(nm, freq_param_prior, starting_freqs[2]);
	    if (state_freq_fixed)
		    freqG_param->fixParameter();
	    parameters.push_back(freqG_param);
	    freq_params.push_back(freqG_param);

	    MCMCUpdaterShPtr freqT_param = MCMCUpdaterShPtr(new StateFreqParam(3));
		freqT_param->setName(freq_name[3]);
	    freqT_param->setTree(t);
	    freqT_param->setStartingValue(starting_freqs[3]);
	    freqT_param->setPrior(freq_param_prior);
        nm = boost::str(boost::format(freq_name[3]));
        jpm->addUnivariateDistribution(nm, freq_param_prior, starting_freqs[3]);
	    if (state_freq_fixed)
		    freqT_param->fixParameter();
	    parameters.push_back(freqT_param);
	    freq_params.push_back(freqT_param);
        }
    else
        {
        nm = boost::str(boost::format("%d_state_freqs") % (subset_index + 1));
        jpm->addMultivariateDistribution(nm, freq_prior, starting_freqs);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from
|	the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
void GTR::calcPMat(double * * pMat, double edgeLength) const
	{
	q_matrix.recalcPMat(pMat, edgeLength);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double GTR::calcUniformizationLambda() const
	{
    PHYCAS_ASSERT(0);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double GTR::calcLMat(double * * lMat) const
	{
    std::cerr << "Error in GTR::calcUMat: q_matrix does not yet have the required recalcLMat function" << std::endl;
    PHYCAS_ASSERT(0);
	//q_matrix.recalcLMat(lMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double GTR::calcUMat(double * * uMat) const
	{
    std::cerr << "Error in GTR::calcLMat: q_matrix does not yet have the required recalcUMat function" << std::endl;
    PHYCAS_ASSERT(0);
	//q_matrix.recalcUMat(uMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `rel_rates_fixed' to true. The fixParameter member function of all GTRRateParam objects is
|	either called immediately (if the `rel_rate_params' vector is not empty) or is called in createParameters (when
|	`rel_rate_params' is built).
*/
void GTR::fixRelRates()
	{
	rel_rates_fixed = true;
	if (!rel_rate_params.empty())
		{
#if 1
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = rel_rate_params.begin(); it != rel_rate_params.end(); ++it)
			(*it)->fixParameter();
#else
		std::for_each(rel_rate_params.begin(), rel_rate_params.end(), boost::lambda::bind(&MCMCUpdater::fixParameter, *boost::lambda::_1));
#endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `rel_rates_fixed' to false. The freeParameter member function of all GTRRateParam objects is
|	called immediately if the `rel_rate_params' vector is not empty.
*/
void GTR::freeRelRates()
	{
	rel_rates_fixed = false;
	if (!rel_rate_params.empty())
		{
#if 1
		//@POL good place to use std::for_each with a boost::lambda functor?
		for (MCMCUpdaterVect::iterator it = rel_rate_params.begin(); it != rel_rate_params.end(); ++it)
			(*it)->freeParameter();
#else
		std::for_each(rel_rate_params.begin(), rel_rate_params.end(), boost::lambda::bind(&MCMCUpdater::freeParameter, *boost::lambda::_1));
#endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a copy of the `rel_rates' vector.
*/
std::vector<double> GTR::getRelRates()
 	{
	return rel_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies the supplied vector `rates' into the data member vector `rel_rates'. Throws XLikelihood exception if any of
|	the values of `rates' is less than 0.0, or if the number of elements in `rates' is not equal to 6. Assumes that
|	`rel_rates' has 6 elements.
*/
void GTR::setRelRates(const std::vector<double> & rates)
 	{
	++time_stamp;
	PHYCAS_ASSERT(rel_rates.size() == 6);

	// Ensure that we are not about to set any negative relative rates
	if (std::find_if(rates.begin(), rates.end(), boost::lambda::_1 < 0.0) != rates.end())
		throw XLikelihood("GTR relative rates cannot be less than 0.0");

	// Ensure that user supplied the correct number of relative rates
	if (rates.size() != 6)
		throw XLikelihood(str(boost::format("The number of GTR relative rates is 6, but %d were supplied") % rates.size()));

	// The six relative rates are always normalized to sum to 1.0
    double sum_rates = std::accumulate(rates.begin(), rates.end(), 0.0);
    std::transform(rates.begin(), rates.end(), rel_rates.begin(), boost::lambda::_1/sum_rates);

	// Make sure QMatrix knows about this change in relative rates
	q_matrix.setRelativeRates(rel_rates);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the state frequency ratios. There are `num_states'-1 state frequency ratios:
|   freqC/freqA, freqG/freqA, and freqT/freqA. The normalized values of all `num_states' state frequencies can be
|   determined from these `num_states'-1 ratios. Throws an XLikelihood exception if `which' is not less than
|   `num_states'-1 or if `value' is negative.
*/
void GTR::setStateFreqRatio(
  unsigned ratio_index,		/**< the 0-based index into the `state_freq_ratio' vector of the element to modify */
  double value)				/**< the new value of `state_freq_ratio'[`ratio_index'] */
	{
	Model::setStateFreqRatio(ratio_index, value);
    if (ratio_index == 2)
        q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Set all of the GTR relative rates using the rate ratios stored in the `rel_rate_ratios' vector. There are 5 relative
|   rate ratios: rAG/rAC, rAT/rAC, rCG/rAC, rCT/rAC and rGT/rAC. The normalized values of all six relative rates is
|   determined from these five ratios, and the resulting relative rates copied into the `rel_rates' vector. Assumes that
|   the `rel_rate_ratios' vector is completely specified.
*/
void GTR::calcRelRatesFromRatios()
	{
	PHYCAS_ASSERT(rel_rate_ratios.size() == 5);

    //std::cerr << "In GTR::calcRelRatesFromRatios:" << std::endl;
    //std::cerr << "  rAG/rAC = " << rel_rate_ratios[0] << std::endl;
    //std::cerr << "  rAT/rAC = " << rel_rate_ratios[1] << std::endl;
    //std::cerr << "  rCG/rAC = " << rel_rate_ratios[2] << std::endl;
    //std::cerr << "  rCT/rAC = " << rel_rate_ratios[3] << std::endl;
    //std::cerr << "  rGT/rAC = " << rel_rate_ratios[4] << std::endl;

    double sum_ratios = std::accumulate(rel_rate_ratios.begin(), rel_rate_ratios.end(), 0.0);
    //std::cerr << "  sum_ratios = " << sum_ratios << std::endl;

    std::copy(rel_rate_ratios.begin(), rel_rate_ratios.end(), rel_rates.begin()+1);
    rel_rates[0] = 1.0;
    double rAC = 1.0/(1.0 + sum_ratios);
    std::transform(rel_rates.begin(), rel_rates.end(), rel_rates.begin(), boost::lambda::_1*rAC);

    //std::cerr << "  rAC = " << rel_rates[0] << std::endl;
    //std::cerr << "  rAG = " << rel_rates[1] << std::endl;
    //std::cerr << "  rAT = " << rel_rates[2] << std::endl;
    //std::cerr << "  rCG = " << rel_rates[3] << std::endl;
    //std::cerr << "  rCT = " << rel_rates[4] << std::endl;
    //std::cerr << "  rGT = " << rel_rates[5] << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the GTR relative rate ratios. There are 5 relative rate ratios: rAG/rAC, rAT/rAC,
|   rCG/rAC, rCT/rAC and rGT/rAC. The normalized values of all six relative rates can be determined from these five
|   ratios. Throws an XLikelihood exception if `which' is not less than 5 or if `value' is negative. Assumes that
|   `rel_rate_ratios' has length 5.
*/
void GTR::setRelRateRatio(
  unsigned which,		/**< the 0-based index into the `rel_rate_ratios' vector of the element to modify */
  double value)         /**< the new value of `rel_rate_ratios'[`which'] */
	{
	++time_stamp;
    if (rel_rate_ratios.size() != 5)
        {
        // if rel_rate_ratios not yet allocated, allocate it and fill with negative values
        // no more negative values indicates that all have been specified and relative rates can be calculated
        rel_rate_ratios.assign(5, -1.0);
        }

	// Ensure that user supplied a correct index
	if (which >= 5)
		throw XLikelihood(str(boost::format("Index of a GTR rate ratio parameter must be less than 5, but %d was supplied") % which));

	// Ensure that we are not about to set a negative relative rate
	if (value < 0.0)
		throw XLikelihood(str(boost::format("GTR rate ratio parameters must be greater than or equal to 0.0, but the value %f was supplied") % value));

	// Ok, seems safe to overwrite value now
	rel_rate_ratios[which] = value;

    if (rel_rate_ratios[0] > 0.0 && rel_rate_ratios[1] > 0.0 && rel_rate_ratios[2] > 0.0 && rel_rate_ratios[3] > 0.0 && rel_rate_ratios[4] > 0.0)
        {
        // all rate ratios have been set, so go ahead and compute relative rates from the ratios and
        // inform q_matrix of the new values, then reset all values to -1.0
        calcRelRatesFromRatios();
        q_matrix.setRelativeRates(rel_rates);
        if (_joint_prior_manager)
            _joint_prior_manager->multivariateModified(boost::str(boost::format("%d_relrates") % (subset_index + 1)), rel_rates);
        rel_rate_ratios.assign(5, -1.0);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the GTR relative rate parameters. Throws an XLikelihood exception if
|	`param_index' is not less than 6 or if `value' is negative. Assumes that `rel_rates' has length 6.
*/
void GTR::setRelRateUnnorm(
  unsigned param_index,		/**< the 0-based index into the `rel_rates' vector of the element to modify */
  double value)				/**< the new value of `rel_rates'[`param_index'] */
	{
	++time_stamp;
	PHYCAS_ASSERT(rel_rates.size() == 6);

	// Ensure that user supplied a correct index
	if (param_index >= 6)
		throw XLikelihood(str(boost::format("Index of a GTR rate parameter must be less than 6, but %d was supplied") % param_index));

	// Ensure that we are not about to set a negative relative rate
	if (value < 0.0)
		throw XLikelihood(str(boost::format("GTR rate parameters must be greater than or equal to 0.0, but the value %f was supplied") % value));

	// Ok, seems safe to overwrite value now
	rel_rates[param_index] = value;	//@POL do we need rel_rates if q_matrix is maintaining a copy?

	// Make sure QMatrix knows about this change in relative rates
	q_matrix.setRelativeRates(rel_rates);	//@POL could speed up by just changing one rate rather than copying all 6 each time
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `rel_rate_prior'.
*/
MultivarProbDistShPtr GTR::getRelRatePrior()
 	{
	return rel_rate_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `rel_rate_param_prior'.
*/
ProbDistShPtr GTR::getRelRateParamPrior()
 	{
	return rel_rate_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rel_rate_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void GTR::setRelRatePrior(MultivarProbDistShPtr d)
 	{
	rel_rate_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `rel_rate_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void GTR::setRelRateParamPrior(ProbDistShPtr d)
 	{
	rel_rate_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all four nucleotide frequency parameters. The base class version is
|	called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member
|	knows about the change in base frequencies.
*/
void GTR::setNucleotideFreqs(
  double freqA,				/**< the new frequency of base A */
  double freqC,				/**< the new frequency of base C */
  double freqG,				/**< the new frequency of base G */
  double freqT)				/**< the new frequency of base T/U */
	{
	Model::setNucleotideFreqs(freqA, freqC, freqG, freqT);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all four nucleotide frequency parameters to 0.25. The base class
|	version is called to do most of the work, and this function is responsible only for ensuring that the `q_matrix'
|	data member knows about the change in base frequencies.
*/
void GTR::setAllFreqsEqual()
	{
	Model::setAllFreqsEqual();
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets one of the state frequency parameters (the unnormalized values that
|	determine the values in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'.
|	The base class version is called to do most of the work, and this function is responsible only for ensuring that the
|	`q_matrix' data member knows about the change in base frequencies.
*/
void GTR::setStateFreqUnnorm(
  unsigned param_index,		/**< the 0-based index into the `state_freq_unnorm' vector of the element to modify */
  double value)				/**< the new value of `state_freq_unnorm'[`param_index'] */
	{
	Model::setStateFreqUnnorm(param_index, value);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. The base class version
|	is called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member
|	knows about the change in base frequencies.
*/
void GTR::setStateFreqsUnnorm(
  const std::vector<double> & values)	/**< the new unnormalized state frequencies */
	{
	Model::setStateFreqsUnnorm(values);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
ProbDistShPtr GTR::getStateFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void GTR::setStateFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_prior'.
*/
MultivarProbDistShPtr GTR::getStateFreqPrior()
 	{
	return freq_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_prior' data member to the supplied MultivariateProbabilityDistribution shared pointer `d'.
*/
void GTR::setStateFreqPrior(MultivarProbDistShPtr d)
 	{
	freq_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for
|	compatability with MrBayes). The GTR model provide additional columns for the six relative rates, the base
|	frequencies, the gamma shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if
|	an invariable sites model is being used).
*/
std::string GTR::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	std::string s;
	string_vect_t::const_iterator it;
	for (it = relrate_name.begin(); it != relrate_name.end(); ++it)
		{
		s += boost::str(boost::format("\t%s") % (*it));
		}
	for (it = freq_name.begin(); it != freq_name.end(); ++it)
		{
		s += boost::str(boost::format("\t%s") % (*it));
		}
	s += Model::paramHeader();
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes). Assumes that
|	the `state_freqs' vector has length 4 and the `rel_rates' vector has length 6.
*/
std::string GTR::paramReport(
  unsigned ndecimals, 						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
	{
	PHYCAS_ASSERT(rel_rates.size() == 6);
	PHYCAS_ASSERT(state_freqs.size() == 4);
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);

	std::string s = str(boost::format(fmt) % rel_rates[0]);
	s += str(boost::format(fmt) % rel_rates[1]);
	s += str(boost::format(fmt) % rel_rates[2]);
	s += str(boost::format(fmt) % rel_rates[3]);
	s += str(boost::format(fmt) % rel_rates[4]);
	s += str(boost::format(fmt) % rel_rates[5]);
	//std::string s = str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f") % rel_rates[0] % rel_rates[1] % rel_rates[2] % rel_rates[3] % rel_rates[4] % rel_rates[5]);

	s += str(boost::format(fmt) % state_freqs[0]);
	s += str(boost::format(fmt) % state_freqs[1]);
	s += str(boost::format(fmt) % state_freqs[2]);
	s += str(boost::format(fmt) % state_freqs[3]);
	//s += str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f") % state_freqs[0] % state_freqs[1] % state_freqs[2] % state_freqs[3]);

	s += Model::paramReport(ndecimals, include_edgelen_hyperparams);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the transition/transversion ratio given the six relative rates stored in the `rel_rates' vector, and the
|	relative base frequencies, which are stored in the `state_freqs' vector. Here are the details of the calculation
|	(for brevity, `state_freqs'[0] has been symbolized piA, `state_freqs'[1] by piC, `state_freqs'[2] by piG and
|	`state_freqs'[3] by piT, `rel_rates'[0] by rAC, `rel_rates'[1] by rAG `rel_rates'[2] by rAT, `rel_rates'[3] by rCG,
|	`rel_rates'[4] by rCT and `rel_rates'[5] by rGT):
|>
|	Pr(any transition | dt)   = Pr(AG) + Pr(CT) + Pr(GA) + Pr(TC)
|	                          = (piA piG rAG dt) + (piC piT rCT dt) + (piG piA rAG dt) + (piT piC rCT dt)
|	                          = 2 dt (piA piG rAG + piT piC rCT dt)
|
|	Pr(any transversion | dt) = Pr(AC) + Pr(AT) + Pr(CA) + Pr(CG) + Pr(GC) + Pr(GT) + Pr(TA) + Pr(TG)
|	                          = (piA piC rAC dt) + (piA piT rAT dt) + (piC piA rAC dt) + (piC piG rCG dt)
|	                            + (piG piC rCG dt) + (piG piT rGT dt) + (piT piA rAT dt) + (piT piG rGT dt)
|	                          = 2 dt (piA piC rAC + piA piT rAT + piC piG rCG +  + piG piT rGT)
|
|	                      2 dt (piA piG rAG + piT piC rCT)
|	TRatio = ------------------------------------------------------------
|	         2 dt (piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT)
|
|	                     piA piG rAG + piT piC rCT
|	       = -----------------------------------------------------
|	         piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT
|<
*/
double GTR::calcTRatio()
	{
	double numerator = state_freqs[0]*state_freqs[2]*rel_rates[1] + state_freqs[3]*state_freqs[1]*rel_rates[4];
    double denominator = state_freqs[0]*state_freqs[1]*rel_rates[0] + state_freqs[0]*state_freqs[3]*rel_rates[2] + state_freqs[1]*state_freqs[2]*rel_rates[3] + state_freqs[2]*state_freqs[3]*rel_rates[5];
	PHYCAS_ASSERT(denominator != 0.0);
	double tratio = numerator/denominator;
	return tratio;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `rel_rates[param_index]'.
*/
double GTR::getRelRateUnnorm(
  unsigned param_index)		/**< the 0-based index into the `gamma_rate_probs' vector of the element to return */
	{
	return rel_rates[param_index];
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns number of free parameters, not including edge lengths.
*/
unsigned GTR::getNumFreeParameters() const
    {
    unsigned n = Model::getNumFreeParameters();
    n += 8; // 5 exchangeabilities and 3 nucleotide frequencies
    return n;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the free parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendTransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void GTR::appendPWKParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendPWKParamNames(names, prefix);
    std::string s;

    s = boost::str(boost::format("log(%srAG/%srAC)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%srAT/%srAC)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%srCG/%srAC)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%srCT/%srAC)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%srGT/%srAC)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%sfreqC/%sfreqA)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%sfreqG/%sfreqA)") % prefix % prefix);
    names.push_back(s);

    s = boost::str(boost::format("log(%sfreqT/%sfreqA)") % prefix % prefix);
    names.push_back(s);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the free parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendTransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void GTR::appendFreeParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendFreeParamNames(names, prefix);
    std::string s;

    s = boost::str(boost::format("%srAG") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srAT") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srCG") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srCT") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srGT") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqC") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqG") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqT") % prefix);
    names.push_back(s);

    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the parameters of this model to the supplied vector `names' in the same order used when
|   untransformed parameter values are appended in the member function appendUntransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void GTR::appendParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendParamNames(names, prefix);
    std::string s;

    s = boost::str(boost::format("%srAC") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srAG") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srAT") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srCG") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srCT") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%srGT") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqA") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqC") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqG") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%sfreqT") % prefix);
    names.push_back(s);

    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the log of the determinant of the log (or log-ratio) transformation used by the
|   appendTransformedParamValues() and setParamValueFromTransformed() methods. If g(theta*) is the density function of
|   the transformed values theta*, and f(theta) is the density of the untransformed values theta, then
|   g(theta*) = f(theta) |J|, where |J| is the Jacobian computed by this function.
*/
double GTR::calcLogDetJacobian() const
	{
    double log_det_jacobian = Model::calcLogDetJacobian();

    // exchangeabilities
    log_det_jacobian += log(rel_rates[0]);
    log_det_jacobian += log(rel_rates[1]);
    log_det_jacobian += log(rel_rates[2]);
    log_det_jacobian += log(rel_rates[3]);
    log_det_jacobian += log(rel_rates[4]);
    log_det_jacobian += log(rel_rates[5]);

    // equilibrium relative state frequencies
    log_det_jacobian += log(state_freqs[0]);
    log_det_jacobian += log(state_freqs[1]);
    log_det_jacobian += log(state_freqs[2]);
    log_det_jacobian += log(state_freqs[3]);

    return log_det_jacobian;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the untransformed values of the parameters of this model to the supplied vector `values' in the same order
|   used when parameter names are appended in the member function appendParamNames(). Derived classes should override
|   this function, calling this version before adding the values of parameters specific to the model encapsulated by the
|   derived class.
*/
void GTR::appendUntransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    Model::appendUntransformedParamValues(values);

    //double sumx = std::accumulate(rel_rates.begin(), rel_rates.end(), 0.0);
    //std::cerr << "sum of exchangeabilities: " << sumx << std::endl;
    //std::cerr << "avg of exchangeabilities: " << (sumx/6.0) << std::endl;

    values.push_back(rel_rates[0]);
    values.push_back(rel_rates[1]);
    values.push_back(rel_rates[2]);
    values.push_back(rel_rates[3]);
    values.push_back(rel_rates[4]);
    values.push_back(rel_rates[5]);

    //sumx = std::accumulate(state_freqs.begin(), state_freqs.end(), 0.0);
    //std::cerr << "sum of state freqs: " << sumx << std::endl;
    //std::cerr << "avg of state freqs: " << (sumx/6.0) << std::endl;

    values.push_back(state_freqs[0]);
    values.push_back(state_freqs[1]);
    values.push_back(state_freqs[2]);
    values.push_back(state_freqs[3]);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the transformed values of the parameters of this model to the supplied vector `values' in the same order
|   used when parameter names are appended in the member function appendFreeParamNames(). Derived classes should override
|   this function, calling this version before adding the values of parameters specific to the model encapsulated by the
|   derived class.
*/
void GTR::appendTransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    Model::appendTransformedParamValues(values);
    double log_value;

    // exchangeabilities
    double log_rAC = log(rel_rates[0]);

    //std::cerr << "In GTR::appendTransformedParamValues: values before adding rel_rates" << std::endl;
    //std::copy(values.begin(), values.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;

    //std::cerr << "length of values = " << values.size() << std::endl;

    //std::cerr << "In GTR::appendTransformedParamValues: rel_rates" << std::endl;
    //std::copy(rel_rates.begin(), rel_rates.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;
    //std::vector<double>::iterator tmp = values.end();

    double log_rAG = log(rel_rates[1]);
    log_value = log_rAG - log_rAC;
    values.push_back(log_value);

    double log_rAT = log(rel_rates[2]);
    log_value = log_rAT - log_rAC;
    values.push_back(log_value);

    double log_rCG = log(rel_rates[3]);
    log_value = log_rCG - log_rAC;
    values.push_back(log_value);

    double log_rCT = log(rel_rates[4]);
    log_value = log_rCT - log_rAC;
    values.push_back(log_value);

    double log_rGT = log(rel_rates[5]);
    log_value = log_rGT - log_rAC;
    values.push_back(log_value);

    //std::cerr << "In GTR::appendTransformedParamValues: values after adding rel_rates" << std::endl;
    //std::copy(values.begin(), values.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;

    //std::cerr << "length of values = " << values.size() << std::endl;

    //std::cerr << "In GTR::appendTransformedParamValues: state_freqs" << std::endl;
    //std::copy(state_freqs.begin(), state_freqs.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;

    // equilibrium relative state frequencies
    double log_freqA = log(state_freqs[0]);

    double log_freqC = log(state_freqs[1]);
    log_value = log_freqC - log_freqA;
    values.push_back(log_value);

    double log_freqG = log(state_freqs[2]);
    log_value = log_freqG - log_freqA;
    values.push_back(log_value);

    double log_freqT = log(state_freqs[3]);
    log_value = log_freqT - log_freqA;
    values.push_back(log_value);

    //std::cerr << "In GTR::appendTransformedParamValues: values after adding state_freqs" << std::endl;
    //std::copy(values.begin(), values.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;

    //std::cerr << "length of values = " << values.size() << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Untransform the value `transformed_value' and use it to set the value of the parameter whose name is supplied in
|   `parameter_name'. Returns true if parameter was found, false if not.
*/
bool GTR::setParamValueFromTransformed(
  std::string parameter_name,   /**< is the name of the parameter to set */
  double transformed_value,     /**< is the transformed value of the parameter to set */
  TreeShPtr tree)               /**< is the tree */
	{
    bool found = Model::setParamValueFromTransformed(parameter_name, transformed_value, tree);
    if (found)
        return true;

    if (parameter_name.compare("rAG") == 0)
        {
        double v = exp(transformed_value);
        setRelRateRatio(0, v);
        return true;
        }

    if (parameter_name.compare("rAT") == 0)
        {
        double v = exp(transformed_value);
        setRelRateRatio(1, v);
        return true;
        }

    if (parameter_name.compare("rCG") == 0)
        {
        double v = exp(transformed_value);
        setRelRateRatio(2, v);
        return true;
        }

    if (parameter_name.compare("rCT") == 0)
        {
        double v = exp(transformed_value);
        setRelRateRatio(3, v);
        return true;
        }

    if (parameter_name.compare("rGT") == 0)
        {
        double v = exp(transformed_value);
        setRelRateRatio(4, v);
        return true;
        }

    if (parameter_name.compare("freqC") == 0)
        {
        double v = exp(transformed_value);
        setStateFreqRatio(0, v);
        return true;
        }

    if (parameter_name.compare("freqG") == 0)
        {
        double v = exp(transformed_value);
        setStateFreqRatio(1, v);
        return true;
        }

    if (parameter_name.compare("freqT") == 0)
        {
        double v = exp(transformed_value);
        setStateFreqRatio(2, v);
        return true;
        }

    return false;
    }
