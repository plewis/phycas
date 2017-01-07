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

#include <cmath>
#include <iostream>
#include "hky.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 4, the base frequencies `piA', `piC', `piG' and `piT' to 0.25, and `kappa'
|	to 1.0.
*/
HKY::HKY()
  : Model(4), kappa_fixed(false), kappa(1.0)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");
	}

HKY::~HKY()
	{
	//std::cerr << "HKY dying..." << std::endl;
	Model::Clear();
	//std::cerr << "  kappa_prior use_count() = " << kappa_prior.use_count() << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "HKY85", "HKY85+G", "HKY85+I" or "HKY85+G+I".
*/
std::string HKY::getModelName() const
	{
	std::string s = "HKY85";
	if (num_gamma_rates > 1)
		s += "+G";
	if (is_pinvar_model)
		s += "+I";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Model::createParameters to create the edge length parameters, the edge length hyperparameter, and any
|	parameters related to rate heterogeneity. This function then adds additional HKY85-specific parameters to the
|	supplied `parameters' vector. This incudes the four base frequencies as well as the transition/transversion rate
|	ratio kappa.
*/
void HKY::createParameters(
  JointPriorManagerShPtr jpm,               /**< is the object that keeps the joint prior density up-to-date */
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams,	/**< is the edge length hyperparameter */
  MCMCUpdaterVect & parameters,				/**< is the vector of model-specific parameters to fill */
  int subset_pos) 							/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added; otherwise, the `edgelens' and `edgelen_hyperparams' vectors returned will be empty */
	{
	Model::createParameters(jpm, t, edgelens, edgelen_hyperparams, parameters, subset_pos);

    subset_index = subset_pos < 0 ? 0 : subset_pos;

	freq_name.clear();
    //if (subset_pos < 0)
    //    {
    //    freq_name.push_back("freqA");
    //    freq_name.push_back("freqC");
    //    freq_name.push_back("freqG");
    //    freq_name.push_back("freqT");
    //    }
    //else
    //    {
    //    unsigned m = subset_pos + 1;
    //    freq_name.push_back(boost::str(boost::format("freqA_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqC_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqG_%d") % m));
    //    freq_name.push_back(boost::str(boost::format("freqT_%d") % m));
    //    }
    unsigned m = subset_index + 1;
    freq_name.push_back(boost::str(boost::format("%d_freqA") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqC") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqG") % m));
    freq_name.push_back(boost::str(boost::format("%d_freqT") % m));

	PHYCAS_ASSERT(!kappa_param);
	kappa_param = MCMCUpdaterShPtr(new KappaParam());
    //std::string nm;
    //if (subset_pos < 0)
    //  nm = "kappa";   //TODO: param names should include subset index even if unpartitioned
    //else
    //  nm = boost::str(boost::format("kappa_%d") % (subset_pos + 1));
    std::string nm = boost::str(boost::format("%d_kappa") % (subset_index + 1));
    kappa_param->setName(nm);
	kappa_param->setStartingValue(4.0);
	kappa_param->setTree(t);
	kappa_param->setPrior(kappa_prior);
    jpm->addUnivariateDistribution(nm, kappa_prior, 4.0);
	if (kappa_fixed)
		kappa_param->fixParameter();
	parameters.push_back(kappa_param);

	PHYCAS_ASSERT(freq_params.empty());

    PHYCAS_ASSERT(freq_param_prior || freq_prior);
    if (freq_param_prior)
        {
        // Only add frequency parameters if freqs will be updated separately
        // The other option is to update the frequencies jointly using the
        // StateFreqMove Metropolis-Hastings move (in which case freq_param
        // will be set and freq_param_prior will be empty)

	    MCMCUpdaterShPtr freqA_param = MCMCUpdaterShPtr(new StateFreqParam(0));
        //if (subset_pos < 0)
        //    nm = "freqA";
        //else
        //    nm = boost::str(boost::format("freqA_%d") % (subset_pos + 1));
        nm = boost::str(boost::format("%d_freqA") % (subset_index + 1));
        freqA_param->setName(nm);
	    freqA_param->setTree(t);
	    freqA_param->setStartingValue(1.0);
	    freqA_param->setPrior(freq_param_prior);
        jpm->addUnivariateDistribution(nm, freq_param_prior, 1.0);
	    if (state_freq_fixed)
		    freqA_param->fixParameter();
	    parameters.push_back(freqA_param);
	    freq_params.push_back(freqA_param);

	    MCMCUpdaterShPtr freqC_param = MCMCUpdaterShPtr(new StateFreqParam(1));
        //if (subset_pos < 0)
        //    nm = "freqC";
        //else
        //    nm = boost::str(boost::format("freqC_%d") % (subset_pos + 1));
        nm = boost::str(boost::format("%d_freqC") % (subset_index + 1));
        freqC_param->setName(nm);
	    freqC_param->setTree(t);
	    freqC_param->setStartingValue(1.0);
	    freqC_param->setPrior(freq_param_prior);
        jpm->addUnivariateDistribution(nm, freq_param_prior, 1.0);
	    if (state_freq_fixed)
		    freqC_param->fixParameter();
	    parameters.push_back(freqC_param);
	    freq_params.push_back(freqC_param);

	    MCMCUpdaterShPtr freqG_param = MCMCUpdaterShPtr(new StateFreqParam(2));
        //if (subset_pos < 0)
        //    nm = "freqG";
        //else
        //    nm = boost::str(boost::format("freqG_%d") % (subset_pos + 1));
        nm = boost::str(boost::format("%d_freqG") % (subset_index + 1));
        freqG_param->setName(nm);
	    freqG_param->setTree(t);
	    freqG_param->setStartingValue(1.0);
	    freqG_param->setPrior(freq_param_prior);
        jpm->addUnivariateDistribution(nm, freq_param_prior, 1.0);
	    if (state_freq_fixed)
		    freqG_param->fixParameter();
	    parameters.push_back(freqG_param);
	    freq_params.push_back(freqG_param);

	    MCMCUpdaterShPtr freqT_param = MCMCUpdaterShPtr(new StateFreqParam(3));
        //if (subset_pos < 0)
        //    nm = "freqT";
        //else
        //    nm = boost::str(boost::format("freqT_%d") % (subset_pos + 1));
        nm = boost::str(boost::format("%d_freqT") % (subset_index + 1));
        freqT_param->setName(nm);
	    freqT_param->setTree(t);
	    freqT_param->setStartingValue(1.0);
	    freqT_param->setPrior(freq_param_prior);
        jpm->addUnivariateDistribution(nm, freq_param_prior, 1.0);
	    if (state_freq_fixed)
		    freqT_param->fixParameter();
	    parameters.push_back(freqT_param);
	    freq_params.push_back(freqT_param);
        }
    else    // freq_prior
        {
        nm = boost::str(boost::format("%d_state_freqs") % (subset_index + 1));
        double_vect_t starting_freqs(4, 0.25);
        jpm->addMultivariateDistribution(nm, freq_prior, starting_freqs);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for
|	compatability with MrBayes). The HKY model provide additional columns for kappa, the base frequencies, the gamma
|	shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if an invariable sites model
|	is being used)
*/
std::string HKY::paramHeader() const
	{
	std::string s;
	s += boost::str(boost::format("\t%s") % kappa_param->getName());
	for (string_vect_t::const_iterator it = freq_name.begin(); it != freq_name.end(); ++it)
		{
		s += boost::str(boost::format("\t%s") % (*it));
		}
	s += Model::paramHeader();
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the pure virtual base class version to generate a string of tab-separated values of model-specific
|	parameters suitable for saving in a sampled parameter file (e.g. like the .p files saved by MrBayes).
*/
std::string HKY::paramReport(
  unsigned ndecimals,						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
	{
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s = str(boost::format(fmt) % kappa);
	s += str(boost::format(fmt) % state_freqs[0]);
	s += str(boost::format(fmt) % state_freqs[1]);
	s += str(boost::format(fmt) % state_freqs[2]);
	s += str(boost::format(fmt) % state_freqs[3]);
	//std::string s = str(boost::format("\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f") % kappa % state_freqs[0] % state_freqs[1] % state_freqs[2] % state_freqs[3]);
	s += Model::paramReport(ndecimals, include_edgelen_hyperparams);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to true. The fixParameter member function of the KappaParam object is either
|	called immediately (if `kappa_param' is a valid pointer) or is called in createParameters (when `kappa_param' is
|	first assigned).
*/
void HKY::fixKappa()
	{
	kappa_fixed = true;
	if (kappa_param)
		kappa_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to false. The freeParameter member function of the KappaParam object is called
|	immediately if `kappa_param' is a valid pointer.
*/
void HKY::freeKappa()
	{
	kappa_fixed = false;
	if (kappa_param)
		kappa_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa'.
*/
double HKY::getKappa()
 	{
	return kappa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' to supplied value `k'. Throws XLikelihood exception if `k' is less than or equal to 0.0.
*/
void HKY::setKappa(double k)
 	{
	++time_stamp;
	if (k <= 0.0)
		throw XLikelihood();
	kappa = k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa_prior'.
*/
ProbDistShPtr HKY::getKappaPrior()
 	{
	return kappa_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void HKY::setKappaPrior(ProbDistShPtr d)
 	{
	kappa_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
ProbDistShPtr HKY::getStateFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void HKY::setStateFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_prior'.
*/
MultivarProbDistShPtr HKY::getStateFreqPrior()
 	{
	return freq_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_prior' data member to the supplied MultivariateProbabilityDistribution shared pointer `d'.
*/
void HKY::setStateFreqPrior(MultivarProbDistShPtr d)
 	{
	freq_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' using supplied value `tratio', which is the transition/transversion ratio, not the ts/tv rate ratio.
|	The quantity `kappa' is related to `tratio' as follows (where piA means `state_freqs'[0], piC means
|	`state_freqs'[1], piC means	`state_freqs'[2] and piT means `state_freqs'[3]):
|>
|	         tratio (piA + piG) (piC + piT)
|	kappa = -------------------------------
|	              (piA piG + piC piT)
|>
|	Throws XLikelihood exception if `tratio' is less than or equal to 0.0.
*/
void HKY::setKappaFromTRatio(double tratio)
 	{
	++time_stamp;
	if (tratio <= 0.0)
		throw XLikelihood();
	double numerator   = tratio*(state_freqs[0] + state_freqs[2])*(state_freqs[1] + state_freqs[3]);
	double denominator = ((state_freqs[0]*state_freqs[2]) + (state_freqs[1]*state_freqs[3]));
	PHYCAS_ASSERT(denominator != 0.0);
	kappa = numerator/denominator;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the transition/transversion ratio given the transition/transversion rate ratio `kappa' and the relative
|	base frequencies, which are stored in the `state_freqs' vector. Here are the details of the calculation (for brevity,
|	`state_freqs'[0] has been symbolized piA, `state_freqs'[1] by piC, `state_freqs'[2] by piG and `state_freqs'[3] by piT):
|>
|	Parameters: b = transversion rate, k = kappa
|	Pr(any transition | dt)   = Pr(AG) + Pr(CT) + Pr(GA) + Pr(TC)
|	                          = (piA piG k b dt) + (piC piT k b dt) + (piG piA k b dt) + (piT piC k b dt)
|	                          = 2 k b dt (piA piG + piC piT)
|
|	Pr(any transversion | dt) = Pr(AC) + Pr(AT) + Pr(CA) + Pr(CG) + Pr(GC) + Pr(GT) + Pr(TA) + Pr(TG)
|	                          = (piA piC b dt) + (piA piT b dt) + (piC piA b dt) + (piC piG b dt)
|	                            + (piG piC b dt) + (piG piT b dt) + (piT piA b dt) + (piT piG b dt)
|	                          = 2 b dt (piA + piG) (piC + piT)
|
|	          2 k b dt (piA piG + piC piT)     k (piA piG + piC piT)
|	TRatio = ------------------------------ = -----------------------
|	         2 b dt (piA + piG) (piC + piT)   (piA + piG) (piC + piT)
|<
*/

double HKY::calcTRatio()
	{
	double numerator   = kappa*((state_freqs[0]*state_freqs[2]) + (state_freqs[1]*state_freqs[3]));
	double denominator = (state_freqs[0] + state_freqs[2])*(state_freqs[1] + state_freqs[3]);
	PHYCAS_ASSERT(denominator != 0.0);
	double tratio = numerator/denominator;
	return tratio;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double HKY::calcUniformizationLambda() const
	{
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];
	double piR = piA + piG;
	double piY = piC + piT;

    // set beta such that edgelen = t
    double beta = 1.0/(2.0*piR*piY + 2.0*kappa*(piA*piG + piC*piT));

    // set lambda to maximim substitution rate + epsilon
    double lambda_epsilon = 0.1;
    double max_rate = beta*(piY + piG*kappa);
    double next_rate = beta*(piR + piT*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piY + piA*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piR + piC*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    return max_rate + lambda_epsilon;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model.
*/
double HKY::calcLMat(double * * lMat) const
	{
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];
	double piR = piA + piG;
	double piY = piC + piT;

    // set beta such that edgelen = t
    double beta = 1.0/(2.0*piR*piY + 2.0*kappa*(piA*piG + piC*piT));

    double lambda = calcUniformizationLambda();

    lMat[0][0] = log(lambda - beta*(piY + piG*kappa));
    lMat[0][1] = log(piC*beta);
    lMat[0][2] = log(piG*kappa*beta);
    lMat[0][3] = log(piT*beta);
    lMat[1][0] = log(piA*beta);
    lMat[1][1] = log(lambda - beta*(piR + piT*kappa));
    lMat[1][2] = log(piG*beta);
    lMat[1][3] = log(piT*kappa*beta);
    lMat[2][0] = log(piA*kappa*beta);
    lMat[2][1] = log(piC*beta);
    lMat[2][2] = log(lambda - beta*(piY + piA*kappa));
    lMat[2][3] = log(piT*beta);
    lMat[3][0] = log(piA*beta);
    lMat[3][1] = log(piC*kappa*beta);
    lMat[3][2] = log(piG*beta);
    lMat[3][3] = log(lambda - beta*(piR + piC*kappa));

    // here are the values on the normal scale in case they are ever needed
    //lMat[0][0] = 1.0 - beta*(piY + piG*kappa)/lambda
    //lMat[0][1] = piC*beta/lambda
    //lMat[0][2] = piG*kappa*beta/lambda
    //lMat[0][3] = piT*beta/lambda
    //lMat[1][0] = piA*beta/lambda
    //lMat[1][1] = 1.0 - beta*(piR + piT*kappa)/lambda
    //lMat[1][2] = piG*beta/lambda
    //lMat[1][3] = piT*kappa*beta/lambda
    //lMat[2][0] = piA*kappa*beta/lambda
    //lMat[2][1] = piC*beta/lambda
    //lMat[2][2] = 1.0 - beta*(piY + piA*kappa)/lambda
    //lMat[2][3] = piT*beta/lambda
    //lMat[3][0] = piA*beta/lambda
    //lMat[3][1] = piC*kappa*beta/lambda
    //lMat[3][2] = piG*beta/lambda
    //lMat[3][3] = 1.0 - beta*(piR + piC*kappa)/lambda

    return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model.
*/
double HKY::calcUMat(double * * uMat) const
	{
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];
	double piR = piA + piG;
	double piY = piC + piT;

    // set beta such that edgelen = t
    double beta = 1.0/(2.0*piR*piY + 2.0*kappa*(piA*piG + piC*piT));

    // set lambda to maximim substitution rate + epsilon
    double lambda_epsilon = 0.1;
    double max_rate = beta*(piY + piG*kappa);
    double next_rate = beta*(piR + piT*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piY + piA*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    next_rate = beta*(piR + piC*kappa);
    if (next_rate > max_rate)
        max_rate = next_rate;
    double lambda = max_rate + lambda_epsilon;

    uMat[0][0] = 1.0 - beta*(piY + piG*kappa)/lambda;
    uMat[0][1] = piC*beta/lambda;
    uMat[0][2] = piG*kappa*beta/lambda;
    uMat[0][3] = piT*beta/lambda;
    uMat[1][0] = piA*beta/lambda;
    uMat[1][1] = 1.0 - beta*(piR + piT*kappa)/lambda;
    uMat[1][2] = piG*beta/lambda;
    uMat[1][3] = piT*kappa*beta/lambda;
    uMat[2][0] = piA*kappa*beta/lambda;
    uMat[2][1] = piC*beta/lambda;
    uMat[2][2] = 1.0 - beta*(piY + piA*kappa)/lambda;
    uMat[2][3] = piT*beta/lambda;
    uMat[3][0] = piA*beta/lambda;
    uMat[3][1] = piC*kappa*beta/lambda;
    uMat[3][2] = piG*beta/lambda;
    uMat[3][3] = 1.0 - beta*(piR + piC*kappa)/lambda;

    return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from
|	the base class Model.
*/
void HKY::calcPMat(double * * pMat, double edgeLength) const
	{
	PHYCAS_ASSERT(state_freqs.size() == 4);
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];

	double PiA = piA + piG;
	double PiC = piC + piT;
	double PiG = piA + piG;
	double PiT = piC + piT;

	double bigPiInvA = 1.0/PiA;
	double bigPiInvC = 1.0/PiC;
	double bigPiInvG = 1.0/PiG;
	double bigPiInvT = 1.0/PiT;

	double t = edgeLength;
    // The next two lines fix the "Rota" bug; see BUGS file for details
    if (t < 1.e-8)
        t = 1.e-8; //TreeNode::edgeLenEpsilon;
	double ta, tb, tc, td, y;
	double denom = ((piA + piG)*(piC + piT) + kappa*((piA*piG) + (piC*piT)));
	double beta = 0.5/denom;
	double x = exp(-beta*t);

	// changes to base A
	td			= -beta*(1 + PiA*(kappa - 1.0));
	y			= exp(t*td);
	ta			= piA*(bigPiInvA - 1.0);
	tb			= (PiA - piA)*bigPiInvA;
	tc			= piA*bigPiInvA;
	pMat[0][0]	= piA + (x*ta) + (y*tb);
	pMat[1][0]	= piA*(1.0 - x);
	pMat[2][0]	= piA + (x*ta) - (y*tc);
	pMat[3][0]	= pMat[1][0];

	// changes to base C
	td = -beta*(1 + PiC*(kappa - 1.0));
	y = exp(t*td);
	ta = piC*(bigPiInvC - 1.0);
	tb = (PiC - piC)*bigPiInvC;
	tc = piC*bigPiInvC;
	pMat[0][1] = piC*(1.0 - x);
	pMat[1][1] = piC + (x*ta) + (y*tb);
	pMat[2][1] = pMat[0][1];
	pMat[3][1] = piC + (x*ta) - (y*tc);

	// changes to base G
	td = -beta*(1 + PiG*(kappa - 1.0));
	y = exp(t*td);
	ta = piG*(bigPiInvG - 1.0);
	tb = (PiG - piG)*bigPiInvG;
	tc = piG*bigPiInvG;
	pMat[0][2] = piG + (x*ta) - (y*tc);
	pMat[1][2] = piG*(1.0 - x);
	pMat[2][2] = piG + (x*ta) + (y*tb);
	pMat[3][2] = pMat[1][2];

	// changes to base T
	td = -beta*(1 + PiT*(kappa - 1.0));
	y = exp(t*td);
	ta = piT*(bigPiInvT - 1.0);
	tb = (PiT - piT)*bigPiInvT;
	tc = piT*bigPiInvT;
	pMat[0][3] = piT*(1.0 - x);
	pMat[1][3] = piT + (x*ta) - (y*tc);
	pMat[2][3] = pMat[0][3];
	pMat[3][3] = piT + (x*ta) + (y*tb);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns number of free parameters, not including edge lengths.
*/
unsigned HKY::getNumFreeParameters() const
    {
    unsigned n = Model::getNumFreeParameters();
    n += 4; // kappa and 3 nucleotide frequencies
    return n;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the free parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendTransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void HKY::appendPWKParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendPWKParamNames(names, prefix);

    std::string s;

    s = boost::str(boost::format("log(%skappa)") % prefix);
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
void HKY::appendFreeParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendFreeParamNames(names, prefix);

    std::string s;

    s = boost::str(boost::format("%skappa") % prefix);
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
void HKY::appendParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendParamNames(names, prefix);

    std::string s;

    s = boost::str(boost::format("%skappa") % prefix);
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
double HKY::calcLogDetJacobian() const
	{
    double log_det_jacobian = Model::calcLogDetJacobian();

    log_det_jacobian += log(kappa);

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
void HKY::appendUntransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    Model::appendUntransformedParamValues(values);

    PHYCAS_ASSERT(kappa > 0.0);
    values.push_back(kappa);

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
void HKY::appendTransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    Model::appendTransformedParamValues(values);
    double log_value;

    PHYCAS_ASSERT(kappa > 0.0);
    log_value = log(kappa);
    values.push_back(log_value);

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
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Untransform the value `transformed_value' and use it to set the value of the parameter whose name is supplied in
|   `parameter_name'. Returns true if parameter was found, false if not.
*/
bool HKY::setParamValueFromTransformed(
  std::string parameter_name,   /**< is the name of the parameter to set */
  double transformed_value,     /**< is the transformed value of the parameter to set */
  TreeShPtr tree)               /**< is the tree */
	{
    bool found = Model::setParamValueFromTransformed(parameter_name, transformed_value, tree);
    if (found)
        return true;

    if (parameter_name.compare("kappa") == 0)
        {
        double v = exp(transformed_value);
        //kappa_param->sendCurrValueToModel(v);
        setKappa(v);
        if (_joint_prior_manager)
            _joint_prior_manager->univariateModified(boost::str(boost::format("%d_kappa") % (subset_index + 1)), v);
        return true;
        }

    //    if (parameter_name.compare("freqA") == 0)
    //        {
    //        double v = exp(transformed_value);
    //        setStateFreqUnnorm(0, v);
    //        if (_joint_prior_manager)
    //            {
    //            PHYCAS_ASSERT(freq_param_prior);
    //            _joint_prior_manager->univariateModified(boost::str(boost::format("%d_freqA") % (subset_index + 1)), v);
    //            }
    //        return true;
    //        }

    if (parameter_name.compare("freqC") == 0)
        {
        double v = exp(transformed_value);
        setStateFreqRatio(0, v);
        if (_joint_prior_manager)
            {
            if (freq_param_prior)
                _joint_prior_manager->univariateModified(boost::str(boost::format("%d_freqC") % (subset_index + 1)), v);
            else
                setStateFreqRatio(0, v);
            }
        return true;
        }

    if (parameter_name.compare("freqG") == 0)
        {
        double v = exp(transformed_value);
        setStateFreqRatio(1, v);
        if (_joint_prior_manager)
            {
            if (freq_param_prior)
                _joint_prior_manager->univariateModified(boost::str(boost::format("%d_freqG") % (subset_index + 1)), v);
            else
                setStateFreqRatio(1, v);
            }
        return true;
        }

    if (parameter_name.compare("freqT") == 0)
        {
        double v = exp(transformed_value);
        setStateFreqRatio(2, v);
        if (_joint_prior_manager)
            {
            if (freq_param_prior)
                _joint_prior_manager->univariateModified(boost::str(boost::format("%d_freqT") % (subset_index + 1)), v);
            else
                setStateFreqRatio(2, v);

            // This should be the last ratio to be encountered
            calcStateFreqsFromRatios(); // calls _joint_prior_manager->multivariateModified
            }
        return true;
        }

    return false;
    }
