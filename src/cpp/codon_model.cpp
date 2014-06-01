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
#include "ncl/nxsallocatematrix.h"
#include "codon_model.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

static const char * codon_code_to_triplet[] =
	{
	"AAA",
	"AAC",
	"AAG",
	"AAT",
	"ACA",
	"ACC",
	"ACG",
	"ACT",
	"AGA",
	"AGC",
	"AGG",
	"AGT",
	"ATA",
	"ATC",
	"ATG",
	"ATT",
	"CAA",
	"CAC",
	"CAG",
	"CAT",
	"CCA",
	"CCC",
	"CCG",
	"CCT",
	"CGA",
	"CGC",
	"CGG",
	"CGT",
	"CTA",
	"CTC",
	"CTG",
	"CTT",
	"GAA",
    "GAC",
	"GAG",
	"GAT",
	"GCA",
	"GCC",
	"GCG",
	"GCT",
	"GGA",
	"GGC",
	"GGG",
	"GGT",
	"GTA",
	"GTC",
	"GTG",
	"GTT",
	"TAC",
	"TAT",
	"TCA",
	"TCC",
	"TCG",
	"TCT",
	"TGC",
	"TGG",
	"TGT",
	"TTA",
	"TTC",
	"TTG",
	"TTT"
	};

#if 0
static char * codon_code_to_aa_name[] =
	{
	"Lys", //  0 AAA
	"Asn", //  1 AAC
	"Lys", //  2 AAG
	"Asn", //  3 AAT
	"Thr", //  4 ACA
	"Thr", //  5 ACC
	"Thr", //  6 ACG
	"Thr", //  7 ACT
	"Arg", //  8 AGA
	"Ser", //  9 AGC
	"Arg", // 10 AGG
	"Ser", // 11 AGT
	"Ile", // 12 ATA
	"Ile", // 13 ATC
	"Met", // 14 ATG
	"Ile", // 15 ATT
	"Gln", // 16 CAA
	"His", // 17 CAC
	"Gln", // 18 CAG
	"His", // 19 CAT
	"Pro", // 20 CCA
	"Pro", // 21 CCC
	"Pro", // 22 CCG
	"Pro", // 23 CCT
	"Arg", // 24 CGA
	"Arg", // 25 CGC
	"Arg", // 26 CGG
	"Arg", // 27 CGT
	"Leu", // 28 CTA
	"Leu", // 29 CTC
	"Leu", // 30 CTG
	"Leu", // 31 CTT
	"Glu", // 32 GAA
	"Asp", // 33 GAC
	"Glu", // 34 GAG
	"Asp", // 35 GAT
	"Ala", // 36 GCA
	"Ala", // 37 GCC
	"Ala", // 38 GCG
	"Ala", // 39 GCT
	"Gly", // 40 GGA
	"Gly", // 41 GGC
	"Gly", // 42 GGG
	"Gly", // 43 GGT
	"Val", // 44 GTA
	"Val", // 45 GTC
	"Val", // 46 GTG
	"Val", // 47 GTT
	"Tyr", // 48 TAC
	"Tyr", // 49 TAT
	"Ser", // 50 TCA
	"Ser", // 51 TCC
	"Ser", // 52 TCG
	"Ser", // 53 TCT
	"Cys", // 54 TGC
	"Trp", // 55 TGG
	"Cys", // 56 TGT
	"Leu", // 57 TTA
	"Phe", // 58 TTC
	"Leu", // 59 TTG
	"Phe"  // 60 TTT
	};
#endif

//    0    Ala     A       Alanine
//    1    Arg     R       Arginine
//    2    Asn     N       Asparagine
//    3    Asp     D       Aspartic acid (Aspartate)
//    4    Cys     C       Cysteine
//    5    Gln     Q       Glutamine
//    6    Glu     E       Glutamic acid (Glutamate)
//    7    Gly     G       Glycine
//    8    His     H       Histidine
//    9    Ile     I       Isoleucine
//   10    Leu     L       Leucine
//   11    Lys     K       Lysine
//   12    Met     M       Methionine
//   13    Phe     F       Phenylalanine
//   14    Pro     P       Proline
//   15    Ser     S       Serine
//   16    Thr     T       Threonine
//   17    Trp     W       Tryptophan
//   18    Tyr     Y       Tyrosine
//   19    Val     V       Valine

static unsigned codon_code_to_aa_code[] =
	{
	11, //  0 AAA "Lys"
	2, //  1 AAC "Asn"
	11, //  2 AAG "Lys"
	2, //  3 AAT "Asn"
	16, //  4 ACA "Thr"
	16, //  5 ACC "Thr"
	16, //  6 ACG "Thr"
	16, //  7 ACT "Thr"
	1, //  8 AGA "Arg"
	15, //  9 AGC "Ser"
	1, // 10 AGG "Arg"
	15, // 11 AGT "Ser"
	9, // 12 ATA "Ile"
	9, // 13 ATC "Ile"
	12, // 14 ATG "Met"
	9, // 15 ATT "Ile"
	5, // 16 CAA "Gln"
	8, // 17 CAC "His"
	5, // 18 CAG "Gln"
	8, // 19 CAT "His"
	14, // 20 CCA "Pro"
	14, // 21 CCC "Pro"
	14, // 22 CCG "Pro"
	14, // 23 CCT "Pro"
	1, // 24 CGA "Arg"
	1, // 25 CGC "Arg"
	1, // 26 CGG "Arg"
	1, // 27 CGT "Arg"
	10, // 28 CTA "Leu"
	10, // 29 CTC "Leu"
	10, // 30 CTG "Leu"
	10, // 31 CTT "Leu"
	6, // 32 GAA "Glu"
    3, // 33 GAC "Asp"
	6, // 34 GAG "Glu"
	3, // 35 GAT "Asp"
	0, // 36 GCA "Ala"
	0, // 37 GCC "Ala"
	0, // 38 GCG "Ala"
	0, // 39 GCT "Ala"
	7, // 40 GGA "Gly"
	7, // 41 GGC "Gly"
	7, // 42 GGG "Gly"
	7, // 43 GGT "Gly"
	19, // 44 GTA "Val"
	19, // 45 GTC "Val"
	19, // 46 GTG "Val"
	19, // 47 GTT "Val"
	18, // 48 TAC "Tyr"
	18, // 49 TAT "Tyr"
	15, // 50 TCA "Ser"
	15, // 51 TCC "Ser"
	15, // 52 TCG "Ser"
	15, // 53 TCT "Ser"
	4, // 54 TGC "Cys"
	17, // 55 TGG "Trp"
	4, // 56 TGT "Cys"
	10, // 57 TTA "Leu"
	13, // 58 TTC "Phe"
	10, // 59 TTG "Leu"
	13  // 60 TTT "Phe"
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 61, the codon frequencies to 1/61, and `kappa' and `omega' both to 1.0.
*/
Codon::Codon()
  : Model(61), kappa_fixed(false), omega_fixed(false), kappa(1.0), omega(1.0)
	{
	is_codon_model = true;
	state_repr.reserve(64);
	const char * bases[] = {"A", "C", "G", "T"};
	for (unsigned i = 0; i < 4; ++i)
		{
		for (unsigned j = 0; j < 4; ++j)
			{
			for (unsigned k = 0; k < 4; ++k)
				{
				std::string s = str(boost::format("%c%c%c") % bases[i] % bases[j] % bases[k]);
				state_repr.push_back(s);
				}
			}
		}

	std::vector<std::string>::const_iterator it;

	// ignore stop codons
	it = state_repr.begin() + 48;
	state_repr.erase(state_repr.begin()+48);

	it = state_repr.begin() + 49;
	state_repr.erase(state_repr.begin()+49);

	it = state_repr.begin() + 54;
	state_repr.erase(state_repr.begin()+54);

    freq_label.clear();
    for (string_vect_t::iterator frqit = state_repr.begin(); frqit != state_repr.end(); ++frqit)
        {
        std::string s = "freq" + (*frqit);
        freq_label.push_back(s);
        }

	setAllFreqsEqual();
	updateQMatrix();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "Codon", "Codon+G", "Codon+I" or "Codon+G+I".
*/
std::string Codon::getModelName() const
	{
	std::string s = "Codon";
	if (num_gamma_rates > 1)
		s += "+G";
	if (is_pinvar_model)
		s += "+I";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Model::createParameters to create the edge length parameters, the edge length hyperparameter, and any
|	parameters related to rate heterogeneity. This function then adds additional codon-model-specific parameters to the
|	supplied `parameters' vector. This incudes the codon frequencies as well as the transition/transversion rate
|	ratio kappa and the nonsynonymous/synonymous rate ratio omega.
*/
void Codon::createParameters(
  JointPriorManagerShPtr jpm,               /**< is the object that keeps the joint prior density up-to-date */
  TreeShPtr t,								/**< is the tree (the nodes of which are needed for creating edge length parameters) */
  MCMCUpdaterVect & edgelens,				/**< is the vector of edge length parameters to fill */
  MCMCUpdaterVect & edgelen_hyperparams,	/**< is the vector of edge length hyperparameters to fill */
  MCMCUpdaterVect & parameters,				/**< is the vector of model-specific parameters to fill */
  int subset_pos)							/**< if 0 (first subset) or -1 (only subset), edge length parameters and hyperparams will be added; otherwise, the `edgelens' and `edgelen_hyperparams' vectors returned will be empty */
	{
	Model::createParameters(jpm, t, edgelens, edgelen_hyperparams, parameters, subset_pos);

    subset_index = subset_pos < 0 ? 0 : subset_pos;

	PHYCAS_ASSERT(!kappa_param);
	kappa_param = MCMCUpdaterShPtr(new KappaParam());
    std::string nm = boost::str(boost::format("%d_kappa") % (subset_index + 1));
    kappa_param->setName(nm);
	kappa_param->setStartingValue(kappa);
	kappa_param->setTree(t);
	kappa_param->setPrior(kappa_prior);
    jpm->addUnivariateDistribution(nm, kappa_prior, kappa);
	if (kappa_fixed)
		kappa_param->fixParameter();
	parameters.push_back(kappa_param);

	PHYCAS_ASSERT(!omega_param);
	omega_param = MCMCUpdaterShPtr(new OmegaParam());
    //if (subset_pos < 0)
    //    nm = "omega";
    //else
    //    nm = boost::str(boost::format("omega_%d") % (subset_pos + 1));
	nm = boost::str(boost::format("%d_omega") % (subset_index + 1));
    omega_param->setName(nm);
	omega_param->setStartingValue(1.0);   // change to omega
	omega_param->setTree(t);
	omega_param->setPrior(omega_prior);
    jpm->addUnivariateDistribution(nm, omega_prior, 0.05);   // change to omega
	if (omega_fixed)
		omega_param->fixParameter();
	parameters.push_back(omega_param);

	PHYCAS_ASSERT(freq_params.empty());
    PHYCAS_ASSERT(freq_param_prior || freq_prior);

    if (freq_param_prior)
        {
        // Only add frequency parameters if freqs will be updated separately
        // The other option is to update the frequencies jointly using the
        // StateFreqMove Metropolis-Hastings move (in which case freq_prior
        // will be set and freq_param_prior will be empty)
		freq_name.clear();
		for (unsigned i = 0; i < 61; ++i)
			{
			MCMCUpdaterShPtr state_freq_param = MCMCUpdaterShPtr(new StateFreqParam(i));
            //if (subset_pos < 0)
            //    {
            //    nm = boost::str(boost::format("freq%s") % state_repr[i]);
            //    freq_name.push_back(nm);
            //    state_freq_param->setName(nm);
            //    }
            //else
            //    {
            //    nm = boost::str(boost::format("freq%s_%d") % state_repr[i] % (subset_pos + 1));
            //    freq_name.push_back(nm);
            //    state_freq_param->setName(nm);
            //    }
            nm = boost::str(boost::format("%d_freq%s") % state_repr[i] % (subset_index + 1));
            freq_name.push_back(nm);
            state_freq_param->setName(nm);
			state_freq_param->setTree(t);
			state_freq_param->setStartingValue(1.0);
			state_freq_param->setPrior(freq_param_prior);
            nm = boost::str(boost::format("%d_freq%s") % (subset_index + 1) % state_repr[i]);
            jpm->addUnivariateDistribution(nm, freq_param_prior, 1.0);
			if (state_freq_fixed)
				state_freq_param->fixParameter();
			parameters.push_back(state_freq_param);
			freq_params.push_back(state_freq_param);
			}
		}
	else
		{
        // Frequencies are being updated jointly using the StateFreqMove Metropolis-Hastings move, but we still
		// must report the frequencies in the param file, and thus will need names of all the state frequencies
		// to use as headers
		freq_name.clear();
		for (unsigned i = 0; i < 61; ++i)
			{
            //if (subset_pos < 0)
            //    {
            //    nm = boost::str(boost::format("freq%s") % state_repr[i]);
            //    freq_name.push_back(nm);
            //    }
            //else
            //    {
            //    nm = boost::str(boost::format("freq%s_%d") % state_repr[i] % (subset_pos + 1));
            //    freq_name.push_back(nm);
            //    }
            nm = boost::str(boost::format("%d_freq%s") % (subset_index + 1) % state_repr[i]);
            freq_name.push_back(nm);
			}
        nm = boost::str(boost::format("%d_state_freqs") % (subset_index + 1));
        double_vect_t starting_freqs(61, 1.0/61.0);
        jpm->addMultivariateDistribution(nm, freq_prior, starting_freqs);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for
|	compatability with MrBayes). The Codon model provide additional columns for kappa, omega, the codon frequencies,
|	the gamma shape parameter (if the number of rates is greater than 1) and the pinvar parameter (if an invariable
|	sites model is being used)
*/
std::string Codon::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	std::string s;
	s += boost::str(boost::format("\t%s") % kappa_param->getName());
	s += boost::str(boost::format("\t%s") % omega_param->getName());
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
std::string Codon::paramReport(
  unsigned ndecimals,						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
	{
    std::string fmt = boost::str(boost::format("%%.%df\t") % ndecimals);
	std::string s = boost::str(boost::format(fmt) % kappa);
	s += boost::str(boost::format(fmt) % omega);
	//std::string s = boost::str(boost::format("\t%.5f\t%.5f\t") % kappa % omega);
	for (unsigned i = 0; i < 61; ++i)
		{
		s += str(boost::format(fmt) % state_freqs[i]);
		//s += str(boost::format("%.5f\t") % state_freqs[i]);
		}
	s += Model::paramReport(ndecimals, include_edgelen_hyperparams);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to true. The fixParameter member function of the KappaParam object is either
|	called immediately (if `kappa_param' is a valid pointer) or is called in createParameters (when `kappa_param' is
|	first assigned).
*/
void Codon::fixKappa()
	{
	kappa_fixed = true;
	if (kappa_param)
		kappa_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `kappa_fixed' to false. The freeParameter member function of the KappaParam object is called
|	immediately if `kappa_param' is a valid pointer.
*/
void Codon::freeKappa()
	{
	kappa_fixed = false;
	if (kappa_param)
		kappa_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `omega_fixed' to true. The fixParameter member function of the OmegaParam object is either
|	called immediately (if `omega_param' is a valid pointer) or is called in createParameters (when `omega_param' is
|	first assigned).
*/
void Codon::fixOmega()
	{
	omega_fixed = true;
	if (omega_param)
		omega_param->fixParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `omega_fixed' to false. The freeParameter member function of the OmegaParam object is called
|	immediately if `omega_param' is a valid pointer.
*/
void Codon::freeOmega()
	{
	omega_fixed = false;
	if (omega_param)
		omega_param->freeParameter();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa'.
*/
double Codon::getKappa()
 	{
	return kappa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa' to supplied value `k'. Throws XLikelihood exception if `k' is less than or equal to 0.0.
*/
void Codon::setKappa(double k)
 	{
	++time_stamp;
	if (k <= 0.0)
		throw XLikelihood();
	kappa = k;
	updateQMatrix();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `omega'.
*/
double Codon::getOmega()
 	{
	return omega;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `omega' to supplied value `w'. Throws XLikelihood exception if `w' is less than or equal to 0.0.
*/
void Codon::setOmega(double w)
 	{
	++time_stamp;
	if (w <= 0.0)
		throw XLikelihood();
	omega = w;
	updateQMatrix();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `kappa_prior'.
*/
ProbDistShPtr Codon::getKappaPrior()
 	{
	return kappa_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `kappa_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Codon::setKappaPrior(ProbDistShPtr d)
 	{
	kappa_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `omega_prior'.
*/
ProbDistShPtr Codon::getOmegaPrior()
 	{
	return omega_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `omega_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Codon::setOmegaPrior(ProbDistShPtr d)
 	{
	omega_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_param_prior'.
*/
ProbDistShPtr Codon::getStateFreqParamPrior()
 	{
	return freq_param_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of data member `freq_prior'.
*/
MultivarProbDistShPtr Codon::getStateFreqPrior()
 	{
	return freq_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_prior' data member to the supplied MultivariateProbabilityDistribution shared pointer `d'.
*/
void Codon::setStateFreqPrior(MultivarProbDistShPtr d)
 	{
	freq_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `freq_param_prior' data member to the supplied ProbabilityDistribution shared pointer `d'.
*/
void Codon::setStateFreqParamPrior(ProbDistShPtr d)
 	{
	freq_param_prior = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from
|	the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
void Codon::calcPMat(double * * pMat, double edgeLength) const
	{
	q_matrix.recalcPMat(pMat, edgeLength);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double Codon::calcUniformizationLambda() const
	{
    PHYCAS_ASSERT(0);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double Codon::calcLMat(double * * lMat) const
	{
    std::cerr << "Error in Codon::calcLMat: q_matrix does not yet have the required recalcLMat function" << std::endl;
    PHYCAS_ASSERT(0);
	//updateQMatrix();
	//q_matrix.recalcLMat(lMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model. Uses the data member `q_matrix' to perform the calculation.
*/
double Codon::calcUMat(double * * uMat) const
	{
    std::cerr << "Error in Codon::calcUMat: q_matrix does not yet have the required recalcUMat function" << std::endl;
    PHYCAS_ASSERT(0);
	//updateQMatrix();
	//q_matrix.recalcUMat(uMat, edgeLength);
    return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets one of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. The member function
|	normalizeFreqs is called automatically to recalculate `state_freqs'. Assumes `param_index' is less than `num_states'
|	(i.e. a valid index into `state_freq_unnorm) and `value' is non-negative. Most of the work is done by the base
|	class version (i.e. Model::setStateFreqUnnorm); however this override is necessary because of the need to inform
|	`q_matrix' of the change.
*/
void Codon::setStateFreqUnnorm(
  unsigned param_index,		/**< the 0-based index into the `state_freq_unnorm' vector of the element to modify */
  double value)				/**< the new value of `state_freq_unnorm'[`param_index'] */
	{
	Model::setStateFreqUnnorm(param_index, value);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets all of the state frequency parameters (the unnormalized values that determine the values
|	in the `state_freqs' vector when normalized) in the data member vector `state_freq_unnorm'. Most of the work is done
|	by the base class version (i.e. Model::setStateFreqsUnnorm); however this override is necessary because of the need
|	to inform `q_matrix' of the change.
*/
void Codon::setStateFreqsUnnorm(
  const std::vector<double> & values)	/**< the new unnormalized state frequencies */
	{
	Model::setStateFreqsUnnorm(values);
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of Model base class function that sets all state frequencies to 1/num_states. The base class version is
|	called to do most of the work, and this function is responsible only for ensuring that the `q_matrix' data member
|	knows about the change in state frequencies.
*/
void Codon::setAllFreqsEqual()
	{
	Model::setAllFreqsEqual();
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the 61 codon state frequencies to the values expected based on the four base frequencies provided (all four
|	values provided should greater than or equal to 0.0, but do not need to sum to 1.0). The frequency of a codon is,
|	almost, the product of the three component base frequencies: the "almost" qualification being needed because the
|	three stop codons are not included, so each codon frequency must be corrected by dividing by the sum of the 61
|	non-stop three-nucleotide products. For example, if the specified base freqencies were 0.1, 0.2, 0.3 and 0.4, then
|	the frequency of the ACT codon would be the product (0.1)*(0.2)*(0.4) = 0.008, divided by the sum of all 61 such
|	products, which in this case is 0.972, yielding 0.008/0.971 = 0.00823. Note state frequencies set using this
|	function will be obliterated unless state frequencies are fixed using the fixStateFreqs method.
*/
inline void Codon::setNucleotideFreqs(
  double freqA,				/**< the new value of `state_freq_unnorm'[0] (i.e. frequency of base A) */
  double freqC,				/**< the new value of `state_freq_unnorm'[1] (i.e. frequency of base C) */
  double freqG,				/**< the new value of `state_freq_unnorm'[2] (i.e. frequency of base G) */
  double freqT)				/**< the new value of `state_freq_unnorm'[3] (i.e. frequency of base T/U) */
	{
	++time_stamp;
	PHYCAS_ASSERT(num_states == 4);
	PHYCAS_ASSERT(freqA >= 0.0);
	PHYCAS_ASSERT(freqC >= 0.0);
	PHYCAS_ASSERT(freqG >= 0.0);
	PHYCAS_ASSERT(freqT >= 0.0);
	PHYCAS_ASSERT(state_freq_unnorm.size() == 61);
	state_freq_unnorm[0]  = freqA*freqA*freqA;	// 0 AAA
	state_freq_unnorm[1]  = freqA*freqA*freqC;	// 1 AAC
	state_freq_unnorm[2]  = freqA*freqA*freqG;	// 2 AAG
	state_freq_unnorm[3]  = freqA*freqA*freqT;	// 3 AAT
	state_freq_unnorm[4]  = freqA*freqC*freqA;	// 4 ACA
	state_freq_unnorm[5]  = freqA*freqC*freqC;	// 5 ACC
	state_freq_unnorm[6]  = freqA*freqC*freqG;	// 6 ACG
	state_freq_unnorm[7]  = freqA*freqC*freqT;	// 7 ACT
	state_freq_unnorm[8]  = freqA*freqG*freqA;	// 8 AGA
	state_freq_unnorm[9]  = freqA*freqG*freqC;	// 9 AGC
	state_freq_unnorm[10] = freqA*freqG*freqG;	// 10 AGG
	state_freq_unnorm[11] = freqA*freqG*freqT;	// 11 AGT
	state_freq_unnorm[12] = freqA*freqT*freqA;	// 12 ATA
	state_freq_unnorm[13] = freqA*freqT*freqC;	// 13 ATC
	state_freq_unnorm[14] = freqA*freqT*freqG;	// 14 ATG
	state_freq_unnorm[15] = freqA*freqT*freqT;	// 15 ATT
	state_freq_unnorm[16] = freqC*freqA*freqA;	// 16 CAA
	state_freq_unnorm[17] = freqC*freqA*freqC;	// 17 CAC
	state_freq_unnorm[18] = freqC*freqA*freqG;	// 18 CAG
	state_freq_unnorm[19] = freqC*freqA*freqT;	// 19 CAT
	state_freq_unnorm[20] = freqC*freqC*freqA;	// 20 CCA
	state_freq_unnorm[21] = freqC*freqC*freqC;	// 21 CCC
	state_freq_unnorm[22] = freqC*freqC*freqG;	// 22 CCG
	state_freq_unnorm[23] = freqC*freqC*freqT;	// 23 CCT
	state_freq_unnorm[24] = freqC*freqG*freqA;	// 24 CGA
	state_freq_unnorm[25] = freqC*freqG*freqC;	// 25 CGC
	state_freq_unnorm[26] = freqC*freqG*freqG;	// 26 CGG
	state_freq_unnorm[27] = freqC*freqG*freqT;	// 27 CGT
	state_freq_unnorm[28] = freqC*freqT*freqA;	// 28 CTA
	state_freq_unnorm[29] = freqC*freqT*freqC;	// 29 CTC
	state_freq_unnorm[30] = freqC*freqT*freqG;	// 30 CTG
	state_freq_unnorm[31] = freqC*freqT*freqT;	// 31 CTT
	state_freq_unnorm[32] = freqG*freqA*freqA;	// 32 GAA
	state_freq_unnorm[33] = freqG*freqA*freqC;	// 33 GAC
	state_freq_unnorm[34] = freqG*freqA*freqG;	// 34 GAG
	state_freq_unnorm[35] = freqG*freqA*freqT;	// 35 GAT
	state_freq_unnorm[36] = freqG*freqC*freqA;	// 36 GCA
	state_freq_unnorm[37] = freqG*freqC*freqC;	// 37 GCC
	state_freq_unnorm[38] = freqG*freqC*freqG;	// 38 GCG
	state_freq_unnorm[39] = freqG*freqC*freqT;	// 39 GCT
	state_freq_unnorm[40] = freqG*freqG*freqA;	// 40 GGA
	state_freq_unnorm[41] = freqG*freqG*freqC;	// 41 GGC
	state_freq_unnorm[42] = freqG*freqG*freqG;	// 42 GGG
	state_freq_unnorm[43] = freqG*freqG*freqT;	// 43 GGT
	state_freq_unnorm[44] = freqG*freqT*freqA;	// 44 GTA
	state_freq_unnorm[45] = freqG*freqT*freqC;	// 45 GTC
	state_freq_unnorm[46] = freqG*freqT*freqG;	// 46 GTG
	state_freq_unnorm[47] = freqG*freqT*freqT;	// 47 GTT
	state_freq_unnorm[48] = freqT*freqA*freqC;	// 48 TAC
	state_freq_unnorm[49] = freqT*freqA*freqT;	// 49 TAT
	state_freq_unnorm[50] = freqT*freqC*freqA;	// 50 TCA
	state_freq_unnorm[51] = freqT*freqC*freqC;	// 51 TCC
	state_freq_unnorm[52] = freqT*freqC*freqG;	// 52 TCG
	state_freq_unnorm[53] = freqT*freqC*freqT;	// 53 TCT
	state_freq_unnorm[54] = freqT*freqG*freqC;	// 54 TGC
	state_freq_unnorm[55] = freqT*freqG*freqG;	// 55 TGG
	state_freq_unnorm[56] = freqT*freqG*freqT;	// 56 TGT
	state_freq_unnorm[57] = freqT*freqT*freqA;	// 57 TTA
	state_freq_unnorm[58] = freqT*freqT*freqC;	// 58 TTC
	state_freq_unnorm[59] = freqT*freqT*freqG;	// 59 TTG
	state_freq_unnorm[60] = freqT*freqT*freqT;	// 60 TTT
	normalizeFreqs();
	q_matrix.setStateFreqs(state_freqs);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Uses current values of `omega' and `kappa' to update the relative rate (rr) vector inside `q_matrix'. This function
|	should only be called if either `omega' or `kappa' change because it precipitates the relatively expensive
|	eigendecomposition.
*/
void Codon::updateQMatrix() const
	{
	std::vector<double> rel_rate_vect;
	rel_rate_vect.reserve(1830);	// 1830 = 61*(61-1)/2, the number of elements in upper (or lower) triangle

	for (unsigned row = 0; row < 60; ++row)
		{
		unsigned from_aa_state = codon_code_to_aa_code[row];
		for (unsigned col = row + 1; col < 61; ++col)
			{
			// first determine number of nucleotide changes associated with the element at row, col
			unsigned to_aa_state = codon_code_to_aa_code[col];
			unsigned nchanges = 0;
			char from_base = '\0';
			char to_base = '\0';
			for (unsigned k = 0; k < 3; ++k)
				{
				if (codon_code_to_triplet[row][k] != codon_code_to_triplet[col][k])
					{
					from_base = codon_code_to_triplet[row][k];
					to_base = codon_code_to_triplet[col][k];
					++nchanges;
					}
				}

			// second determine the relative rate of the element at row, col
			double rel_rate = 0.0;			// this rate will be used if nchanges > 1
            PHYCAS_ASSERT(nchanges > 0);	// all off-diagonal elements should require at least 1 nucleotide change
			if (nchanges == 1)
				{
				// Determine transition vs. transversion
				bool transversion = true;
				if ((from_base == 'A' && to_base == 'G') || (from_base == 'C' && to_base == 'T') || (from_base == 'G' && to_base == 'A') || (from_base == 'T' && to_base == 'C'))
					transversion = false;

				if (from_aa_state == to_aa_state)
					{
					// synonymous substitution
					rel_rate = (transversion ? 1.0 : kappa);
					}
				else
					{
					// nonsynonymous substitution
					rel_rate = (transversion ? omega : omega*kappa);
					}
				}
			rel_rate_vect.push_back(rel_rate);
			}
		}

	q_matrix.setRelativeRates(rel_rate_vect);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns number of free parameters, not including edge lengths.
*/
unsigned Codon::getNumFreeParameters() const
    {
    unsigned n = Model::getNumFreeParameters();
    n += 62; // kappa, omega, and 60 codon frequencies
    return n;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendTransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void Codon::appendFreeParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendFreeParamNames(names, prefix);
    std::string s;

    s = boost::str(boost::format("%skappa") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%somega") % prefix);
    names.push_back(s);

    string_vect_t::const_iterator it = freq_label.begin();
    for (++it; it != freq_label.end(); ++it)
        {
        s = boost::str(boost::format("%s%s") % prefix % (*it).c_str());
        names.push_back(s);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the parameters of this model to the supplied vector `names' in the same order used when
|   untransformed parameter values are appended in the member function appendUntransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void Codon::appendParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    Model::appendParamNames(names, prefix);
    std::string s;

    s = boost::str(boost::format("%skappa") % prefix);
    names.push_back(s);

    s = boost::str(boost::format("%somega") % prefix);
    names.push_back(s);

    for (string_vect_t::const_iterator it = freq_label.begin(); it != freq_label.end(); ++it)
        {
        s = boost::str(boost::format("%s%s") % prefix % (*it).c_str());
        names.push_back(s);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the log of the determinant of the log (or log-ratio) transformation used by the
|   appendTransformedParamValues() and setParamValueFromTransformed() methods. If g(theta*) is the density function of
|   the transformed values theta*, and f(theta) is the density of the untransformed values theta, then
|   g(theta*) = f(theta) |J|, where |J| is the Jacobian computed by this function.
*/
double Codon::calcLogDetJacobian() const
	{
    double log_det_jacobian = Model::calcLogDetJacobian();

    log_det_jacobian += log(kappa);
    log_det_jacobian += log(omega);

    // equilibrium relative state frequencies
    for (unsigned i = 0; i < 61; ++i)
        log_det_jacobian += log(state_freqs[i]);

    return log_det_jacobian;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the untransformed values of the parameters of this model to the supplied vector `values' in the same order
|   used when parameter names are appended in the member function appendParamNames(). Derived classes should override
|   this function, calling this version before adding the values of parameters specific to the model encapsulated by the
|   derived class.
*/
void Codon::appendUntransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    Model::appendUntransformedParamValues(values);

    values.push_back(kappa);
    values.push_back(omega);

    for (unsigned i = 0; i < 61; ++i)
        values.push_back(state_freqs[i]);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the transformed values of the parameters of this model to the supplied vector `values' in the same order
|   used when parameter names are appended in the member function appendFreeParamNames(). Derived classes should override
|   this function, calling this version before adding the values of parameters specific to the model encapsulated by the
|   derived class.
*/
void Codon::appendTransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    Model::appendTransformedParamValues(values);
    double log_value;

    // equilibrium relative state frequencies
    double log_freqA = log(state_freqs[0]);

    for (unsigned i = 1; i < 61; ++i)
        {
        log_value = log(state_freqs[i]) - log_freqA;
        values.push_back(log_value);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Untransform the value `transformed_value' and use it to set the value of the parameter whose name is supplied in
|   `parameter_name'. Returns true if parameter was found, false if not.
*/
bool Codon::setParamValueFromTransformed(
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

    if (parameter_name.compare("omega") == 0)
        {
        double v = exp(transformed_value);
        //kappa_param->sendCurrValueToModel(v);
        setOmega(v);
        if (_joint_prior_manager)
            _joint_prior_manager->univariateModified(boost::str(boost::format("%d_omega") % (subset_index + 1)), v);
        return true;
        }

    string_vect_t::iterator it = std::find(freq_label.begin(), freq_label.end(), parameter_name);
    if (it != freq_label.end())
        {
        double v = exp(transformed_value);
        unsigned i = (unsigned)(it - freq_label.begin());
        setStateFreqRatio(i, v);
        if (_joint_prior_manager)
            _joint_prior_manager->univariateModified(boost::str(boost::format("%d_%s") % (subset_index + 1) % parameter_name), v);
        return true;
        }

    return false;
    }
