/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include "phycas/force_include.h"
#include <algorithm>
#include <numeric>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include "char_super_matrix.hpp"
#include "basic_tree.hpp"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "tip_data.hpp"
#include "internal_data.hpp"
#include "sim_data.hpp"
#include "phycas_string.hpp"
#include "basic_lot.hpp"
#include "edge_iterators.hpp"
#include "partition_model.hpp"
#include "char_super_matrix.hpp"
#include "edge_endpoints.hpp"

static int8_t codon_state_codes[] =
	{
	0,	// 0 AAA
	1,	// 1 AAC
	2,	// 2 AAG
	3,	// 3 AAT
	4,	// 4 ACA
	5,	// 5 ACC
	6,	// 6 ACG
	7,	// 7 ACT
	8,	// 8 AGA
	9,	// 9 AGC
	10, // 10 AGG
	11, // 11 AGT
	12, // 12 ATA
	13, // 13 ATC
	14, // 14 ATG
	15, // 15 ATT
	16, // 16 CAA
	17, // 17 CAC
	18, // 18 CAG
	19, // 19 CAT
	20, // 20 CCA
	21, // 21 CCC
	22, // 22 CCG
	23, // 23 CCT
	24, // 24 CGA
	25, // 25 CGC
	26, // 26 CGG
	27, // 27 CGT
	28, // 28 CTA
	29, // 29 CTC
	30, // 30 CTG
	31, // 31 CTT
	32, // 32 GAA
	33, // 33 GAC
	34, // 34 GAG
	35, // 35 GAT
	36, // 36 GCA
	37, // 37 GCC
	38, // 38 GCG
	39, // 39 GCT
	40, // 40 GGA
	41, // 41 GGC
	42, // 42 GGG
	43, // 43 GGT
	44, // 44 GTA
	45, // 45 GTC
	46, // 46 GTG
	47, // 47 GTT
	61, // 48 TAA stop
	48, // 49 TAC
	61, // 50 TAG stop
	49, // 51 TAT
	50, // 52 TCA
	51, // 53 TCC
	52, // 54 TCG
	53, // 55 TCT
	61, // 56 TGA stop
	54, // 57 TGC
	55, // 58 TGG
	56, // 59 TGT
	57, // 60 TTA
	58, // 61 TTC
	59, // 62 TTG
	60	// 63 TTT
	};

/*
static char * codon_labels[] =
	{
	"AAA Lys",	// 0
	"AAC Asn",	// 1
	"AAG Lys",	// 2
	"AAT Asn",	// 3
	"ACA Thr",	// 4
	"ACC Thr",	// 5
	"ACG Thr",	// 6
	"ACT Thr",	// 7
	"AGA Arg",	// 8
	"AGC Ser",	// 9
	"AGG Arg",	// 10
	"AGT Ser",	// 11
	"ATA Leu",	// 12
	"ATC Leu",	// 13
	"ATG Met",	// 14
	"ATT Leu",	// 15
	"CAA Gln",	// 16
	"CAC His",	// 17
	"CAG Gln",	// 18
	"CAT His",	// 19
	"CCA Pro",	// 20
	"CCC Pro",	// 21
	"CCG Pro",	// 22
	"CCT Pro",	// 23
	"CGA Arg",	// 24
	"CGC Arg",	// 25
	"CGG Arg",	// 26
	"CGT Arg",	// 27
	"CTA Leu",	// 28
	"CTC Leu",	// 29
	"CTG Leu",	// 30
	"CTT Leu",	// 31
	"GAA Glu",	// 32
	"GAC Asp",	// 33
	"GAG Glu",	// 34
	"GAT Asp",	// 35
	"GCA Ala",	// 36
	"GCC Ala",	// 37
	"GCG Ala",	// 38
	"GCT Ala",	// 39
	"GGA Gly",	// 40
	"GGC Gly",	// 41
	"GGG Gly",	// 42
	"GGT Gly",	// 43
	"GTA Val",	// 44
	"GTC Val",	// 45
	"GTG Val",	// 46
	"GTT Val",	// 47
	"TAC Tyr",	// 48
	"TAT Tyr",	// 49
	"TCA Ser",	// 50
	"TCC Ser",	// 51
	"TCG Ser",	// 52
	"TCT Ser",	// 53
	"TGC Cys",	// 54
	"TGG Trp",	// 55
	"TGT Cys",	// 56
	"TTA Leu",	// 57
	"TTC Phe",	// 58
	"TTG Leu",	// 59
	"TTT Phe"	// 60
	};
*/

namespace phycas
{

// **************************************************************************************
// ***** Former TreeLikelihood inlines (begin) ******************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	TreeLikelihood constructor.
*/
TreeLikelihood::TreeLikelihood(
  PartitionModelShPtr mod)		/**< is the partition model */
  :
  likelihood_root(0),
  store_site_likes(false),
  no_data(false),
  partition_model(mod),
  debugging_now(false),
  nevals(0),
  _num_taxa(0),
  _num_char(0)
    {
    unsigned num_subsets = partition_model->getNumSubsets();
    rate_means.resize(num_subsets);
    rate_probs.resize(num_subsets);
    for (unsigned i = 0; i < num_subsets; ++i)
    	{
    	rate_means[i].assign(partition_model->subset_num_rates[i], 1.0);
    	rate_probs[i].assign(partition_model->subset_num_rates[i], 1.0);
    	partition_model->subset_model[i]->recalcRatesAndProbs(rate_means[i], rate_probs[i]);
		//std::cerr << "Rate means and probs for model " << i << "(" << partition_model->subset_model[i]->getModelName() << ") in TreeLikelihood::TreeLikelihood:" << std::endl;//temp
		//std::copy(rate_means[i].begin(), rate_means[i].end(), std::ostream_iterator<double>(std::cerr," "));//temp
		//std::copy(rate_probs[i].begin(), rate_probs[i].end(), std::ostream_iterator<double>(std::cerr," "));//temp
		//std::cerr << std::endl;
    	}
	cla_pool = CondLikelihoodStorageShPtr(new CondLikelihoodStorage());
	underflow_manager.setTriggerSensitivity(50);
	underflow_manager.setCorrectToValue(10000.0);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	TreeLikelihood destructor.
*/
TreeLikelihood::~TreeLikelihood()
	{
	//std::cerr << "\n>>>>> TreeLikelihood dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Return the `_tl_prior' shared pointer.
*/
TreeLengthDistributionShPtr TreeLikelihood::getTreeLengthPrior()
    {
	return _tl_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Return the `_tl_ref_dist' shared pointer.
*/
TreeLengthDistributionShPtr TreeLikelihood::getTreeLengthRefDist()
    {
	return _tl_ref_dist;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Sets the `_tl_prior' shared pointer.
*/
void TreeLikelihood::setTreeLengthPrior(
  TreeLengthDistributionShPtr tl_prior) /**< is the new tree length prior */
    {
	_tl_prior = tl_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Sets the `_tl_ref_dist' shared pointer.
*/
void TreeLikelihood::setTreeLengthRefDist(
  TreeLengthDistributionShPtr tl_ref_dist) /**< is the new tree length reference distribution */
    {
    _tl_ref_dist = tl_ref_dist;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	NEEDS WORK Use samples in `fitting_sample' to parameterize `ref_dist'. This function is called during generalized
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `ref_dist'.
*/
//void TreeLikelihood::finalizeTreeLenWorkingPrior(
//  LotShPtr lot)     /**< is the random number generator to provide to the newly-created TreeLengthDistribution object */
//    {
//    if (_tl_ref_dist)
//        {
//        // working prior has already been calculated
//        return;
//        }
//
//    // To fit the tree length prior, we need to first build a fitting sample of tree lengths using tl_fitting_sample,
//    // which by now contains edge lengths from all EdgeLenParam objects (not just this one). tl_fitting_sample should
//    // have length equal to the number of samples, and tl_fitting_sample[j] is a vector of edge lengths, one from
//    // each EdgeLenParam objects in existence (i.e. one for each edge). Thus, the sum of each row of tl_fitting_sample
//    // represents a sampled tree length.
//    //
//    // In addition to the tree length, we also need to estimate the magnitude of the Dirichlet parameters for external
//    // edge length proportions and the relative lengths of internal vs. external edges
//
//    // First, summarize the sampled tree lengths
//    double n = 0.0;
//    double sum_tl = 0.0;
//    double ss_tl = 0.0;
//    for (double_vect_vect_t::iterator it = tl_fitting_sample.begin(); it != tl_fitting_sample.end(); ++it)
//        {
//        PHYCAS_ASSERT(std::find(it->begin(), it->end(), -1.0) == it->end());
//        double tl = std::accumulate(it->begin(), it->end(), 0.0);
//        sum_tl += tl;
//        ss_tl += pow(tl, 2.0);
//        n += 1.0;
//        }
//    PHYCAS_ASSERT(n > 1.0);
//    double mean_tl = sum_tl/n;
//    double var_tl = (ss_tl - n*pow(mean_tl, 2.0))/(n - 1.0);
//    double betaT = mean_tl/var_tl;
//    double alphaT = betaT*mean_tl;
//
//    // still working on estimating these...
//    double alpha = 1.0;
//    double c = 1.0;
//    _tl_ref_dist.reset(new TreeLengthDistribution(alphaT, betaT, alpha, c));
//    _tl_ref_dist->SetLot(lot.get());
//    }

/*----------------------------------------------------------------------------------------------------------------------
|	Function used by EdgeLenParams to report an edge length to be added to the tree length fitting sample. Only used if
|   the Rannala-Zhu-Yang tree length prior is in use.
*/
//void TreeLikelihood::educateTreeLenWorkingPrior(
//  TreeShPtr tree,                   /**< is the tree */
//  double edge_length,               /**< is the edge length to add to the tree length fitting sample */
//  unsigned edge_len_index)          /**< is the index of the reporting slave */
//	{
//    unsigned nedges = tree->GetNNodes() - 1;
//    PHYCAS_ASSERT(edge_len_index < nedges);
//    if (tl_fitting_sample.empty())
//        {
//        // Create first row
//        double_vect_t tmp(nedges, -1.0);
//        tmp[edge_len_index] = edge_length;
//        tl_fitting_sample.push_back(tmp);
//        }
//    else
//        {
//        double_vect_vect_t::reverse_iterator last_row = tl_fitting_sample.rbegin();
//        double_vect_t::iterator it = std::find(last_row->begin(), last_row->end(), -1.0);
//        if (it == last_row->end())
//            {
//            // This row has been filled, start a new row
//            double_vect_t tmp(nedges, -1.0);
//            tmp[edge_len_index] = edge_length;
//            tl_fitting_sample.push_back(tmp);
//            }
//        else
//            {
//            // This row has NOT yet been filled
//            PHYCAS_ASSERT((*last_row)[edge_len_index] == -1.0);
//            (*last_row)[edge_len_index] = edge_length;
//            }
//        }
//	}
//

/*----------------------------------------------------------------------------------------------------------------------
|   Resets the `partition_model' shared pointer.
*/
void TreeLikelihood::releasePartitionModel()
    {
	partition_model.reset();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the vector of site likelihoods (data member `site_likelihood') computed in
|   TreeLikelihood::harvestLnLFromValidEdge.
*/
const std::vector<double> & TreeLikelihood::getSiteLikelihoods() const
    {
    return site_likelihood;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the vector of site likelihood underflow correction factors (data member `site_uf') computed
|   TreeLikelihood::harvestLnLFromValidEdge.
*/
const std::vector<double> & TreeLikelihood::getSiteUF() const
    {
    return site_uf;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the vector of pattern counts (data member `pattern_counts').
*/
const std::vector<double> & TreeLikelihood::getPatternCounts() const
    {
    return pattern_counts;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the data member `charIndexToPatternIndex'.
*/
const std::vector<unsigned> & TreeLikelihood::getCharIndexToPatternIndex() const
    {
    return charIndexToPatternIndex;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `store_site_likes'.
*/
bool TreeLikelihood::storingSiteLikelihoods() const
    {
    return store_site_likes;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the current value of the data member `store_site_likes' to true if `yes' is true, and false if `yes' is false.
*/
void TreeLikelihood::storeSiteLikelihoods(
  bool yes)    /**< is the value to which `store_site_likes' should be set */
    {
    store_site_likes = yes;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps InternalData data member `state_time' and edge lengths for the supplied nodes `nd1' and `nd2'. Assumes
|	`nd1' and `nd2' are both internal nodes.
*/
void TreeLikelihood::swapInternalDataAndEdgeLen(
  TreeNode * nd1,	/**< is the first of two nodes whose InternalData structures and edge lengths are to be swapped */
  TreeNode * nd2)	/**< is the second of two nodes whose InternalData structures and edge lengths are to be swapped */
	{
	PHYCAS_ASSERT(nd1->IsInternal());
	PHYCAS_ASSERT(no_data || nd1->GetInternalData() != NULL);
	PHYCAS_ASSERT(nd2->IsInternal());
	PHYCAS_ASSERT(no_data || nd2->GetInternalData() != NULL);

	// Swap edge lengths
	double edgelen = nd1->GetEdgeLen();
	nd1->SetEdgeLen(nd2->GetEdgeLen());
	nd2->SetEdgeLen(edgelen);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the setTriggerSensitivity function of the data member `underflow_manager' to set the number of edges that must
|	be traversed before taking action to prevent underflow.
*/
void TreeLikelihood::setUFNumEdges(
  unsigned nedges)	/**< is the number of edges to traverse before taking action to prevent underflow */
	{
	underflow_manager.setTriggerSensitivity(nedges);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of bytes allocated for each CLA. This equals sizeof(LikeFltType) times the product of the number of
|	patterns, number of rates and number of states. Calls corresponding function of data member `cla_pool' to get the
|	value returned.
*/
unsigned TreeLikelihood::bytesPerCLA() const
	{
	return cla_pool->bytesPerCLA();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of CondLikelihood objects created since the `cla_pool' data member was constructed, or since the
|	last call to the function clearStack of `cla_pool', which resets the value to zero. Calls corresponding function of
|	data member `cla_pool' to get the value returned.
*/
unsigned TreeLikelihood::numCLAsCreated() const
	{
	return cla_pool->numCLAsCreated();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of CondLikelihood objects stored in `cla_pool'. The total number of CLAs currently
|	checked out to the tree can be obtained as TreeLikelihood::numCLAsCreated() minus TreeLikelihood::numCLAsStored().
|	Calls corresponding function of data member `cla_pool' to get the value returned.
*/
unsigned TreeLikelihood::numCLAsStored() const
	{
	return cla_pool->numCLAsStored();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of `likelihood_root'. See TreeLikelihood::useAsLikelihoodRoot for more information about
|	the meaning of the likelihood root.
*/
TreeNode * TreeLikelihood::getLikelihoodRoot()
	{
	return likelihood_root;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of the node currently serving as the likelihood root. If `likelihood_root' is NULL, returns -1
|	instead to indicate that no node is currently designated as the likelihood root. This function was written primarily
|	for use by the TreeViewer.py application for debugging purposes.
*/
int TreeLikelihood::getLikelihoodRootNodeNum() const
	{
	return (likelihood_root ? (int)likelihood_root->GetNodeNumber() : -1);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to `focal_nd' when `avoid' is the likelihood root. This function obtains the correct
|	shared pointer, but if that pointer does not point to an object, it goes to `cla_pool' to get another one.
*/
CondLikelihoodShPtr getCondLikePtr(
  TreeNode * focal_nd,	/**< is the focal node */
  TreeNode * avoid)		/**< is the focal node neighbor (node closer to likelihood root) */
	{
	EdgeEndpoints e(focal_nd, avoid);
	return getCondLikePtr(e);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to `focal_nd' when `avoid' is the likelihood root.
*/
ConstCondLikelihoodShPtr getValidCondLikePtr(
  const TreeNode * focal_nd,	/**< is the focal node */
  const TreeNode * avoid)		/**< is the focal node neighbor (node closer to likelihood root) */
	{
	ConstEdgeEndpoints e(focal_nd, avoid);
	return getValidCondLikePtr(e);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to the focal node of `edge' when `edge' focal neighbor is the likelihood root. This
|	function obtain the correct shared pointer, but if that pointer does not point to a CondLikelihood object, it goes
|	to `cla_pool' to get one.
*/
CondLikelihoodShPtr getCondLikePtr(
  EdgeEndpoints edge) /**< is the edge specifying the focal node and the focal node neighbor */
	{
	TreeNode * actual_child = edge.getActualChild();
	if (actual_child == edge.getFocalNode())
		{
		// focal node F is a child of N, the focal neighbor (likelihood root is somewhere below N)
		//
		//		\	/
		//		 \ /
		//		  F
		//	 \	 / <-- filial CLA of focal node is returned
		//	  \ /
		//	   N
		//	  /
		//
		PHYCAS_ASSERT(actual_child->IsInternal());
		InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getChildCondLikePtr();
		}

	// focal neighbor N is a child of the focal node F
	PHYCAS_ASSERT(actual_child == edge.getFocalNeighbor());
	if (actual_child->IsInternal())
		{
		// focal neighbor N is an internal node (likelihood root is somewhere above N)
		//
		//		\	/
		//		 \ /
		//		  N
		//	 \	 /
		//	  \ / <-- parental CLA of focal neighbor is returned
		//	   F
		//	  /
		//
		InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getParentalCondLikePtr();
		}

	// focal neighbor N is a tip node (N equals the likelihood root in this case)
	//
	//		  N
	//	 \	 /
	//	  \ / <-- parental CLA of focal neighbor is returned
	//	   F
	//	  /
	//
	TipData * child_tip_data = actual_child->GetTipData();
	PHYCAS_ASSERT(child_tip_data != NULL);
	return child_tip_data->getParentalCondLikePtr();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to the focal node of `edge' when `edge' focal neighbor is the likelihood root.
|	This is a const version of the corresponding getCondLikePtr function, to be used when it can be assumed that the
|	conditional likelihood arrays are up-to-date.
*/
ConstCondLikelihoodShPtr getValidCondLikePtr(
  ConstEdgeEndpoints edge) /**< is the edge comprising the focal node and focal node neighbor */
	{
#if 1
	const TreeNode * actual_child = edge.getActualChild();
	if (actual_child == edge.getFocalNode())
		{
		// focal node F is a child of N, the focal neighbor (likelihood root is somewhere below N)
		//
		//		\	/
		//		 \ /
		//		  F
		//	 \	 / <-- filial CLA of focal node is returned
		//	  \ /
		//	   N
		//	  /
		//
		PHYCAS_ASSERT(actual_child->IsInternal());
		const InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getValidChildCondLikePtr();
		}

	// focal neighbor N is a child of the focal node F
	PHYCAS_ASSERT(actual_child == edge.getFocalNeighbor());
	if (actual_child->IsInternal())
		{
		// focal neighbor N is an internal node (likelihood root is somewhere above N)
		//
		//		\	/
		//		 \ /
		//		  N
		//	 \	 /
		//	  \ / <-- parental CLA of focal neighbor is returned
		//	   F
		//	  /
		//
		const InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getValidParentalCondLikePtr();
		}

	// focal neighbor N is a tip node (N equals the likelihood root in this case)
	//
	//		  N
	//	 \	 /
	//	  \ / <-- parental CLA of focal neighbor is returned
	//	   F
	//	  /
	//
	const TipData * child_tip_data = actual_child->GetTipData();
	PHYCAS_ASSERT(child_tip_data != NULL);
	return child_tip_data->getValidParentalCondLikePtr();
#else
	const TreeNode * c = edge.getActualChild();
	if (edge.getFocalNode() == c)
		{
		PHYCAS_ASSERT(c->IsInternal());
		const InternalData * childInternalData = c->GetInternalData();
		PHYCAS_ASSERT(childInternalData != NULL);
		return childInternalData->getValidChildCondLike();
		}
	// moving up the tree in calculations (root to leaves).
	PHYCAS_ASSERT(c == edge.getFocalNeighbor());
	if (c->IsInternal())
		{
		const InternalData * childInternalData = c->GetInternalData();
		PHYCAS_ASSERT(childInternalData != NULL);
		return childInternalData->getValidParentalCondLike();
		}
	const TipData * childTipData = c->GetTipData();
	PHYCAS_ASSERT(childTipData != NULL);
	return childTipData->getValidParentalCondLike();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to true.
*/
void TreeLikelihood::setNoData()
	{
	no_data = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to false.
*/
void TreeLikelihood::setHaveData()
	{
	no_data = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the `partition_model' data member, replacing the partition model defined in the
|   constructor.
*/
void TreeLikelihood::replacePartitionModel(
  PartitionModelShPtr m)
	{
	partition_model = m;
	recalcRelativeRates();

    // invalidate all CLAs
    useAsLikelihoodRoot(NULL);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the values of `subset_num_states' and `subset_num_rates' data members according to each `subset_model', then
|   calls the recalcRatesAndProbs function of the model to force recalculation of its `rate_means' and `rate_probs'
|   vectors. Should be called after changing the number of rate categories, the gamma shape parameter, or the pinvar
|   parameter of any subset model. Note that if the number of rate categories changes, trees on which likelihoods need
|   to be calculated also need to be re-equipped by calling prepareForLikelihood.
*/
void TreeLikelihood::recalcRelativeRates()
	{
    //std::cerr << "TreeLikelihood::recalcRelativeRates:" << std::endl;
	for (unsigned i = 0; i < partition_model->getNumSubsets(); ++i)
	    {
		PHYCAS_ASSERT(partition_model->subset_num_states[i] == partition_model->subset_model[i]->getNumStates());
		PHYCAS_ASSERT(partition_model->subset_num_rates[i] == partition_model->subset_model[i]->getNRatesTotal());
        partition_model->subset_model[i]->recalcRatesAndProbs(rate_means[i], rate_probs[i]); //POL_BOOKMARK recalcRatesAndProbs call

        //std::cerr << "  pinvar[" << i << "] = " << partition_model->subset_model[i]->getPinvar();
        //std::cerr << "\n  shape[" << i << "] = " << partition_model->subset_model[i]->getShape();
        //std::cerr << "\n  rate_means[" << i << "] = ";
        //std::copy(rate_means[i].begin(), rate_means[i].end(), std::ostream_iterator<double>(std::cerr, " "));
        //std::cerr << "\n  rate_probs[" << i << "] = ";
        //std::copy(rate_probs[i].begin(), rate_probs[i].end(), std::ostream_iterator<double>(std::cerr, " "));
        //std::cerr << "\n";

	    }
    //std::cerr << std::endl;
	if (!no_data)
		cla_pool->setCondLikeDimensions(partition_model->subset_num_patterns, partition_model->subset_num_rates, partition_model->subset_num_states);
	underflow_manager.setDimensions(partition_model->subset_num_patterns, partition_model->subset_num_rates, partition_model->subset_num_states);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `cla_pool' data member.
*/
const CondLikelihoodStorageShPtr TreeLikelihood::getCLAStorage() const
	{
	return cla_pool;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Utility function that returns the sum of pattern counts stored in the `pattern_counts' data member. Note that this
|	function casts the sum of pattern counts (which are doubles) to an unsigned value. This is only strictly appropriate
|	if the pattern counts are actually whole numbers, but in many (all?) cases where fractional pattern counts are used
|	e.g. Gelfand-Ghosh method), the sum of the pattern counts should be a whole number even though individual pattern
|	counts are not.
*/
unsigned TreeLikelihood::sumPatternCounts() const
	{
	return (unsigned)std::accumulate(pattern_counts.begin(), pattern_counts.end(), 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the length of the `pattern_counts' vector, which is equal to the number of stored
|	patterns.
*/
unsigned TreeLikelihood::getNumPatterns() const
	{
	return (unsigned)pattern_counts.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `_num_taxa' data member.
*/
unsigned TreeLikelihood::getNTaxa() const
	{
	return _num_taxa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `_num_char' data member.
*/
unsigned TreeLikelihood::getNChar() const
	{
	return _num_char;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the number of rates for subset i. Assumes that the length of the `rate_means' vector
|	for subset i equals the value of the `partition_model->subset_num_rates' data member for subset i. Also assumes that
|	the number of subsets in the partition is greater than i.
*/
unsigned TreeLikelihood::getNRatesTotal(
  unsigned i) 	/**< is the subset of interest */
  const
	{
	PHYCAS_ASSERT(partition_model->getNumSubsets() > i);
	PHYCAS_ASSERT(rate_means[i].size() == partition_model->subset_num_rates[i]);
	return partition_model->subset_num_rates[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the number of states for subset i.
*/
unsigned TreeLikelihood::getNumStates(
  unsigned i) 	/**< is the subset of interest */
  const
	{
	PHYCAS_ASSERT(partition_model->getNumSubsets() > i);
	return partition_model->subset_num_states[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns a copy of the (shared_ptr) data member `partition_model'.
*/
PartitionModelShPtr TreeLikelihood::getPartitionModel() const
	{
	return partition_model;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `state_list'.
*/
const state_list_vect_t & TreeLikelihood::getStateList() const
	{
	return state_list;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `state_list_pos'.
*/
const state_list_pos_vect_t & TreeLikelihood::getStateListPos() const
	{
	return state_list_pos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns `rate_means[i]'. Assumes that the number of subsets in the partition is greater than
|	i.
*/
const std::vector<double> & TreeLikelihood::getRateMeans(
  unsigned i) const	/**< is the subset for which rate means are desired */
	{
	PHYCAS_ASSERT(partition_model->getNumSubsets() > i);
	return rate_means[i]; //POL_BOOKMARK TreeLikelihood::getRateMeans
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns `rate_probs[i]'. Assumes that the number of subsets in the partition is greater than
|	i.
*/
const std::vector<double> & TreeLikelihood::getRateProbs(
  unsigned i) const	/**< is the subset for which rate probabilities are desired */
	{
	PHYCAS_ASSERT(partition_model->getNumSubsets() > i);
	return rate_probs[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a vector of lower rate category boundaries for subset i, where i=0,1,.... Assumes that the number of subsets
|   in the partition is at least i+1.
*/
std::vector<double> TreeLikelihood::getCategoryLowerBoundaries(
  unsigned i    /**< is the subset for which category lower boundaries are desired */
  ) const
	{
	PHYCAS_ASSERT(partition_model->getNumSubsets() > i);
	std::vector<double> tmp_means;
	std::vector<double> returned_boundaries;
	partition_model->subset_model[i]->recalcGammaRatesAndBoundaries(tmp_means, returned_boundaries);
	return returned_boundaries;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls TreeLikelihood::simulateImpl specifying that the transition probabilities should be recalculated (or
|	calculated the first time) before beginning simulations. After TreeLikelihood::simulateFirst is called once,
|	TreeLikelihood::simulate can be called many more times to generate more data sets using the same transition
|	probabilities.
*/
void TreeLikelihood::simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar)
	{
	simulateImpl(sim_data, t, rng, nchar, true);	// true means recalculate transition probabilities
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls TreeLikelihood::simulateImpl specifying that the transition probabilities should NOT be recalculated. Only
|	call this function after calling TreeLikelihood::simulateFirst at least once, otherwise the transition probabilities
|	will contain garbage.
*/
void TreeLikelihood::simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar)
	{
    // this function is obsolete, transition probs should always be recalculated at least until the consequences
    // of data partitioning are considered carefully
    PHYCAS_ASSERT(0);
	//simulateImpl(sim_data, t, rng, nchar, false);	// false means do not recalculate transition probabilities
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the (internal) node to use as the likelihood root (it will be stored in the `likelihood_root' data member).
|	The likelihood root is separate from the actual root of the tree, and specifies the node used when harvesting the
|	log-likelihood from surrounding conditional likelihood arrays (CLAs). It behooves one to set the likelihood root to
|	that node requiring the fewest CLA recalculations. Specifying NULL for the likelihood root will result in the
|	unconditional recalculation of all CLAs in the entire tree, and subsequently the likelihood root will be set to the
|	subroot node of the tree (the only child of the root).
*/
void TreeLikelihood::useAsLikelihoodRoot(
  TreeNode * nd)
	{
	likelihood_root = nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Updates the conditional likelihood array for `nd' in a direction away from `avoid'. All adjacent nodes (other than
|	`avoid') are assumed to be valid.
*/
void TreeLikelihood::refreshCLA(
  TreeNode & nd,            /**< is the focal node */
  const TreeNode * avoid)   /**< is the node to avoid */
	{
	if (nd.IsTip())
		return;

	//if (const_cast<TreeNode *>(avoid)->GetParent() == &nd)
	//	{
	//	std::cerr << "@@@@@@@@@@ avoid = " << avoid->GetNodeNumber() << std::endl;
	//	std::cerr << "@@@@@@@@@@ nd    = " << nd.GetNodeNumber() << std::endl;
	//	}

//printf("ENTERING refreshCLA with nd=%d avoid=%d\n", nd.GetNodeNumber(), (avoid==NULL) ? - 1 : avoid->GetNodeNumber());
    unsigned num_subsets = partition_model->getNumSubsets();

	PHYCAS_ASSERT(avoid != NULL);
	//@MTH-NESCENT this const_cast should be safe because we aren't relying on const-ness of CLA's but it needs to be revisited
	//@POL-NESCENT Mark, if we are going to immediately cast away the const on avoid, why do we need to make it const in the first place?
	CondLikelihoodShPtr ndCondLike = getCondLikePtr(&nd, const_cast<TreeNode *>(avoid));
	TreeNode * parent = nd.GetParent();
	TreeNode * lChild = nd.GetLeftChild();
	PHYCAS_ASSERT(parent != NULL);
	PHYCAS_ASSERT(lChild != NULL);

	// The first neighbor can either be the parent or the leftmost child. The second neighbor must be a child, but
	// which child depends on the first neighbor. Third and subsequent neighbors are always next sibs.
	TreeNode * firstNeighbor = NULL;
	TreeNode * secondNeighbor = NULL;
	double firstEdgeLen = 0.0;
	const bool movingTowardLeaves = !(parent == avoid);
	if (movingTowardLeaves)
		{
		// A child of nd is closer to the likelihood root than nd
		firstNeighbor = parent;
		firstEdgeLen = nd.GetEdgeLen();
		secondNeighbor = (lChild == avoid ? lChild->GetRightSib() : lChild);
		}
	else
		{
		// The parent of nd is closer to the likelihood root than nd
		firstNeighbor = lChild;
		firstEdgeLen = lChild->GetEdgeLen();
		secondNeighbor = lChild->GetRightSib();
		}

	PHYCAS_ASSERT(firstNeighbor != NULL);
	PHYCAS_ASSERT(secondNeighbor != NULL);

	if (firstNeighbor->IsTip())
		{
		TipData & firstTD = *(firstNeighbor->GetTipData());
        for (unsigned i = 0; i < num_subsets; ++i)
    		calcPMatTranspose(i, firstTD.getTransposedPMatrices(i), firstTD.getConstStateListPos(i), firstEdgeLen);
		if (secondNeighbor->IsTip())
			{
			// 1. both neighbors are tips
			TipData & secondTD = *(secondNeighbor->GetTipData());
            for (unsigned i = 0; i < num_subsets; ++i)
                calcPMatTranspose(i, secondTD.getTransposedPMatrices(i), secondTD.getConstStateListPos(i), secondNeighbor->GetEdgeLen());
			calcCLATwoTips(*ndCondLike, firstTD, secondTD);
			}
		else
			{
			// 2. first neighbor is a tip, but second is an internal node
			InternalData & secondID = *(secondNeighbor->GetInternalData());
            for (unsigned i = 0; i < num_subsets; ++i)
			    calcPMat(i, secondID.getPMatrices(i), secondNeighbor->GetEdgeLen());
			CondLikelihoodShPtr secCL = getCondLikePtr(secondNeighbor, &nd);
			calcCLAOneTip(*ndCondLike, firstTD, secondID, *secCL);
			}
		}
	else
		{
		InternalData & firstID = *(firstNeighbor->GetInternalData());
        for (unsigned i = 0; i < num_subsets; ++i)
    		calcPMat(i, firstID.getPMatrices(i), firstEdgeLen);
		const CondLikelihood & firCL = *getCondLikePtr(firstNeighbor, &nd);
		if (secondNeighbor->IsTip())
			{
			// 3. first neighbor internal node, but second is a tip
			TipData & secondTD = *(secondNeighbor->GetTipData());
            for (unsigned i = 0; i < num_subsets; ++i)
	    		calcPMatTranspose(i, secondTD.getTransposedPMatrices(i), secondTD.getConstStateListPos(i), secondNeighbor->GetEdgeLen());
			calcCLAOneTip(*ndCondLike, secondTD, firstID, firCL);
			}
		else
			{
			// 4. both neighbors are internal nodes
			InternalData & secondID = *(secondNeighbor->GetInternalData());
            for (unsigned i = 0; i < num_subsets; ++i)
    			calcPMat(i, secondID.getPMatrices(i), secondNeighbor->GetEdgeLen());
			const CondLikelihood & secCL = *getCondLikePtr(secondNeighbor, &nd);
			calcCLANoTips(*ndCondLike, firstID, firCL, secondID, secCL);
			}
		}

	// Deal with possible polytomy in which secondNeighbor has siblings
	for (TreeNode * currNd = secondNeighbor->GetRightSib(); currNd != NULL; currNd = currNd->GetRightSib())
		{
		if (currNd != avoid)
			{
			if (currNd->IsTip())
				{
				TipData & currTD = *(currNd->GetTipData());
                for (unsigned i = 0; i < num_subsets; ++i)
		    		calcPMatTranspose(i, currTD.getTransposedPMatrices(i), currTD.getConstStateListPos(i), currNd->GetEdgeLen());
				conditionOnAdditionalTip(*ndCondLike, currTD);
				}
			else
				{
				InternalData & currID = *(currNd->GetInternalData());
                for (unsigned i = 0; i < num_subsets; ++i)
		    		calcPMat(i, currID.getPMatrices(i), currNd->GetEdgeLen());
				const CondLikelihood & currCL = *getCondLikePtr(currNd, &nd);
				conditionOnAdditionalInternal(*ndCondLike, currID, currCL);
				}
			}
		}

#if 0
	// Turn this section on for debugging purposes only! Walks through conditional likelihood arrays
	// for this node looking for any conditional likelihood that is negative. Negative cond. likes
	// have shown up in the past (18 Oct 2007) for the GTR model, which results in a NaN (represented
	// as -1.#IND in the VC compiler when the log of the site likelihood is taken.
	LikeFltType * cla = ndCondLike->getCLA();
	for (unsigned i = 0; i < num_rates; ++i)
		{
		for (unsigned j = 0; j < num_patterns; ++j)
			{
			for (unsigned k = 0; k < num_states; ++k)
				{
				if (*cla++ < 0.0)
					{
					std::cerr << "*** Doh! cla negative for rate = " << i << ", pattern = " << j << ", state = " << k << " ***" << std::endl;
					std::exit(0);
					}
				}
			}
		}

#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the conditional likelihood of refNd is up to date for calculations centered at some effective root
|	node (neighborCloserToEffectiveRoot will be a node adjacent to refNd, but closer than refNd to the the effective
|	root). Called in a context in which neighborCloserToEffectiveRoot is requesting that all of its neighbors update
|	their likelihood temporaries.
*/
bool TreeLikelihood::isValid(const TreeNode * refNd, const TreeNode * neighborCloserToEffectiveRoot)
	{
	PHYCAS_ASSERT(refNd != NULL);
	PHYCAS_ASSERT(neighborCloserToEffectiveRoot != NULL);
	PHYCAS_ASSERT(refNd != neighborCloserToEffectiveRoot);
	if (refNd->IsTip())
		{
		// tip nodes always return true because they must be the child of neighborCloserToEffectiveRoot
		// and they do not hold filial CLAs
		return true;
		}
	else
		{
		// refNd is internal
		if (refNd->GetParentConst() == neighborCloserToEffectiveRoot)
			{
			// refNd is the child of neighborCloserToEffectiveRoot
			// If refNd has a filial CLA, then return true because that means refNd is valid
			const InternalData * id = refNd->GetInternalData();
			return (bool)(id->childWorkingCLA);
			}
		else
			{
			// neighborCloserToEffectiveRoot is the child of refNd
			// If neighborCloserToEffectiveRoot has a parental CLA, then return true because
			// that means refNd is valid
			const InternalData * id = neighborCloserToEffectiveRoot->GetInternalData();
			return (bool)(id->parWorkingCLA);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally invalidates parental (and, if ref_nd is internal, filial) conditional likelihood arrays and also
|	removes cached conditional likelihood arrays if present. This function always returns false so it can be used in
|	conjunction with effective_postorder_edge_iterator to invalidate every CLA in the entire tree and ensure that there
|	are also no cached CLAs as well.
*/
bool TreeLikelihood::invalidateBothEndsDiscardCache(TreeNode * ref_nd, TreeNode * /* unused */)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();
		if (td != NULL)
			{
			// Invalidate the parental CLAs if they exist
			if (td->parWorkingCLA)
				{
				cla_pool->putCondLikelihood(td->parWorkingCLA);
				td->parWorkingCLA.reset();
				}
			// Remove cached parental CLAs if they exist
			if (td->parCachedCLA)
				{
				cla_pool->putCondLikelihood(td->parCachedCLA);
				td->parCachedCLA.reset();
				}
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();
		if (id != NULL)
			{
			// Invalidate the parental CLAs if they exist
			if (id->parWorkingCLA)
				{
				cla_pool->putCondLikelihood(id->parWorkingCLA);
				id->parWorkingCLA.reset();
				}
			// Remove cached parental CLAs if they exist
			if (id->parCachedCLA)
				{
				cla_pool->putCondLikelihood(id->parCachedCLA);
				id->parCachedCLA.reset();
				}

			// Invalidate the filial CLAs if they exist
			if (id->childWorkingCLA)
				{
				cla_pool->putCondLikelihood(id->childWorkingCLA);
				id->childWorkingCLA.reset();
				}
			// Remove cached filial CLAs if they exist
			if (id->childCachedCLA)
				{
				cla_pool->putCondLikelihood(id->childCachedCLA);
				id->childCachedCLA.reset();
				}
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the data member `cla_pool'.
*/
CondLikelihoodStorageShPtr TreeLikelihood::getCondLikelihoodStorage()
	{
	return cla_pool;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally invalidate parental (and, if ref_nd is internal, filial) conditional likelihood arrays. Always
|	returns false so it can be used in conjunction with effective_postorder_edge_iterator to invalidate every CLA in
|	the entire tree.
*/
bool TreeLikelihood::invalidateBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Invalidate the parental CLAs if they exist
		if (td->parWorkingCLA)
			{
			if (td->parCachedCLA)
				cla_pool->putCondLikelihood(td->parCachedCLA);
			td->parCachedCLA = td->parWorkingCLA;
			td->parWorkingCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Invalidate the parental CLAs if they exist
		if (id->parWorkingCLA)
			{
			if (id->parCachedCLA)
				cla_pool->putCondLikelihood(id->parCachedCLA);
			id->parCachedCLA = id->parWorkingCLA;
			id->parWorkingCLA.reset();
			}

		// Invalidate the filial CLAs if they exist
		if (id->childWorkingCLA)
			{
			if (id->childCachedCLA)
				cla_pool->putCondLikelihood(id->childCachedCLA);
			id->childCachedCLA = id->childWorkingCLA;
			id->childWorkingCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function used in conjunction with effective_postorder_edge_iterator to invalidate the appropriate conditional
|	likelihood arrays (CLAs) of all nodes starting with a focal node. For nodes that are descendants of the focal node
|	(descendant means that a node can be found using only left child and right sib pointers starting from the focal
|	node), it is the parental CLAs that are invalidated. For a node that is an ancestor of the focal node, it is the
|	filial CLAs that are invalidated. For all other nodes (e.g. on independent lineages derived from an ancestral node),
|	it is the parental CLAs that are invalidated. This pattern of invalidation ensures that the likelihood will be
|	correctly computed using any node in the tree as the likelihood root. For any likelihood root n, it now appears as
|	though an edge length just on the other side of the focal node from n has been changed, necessitating the
|	recalculation of all CLAs starting from that point back toward the likelihood root. The return value indicates
|	whether or not the iterator should continue into the subtree defined by `refNd' (away from
|	`neighborCloserToEffectiveRoot'). This function always returns false because the goal is to invalidate every node
|	that needs to be invalidated.
*/
bool TreeLikelihood::invalidateNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip() && !ref_nd->IsTipRoot())
		{
		TipData * td = ref_nd->GetTipData();
		if (!td->parWorkingCLA)
			return false;
		if (td->parCachedCLA)
			cla_pool->putCondLikelihood(td->parCachedCLA);
		td->parCachedCLA = td->parWorkingCLA;
		td->parWorkingCLA.reset();
		}
	else
		{
		if (ref_nd->GetParent() == neighbor_closer_to_likelihood_root)
			{
			// ref_nd is the actual child
			InternalData * id = ref_nd->GetInternalData();
			if (!id->parWorkingCLA)
				return false;
			if (id->parCachedCLA)
				cla_pool->putCondLikelihood(id->parCachedCLA);
			id->parCachedCLA = id->parWorkingCLA;
			id->parWorkingCLA.reset();
			}
		else
			{
			// neighbor_closer_to_likelihood_root is the actual child
			if (neighbor_closer_to_likelihood_root->IsTip())
				return false;
			// Either ref_nd is not a tip, or it is the tip serving as the root node
			InternalData * id = neighbor_closer_to_likelihood_root->GetInternalData();
			if (!id->childWorkingCLA)
				return false;
			if (id->childCachedCLA)
				cla_pool->putCondLikelihood(id->childCachedCLA);
			id->childCachedCLA = id->childWorkingCLA;
			id->childWorkingCLA.reset();
			}
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Invalidates all conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node'. That is, the
|	parental CLA of a child of `focal_node' would be invalidated, as would be the filial CLA for the parent of the
|	`focal_node'. The invalidation propagates throughout the tree from `focal_node' so that regardless of the node
|	serving as the likelihood root, the correct CLAs will be recalculated the next time the log-likelihood is
|	recomputed. If the edge of `focal_node' has been changed, the function invalidateBothEnds(`focal_node') should also
|	be called to invalidate the one remaining CLA not invalidated by this function.
*/
void TreeLikelihood::invalidateAwayFromNode(
  TreeNode & focal_node)		/**< invalidation of conditional likelihood arrays will proceed outward from this node */
	{
    if (no_data)
        return;

	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateNode, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Unconditionally discards cached parental (and, if ref_nd is internal, filial) conditional likelihood arrays if
|	present. This function always returns false so it can be used in conjunction with effective_postorder_edge_iterator
|	to discard every cached CLA in the entire tree (something that needs to be done when any move is accepted).
*/
bool TreeLikelihood::discardCacheBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Remove cached parental CLAs if they exist
		if (td->parCachedCLA)
			{
			cla_pool->putCondLikelihood(td->parCachedCLA);
			td->parCachedCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Remove cached parental CLAs if they exist
		if (id->parCachedCLA)
			{
			cla_pool->putCondLikelihood(id->parCachedCLA);
			id->parCachedCLA.reset();
			}

		// Remove cached filial CLAs if they exist
		if (id->childCachedCLA)
			{
			cla_pool->putCondLikelihood(id->childCachedCLA);
			id->childCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores parental conditional likelihood arrays from cache and stores the existing conditional likelihood arrays
|	for later use. This function always returns false so it can be used in conjunction with
|	effective_postorder_edge_iterator to restore every parental CLA from cache in the entire tree. If `ref_nd'
|	has a working parental CLA but no corresponding cached parental CLA, then the working parental CLA is discarded.
|	This is done because this working parental CLA must have been calculated while the tree was in a proposed state,
|	and thus would be invalid if the tree were reverted.
*/
bool TreeLikelihood::restoreFromCacheParentalOnly(TreeNode * ref_nd, TreeNode * /* unused */)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Store working CLA in any case
		if (td->parWorkingCLA)
			cla_pool->putCondLikelihood(td->parWorkingCLA);
		td->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		InternalData * id = ref_nd->GetInternalData();

		// Store working CLA in any case
		if (id->parWorkingCLA)
			cla_pool->putCondLikelihood(id->parWorkingCLA);
		id->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->parCachedCLA)
			{
			id->parWorkingCLA = id->parCachedCLA;
			id->parCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores parental (and, if ref_nd is internal, filial) conditional likelihood arrays from cache and stores the
|	existing conditional likelihood arrays for later use. This function always returns false so it can be used in
|	conjunction with effective_postorder_edge_iterator to restore every CLA from cache in the entire tree. If `ref_nd'
|	has a working CLA but no corresponding cached CLA, then the working CLA is discarded. This is done because this
|	working CLA must have been calculated while the tree was in a proposed state, and thus would be invalid if the tree
|	were reverted.
*/
bool TreeLikelihood::restoreFromCacheBothEnds(TreeNode * ref_nd, TreeNode * /* unused */)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip())
		{
		// Tip nodes have only parental CLAs
		TipData * td = ref_nd->GetTipData();

		// Store working CLA in any case
		if (td->parWorkingCLA)
			cla_pool->putCondLikelihood(td->parWorkingCLA);
		td->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		// Internal nodes have both parental and filial CLAs
		InternalData * id = ref_nd->GetInternalData();

		// Restore parental CLA from cache

		// Store working CLA in any case
		if (id->parWorkingCLA)
			cla_pool->putCondLikelihood(id->parWorkingCLA);
		id->parWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->parCachedCLA)
			{
			id->parWorkingCLA = id->parCachedCLA;
			id->parCachedCLA.reset();
			}

		// Restore filial CLA from cache

		// Store working CLA in any case
		if (id->childWorkingCLA)
			cla_pool->putCondLikelihood(id->childWorkingCLA);
		id->childWorkingCLA.reset();

		// Move cached to working if cached exists
		if (id->childCachedCLA)
			{
			id->childWorkingCLA = id->childCachedCLA;
			id->childCachedCLA.reset();
			}
		}

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function used in conjunction with effective_postorder_edge_iterator to restore the conditional likelihood arrays
|	(CLAs) of all nodes starting with the supplied focal node `ref_nd'. For nodes that are descendants of the focal node
|	(descendant means that a node can be found using only left child and right sib pointers starting from the focal
|	node), it is the parental CLAs that are restored. For a node that is an ancestor of the focal node, it is the
|	filial CLAs that are restored. For all other nodes (e.g. on independent lineages derived from an ancestral node),
|	it is the parental CLAs that are restored. This function always returns false because the goal is to restore every
|	node that needs to be invalidated. If a node is has a working CLA but no corresponding cached CLA, then the working
|	CLA is discarded. This is done because this working CLA must have been calculated while the tree was in a proposed
|	state, and thus would be invalid if the tree were reverted.
*/
bool TreeLikelihood::restoreFromCacheNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root)
	{
    if (no_data)
        return false;

	if (ref_nd->IsTip() && !ref_nd->IsTipRoot())
		{
		TipData * td = ref_nd->GetTipData();

		// Store the working CLA in any case
		if (td->parWorkingCLA)
			cla_pool->putCondLikelihood(td->parWorkingCLA);
		td->parWorkingCLA.reset();

		// Move cached to working if there is a cached CLA
		if (td->parCachedCLA)
			{
			td->parWorkingCLA = td->parCachedCLA;
			td->parCachedCLA.reset();
			}
		}
	else
		{
		if (ref_nd->GetParent() == neighbor_closer_to_likelihood_root)
			{
			InternalData * id = ref_nd->GetInternalData();

			// Store the working CLA in any case
			if (id->parWorkingCLA)
				cla_pool->putCondLikelihood(id->parWorkingCLA);
			id->parWorkingCLA.reset();

			// Move cached to working if there is a cached CLA
			if (id->parCachedCLA)
				{
				id->parWorkingCLA = id->parCachedCLA;
				id->parCachedCLA.reset();
				}
			}
		else
			{
			// neighbor_closer_to_likelihood_root is the actual child
			if (neighbor_closer_to_likelihood_root->IsTip())
				return false;

			// Either ref_nd is not a tip, or it is the tip serving as the root node
			InternalData * id = neighbor_closer_to_likelihood_root->GetInternalData();

			// Store the working CLA in any case
			if (id->childWorkingCLA)
				cla_pool->putCondLikelihood(id->childWorkingCLA);
			id->childWorkingCLA.reset();

			// Move cached to working if there is a cached CLA
			if (id->childCachedCLA)
				{
				id->childWorkingCLA = id->childCachedCLA;
				id->childCachedCLA.reset();
				}
			}
		}
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores all conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node' from cache. That is, the
|	parental CLA of a child of `focal_node' would be restored, as would be the filial CLA for the parent of the
|	`focal_node'. It may be necessary to call the function invalidateBothEnds(`focal_node') also to restore the one
|	remaining CLA not restored by this function.
*/
void TreeLikelihood::restoreFromCacheAwayFromNode(
  TreeNode & focal_node)		/**< restoration of cached conditional likelihood arrays will proceed outward from this node */
	{
    if (no_data)
        return;

	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::restoreFromCacheNode, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Discards all cached conditional likelihood arrays (CLAs) pointing toward the supplied `focal_node'. That is, the
|	cached parental CLA of a child of `focal_node' would be discarded, as would be the cached filial CLA for the parent
|	of the `focal_node'. None of the working CLAs are affected. It may be necessary to call the function
|	discardCacheBothEnds(`focal_node') also to discard the cache from the two remaining CLAs not affected by this
|	function.
*/
void TreeLikelihood::discardCacheAwayFromNode(
  TreeNode & focal_node)		/**< restoration of cached conditional likelihood arrays will proceed outward from this node */
	{
    if (no_data)
        return;

	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::discardCacheBothEnds, this, _1, _2);
	effective_postorder_edge_iterator(&focal_node, validFunctor); // need only construct unnamed iterator object
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Simulates data for `nchar' characters using the current model and edge lengths. Transition matrices are recomputed
|	if `refresh_probs' is true; if `refresh_probs' is false, assumes transition probability matrices are up-to-date.
|	Only the `num_states' primary states are generated, and the state created for each node is stored initially in the
|	`state' data member of the TipData or InternalData structure associated with the tip or internal node, respectively.
|	This function expects these data structures to be in place (the TreeLikelihood::prepareForSimulation and
|	TreeLikelihood::prepareForLikelihood functions both perform this task). As states are generated for tip nodes, they
|	are copied into a temporary pattern inside `sim_data', and this pattern is then inserted into a pattern map inside
|	`sim_data' when it is completed. Assumes that `nchar' is greater than zero and that the shared pointers `sim_data'
|	and `rng' actually point to real objects.
|
|	The following makes use of T matrices, which are transposed, augmented transition matrices. These are transposed
|	because the "from" states form the columns rather than the rows. They are augmented because, ordinarily, there are
|	additional rows corresponding to ambiguities observed in some tip nodes. With simulated data, however, there are
|	never any ambiguities, so in this case T matrices are nothing more than transposed transition matrices.
|	Important: note that T matrices are only used in the TipData structures; InternalData structures store normal
|	untransposed transition probability matrices in which the rows form the "from" states.
|
|	Here is an example of a T matrix created using the JC model and an edge length equal to 0.1:
|>
|				|--------------- from state -------------|
|					0		   1		  2			 3
|					A		   C		  G			 T
|	t  0   A	 0.90638	0.03121	   0.03121	  0.03121
|	o  1   C	 0.03121	0.90638	   0.03121	  0.03121
|	   2   G	 0.03121	0.03121	   0.90638	  0.03121
|	s  3   T	 0.03121	0.03121	   0.03121	  0.90638
|	t  4   N	 1.00000	1.00000	   1.00000	  1.00000 \
|	a  5 {GT}	 0.06241	0.06241	   0.93759	  0.93759  | These rows not present if prepareForSimulation function
|	t  6 {ACT}	 0.96879	0.96879	   0.09362	  0.96879  | was used to create the TipData structures
|	e  7 {AG}	 0.93757	0.06241	   0.93759	  0.06241 /
|>
|	The `pMatrixTranspose' data member in TipData structures holds the array of T matrices (one T matrix for each
|	rate category).
*/
void TreeLikelihood::simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs) //POLSIM
    {
    // formerly DISABLED_UNTIL_SIMULATION_WORKING_WITH_PARTITIONING
	PHYCAS_ASSERT(sim_data);
	PHYCAS_ASSERT(rng);
	PHYCAS_ASSERT(nchar > 0);

    // Determine number of subsets in the partition.
    unsigned num_partition_subsets = partition_model->getNumSubsets();

    // Get assignments of models to sites
    uint_vect_t default_site_assignments;
    uint_vect_t & site_assignments = default_site_assignments;
    if (num_partition_subsets > 1)
        site_assignments = partition_model->getSiteAssignments();   // Note: nchar is ignored if more than 1 subset
    else
        default_site_assignments.assign(nchar, 0);

    // Make sure SimData::tmp_pattern is the correct length to store a pattern
    sim_data->resetPatternLength(t->GetNTips());    // warning: this calls SimData::clear() so don't move this line down

    // Make sure sim_data's sim_pattern_vect is long enough to save all simulated patterns
    nchar = (unsigned)site_assignments.size();
    sim_data->resizePatternVect(nchar);

    // Recalculate transition probabilities if requested
    if (refresh_probs)
        {
        preorder_iterator nd = t->begin();

        // First preorder node is the root node and represents a special case
        // Its transition matrices must be computed using the "subroot" node's edge length
        // The subroot node's transition matrices need not be calculated
        TipData & ndTD = *(nd->GetTipData());
        TreeNode * subroot = nd->GetLeftChild();
        PHYCAS_ASSERT(subroot);
        PHYCAS_ASSERT(!subroot->GetRightSib()); //@POL need to create a IsSubroot() member function for TreeNode
        for (unsigned s = 0; s < num_partition_subsets; ++s)
            {
            calcTMatForSim(s, ndTD, subroot->GetEdgeLen()); // note: calcTMatForSim correctly accounts for subset relative rate
            }
        ++nd;

        // Skip subroot node as its transition matrices are never used and thus do not need to be computed
        ++nd;

        // Process the remaining nodes in the tree
        for (; nd != t->end(); ++nd)
            {
            if (nd->IsTip())
                {
                TipData & ndTD = *(nd->GetTipData());
                for (unsigned s = 0; s < num_partition_subsets; ++s)
                    {
                    calcTMatForSim(s, ndTD, nd->GetEdgeLen()); // note: calcTMatForSim correctly accounts for subset relative rate
                    }
                }
            else
                {
                InternalData & ndID = *(nd->GetInternalData());
                for (unsigned s = 0; s < num_partition_subsets; ++s)
                    {
                    calcPMat(s, ndID.getPMatrices(s), nd->GetEdgeLen()); // note: calcTMatForSim correctly accounts for subset relative rate
                    }
                }
            }
        }

    // A temporary vector containing just one entry used to specify the sites (or in this case site)
    // to which the nascent pattern belongs. We will reuse tmplist for each site.
    uint_vect_t tmplist;
    tmplist.push_back(0);

    for (unsigned subset = 0; subset < num_partition_subsets; subset++)
        {
        //unsigned num_patterns = partition_model->subset_num_patterns[subset];  //POLSIM
        unsigned num_states = partition_model->subset_num_states[subset];  //POLSIM
        unsigned num_rates = partition_model->subset_num_rates[subset];  //POLSIM

        // Create a vector of cumulative state frequencies to use in choosing starting states
        const std::vector<double> & freqs = partition_model->subset_model[subset]->getStateFreqs();  //POLSIM
        std::vector<double> cum_freqs(num_states, 0.0);
        std::partial_sum(freqs.begin(), freqs.end(), cum_freqs.begin());

        // Create a vector of cumulative rate probabilities to use in choosing relative rates
        std::vector<double> cum_rate_probs(num_rates, 0.0);
        std::partial_sum(rate_probs[subset].begin(), rate_probs[subset].end(), cum_rate_probs.begin()); //POLSIM

        sim_data->wipePattern();

        unsigned character = 0;
        for (std::vector<unsigned>::iterator site = site_assignments.begin(); site != site_assignments.end(); ++site, ++character)
            {
            if (*site != subset)
                continue;

            // Choose a rate for this character (actually, choose index, the actual rate is rate_means[r])
            unsigned r = 0;
            if (num_rates > 1)
                {
                // warning: removing the if statement will invalidate all examples involving simulated data with rate
                // homogeneity because of the call to rng->Uniform here!
                r = (unsigned)(std::lower_bound(cum_rate_probs.begin(), cum_rate_probs.end(), rng->Uniform()) - cum_rate_probs.begin());
                }

            // Generate the starting state
            int8_t j = (unsigned)(std::lower_bound(cum_freqs.begin(), cum_freqs.end(), rng->Uniform()) - cum_freqs.begin());

            // Assign starting state to the tip node currently serving as the root of the tree
            preorder_iterator nd = t->begin();
            TipData & rootTD = *(nd->GetTipData());
            rootTD.state = j;

            sim_data->setState(nd->GetNodeNumber(), j);

            // Go ahead and generate the state for the (only) descendant of the root node (the "subroot" node)
            // Note that the root node's T matrix is used for this calculation; the P matrix of the subroot node
            // is never computed
            unsigned parent_state = (unsigned)j;

            // Get the T matrix for the tip node serving as the root
            //double * * Tmatrix = rootTD.pMatrixTranspose[r];
            double * * Tmatrix = rootTD.getTransposedPMatrices(subset)[r];   //POLSIM

            // Choose a uniform random deviate
            double u = rng->Uniform();

            // Spin the roulette wheel to choose a state for the subroot node
            double cum = 0.0;
            unsigned i = 0;
            for (; i < num_states; ++i)
                {
                double pr = Tmatrix[i][parent_state];
                //std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
                cum += pr;
                if (u < cum)
                    break;
                }

            // Increment iterator so that nd now refers to the subroot (sole descendant of the root)
            ++nd;

            // Assign the new state to the subroot node
            InternalData & ndID = *(nd->GetInternalData());
            ndID.state = (int8_t)i;
            //std::cerr << "  Assigning state " << i << " to node " << nd->GetNodeNumber() << std::endl;

            // Walk the remainder of the tree using the preorder sequence, generating data for each node along the way
            for (++nd; nd != t->end(); ++nd)
                {
                // Get state of parent of nd
                TreeNode * parent = nd->GetParent();
                parent_state = UINT_MAX;
                if (parent->IsTip())
                    {
                    TipData * parentTD = parent->GetTipData();
                    parent_state = (unsigned)parentTD->state;
                    }
                else
                    {
                    InternalData * parentID = parent->GetInternalData();
                    parent_state = (unsigned)parentID->state;
                    }
                PHYCAS_ASSERT(parent_state < num_states);

                if (nd->IsTip())
                    {
                    // Get the T matrix
                    TipData & ndTD = *(nd->GetTipData());
                    //double * * Tmatrix = ndTD.pMatrixTranspose[r];
                    double * * Tmatrix = rootTD.getTransposedPMatrices(subset)[r];   //POLSIM

                    // Choose a uniform random deviate
                    double u = rng->Uniform();

                    // Spin the roulette wheel and assign a state to nd
                    double cum = 0.0;
                    unsigned i = 0;
                    for (; i < num_states; ++i)
                        {
                        double pr = Tmatrix[i][parent_state];
                        //std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
                        cum += pr;
                        if (u < cum)
                            break;
                        }
                    ndTD.state = (int8_t)i;
                    sim_data->setState(nd->GetNodeNumber(), (int8_t)i);
                    }
                else
                    {
                    // Get the T matrix
                    InternalData & ndID = *(nd->GetInternalData());
                    //double * * Pmatrix = ndID.pMatrices[r];
                    double * * Pmatrix = ndID.getPMatrices(subset)[r]; //POLSIM

                    // Choose a uniform random deviate
                    double u = rng->Uniform();

                    // Spin the roulette wheel and assign a state to nd
                    double cum = 0.0;
                    unsigned i = 0;
                    for (; i < num_states; ++i)
                        {
                        double pr = Pmatrix[parent_state][i];
                        //std::cerr << str(boost::format("Pmatrix[%d][%d] = %f") % parent_state % i % pr) << std::endl;
                        cum += pr;
                        if (u < cum)
                            break;
                        }
                    ndID.state = (int8_t)i;
                    }
                }

            // We are now finished simulating data for one character, so insert the pattern just generated
            // into the pattern map maintained by sim_data
            // tmplist is a temporary list containing a single entry (the index of this site)
            // The 1.0 means that the count for this pattern should be incremented by 1
            tmplist[0] = character;
            sim_data->insertPattern(tmplist, 1.0);
            }
        }   // subset loop
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores all CondLikelihood objects stored in TipData or InternalData structures. This will require all CLAs to be
|	recomputed the next time the likelihood needs to be computed. Returns the subroot node (only child of root node).
*/
TreeNode * TreeLikelihood::storeAllCLAs(
  TreeShPtr t)				/**< is the tree from which all CLAs are to be removed */
	{
    // temporary!
    //std::cerr << t->DebugWalkTree() << std::endl;
	//bool ok = t->DebugCheckTree(false, true, 1);

	// Start at the root node
	TreeNode * nd = t->GetFirstPreorder();
	PHYCAS_ASSERT(nd);

	// Move to the subroot node
	nd = nd->GetNextPreorder();
	PHYCAS_ASSERT(nd);

    if (no_data)
        return nd;

	// Invalidate (and do not cache) all CLAs from the tree. This will require all CLAs to be recomputed
	// when the likelihood is computed using the subroot as the likelihood root. This path should be taken
	// if a parameter is changed that invalidates the entire tree.
	NodeValidityChecker validFunctor = boost::bind(&TreeLikelihood::invalidateBothEndsDiscardCache, this, _1, _2);
	effective_postorder_edge_iterator(nd, validFunctor); // constructor does all the work we need
	invalidateBothEndsDiscardCache(nd);

	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sanity check to make sure TreeLikelihood::storeAllCLAs really has stored all CLAs. Returns true if any CLAs are
|	found, false if no CLAs are found (not even in cache positions).
*/
bool TreeLikelihood::debugCheckCLAsRemainInTree(
  TreeShPtr t) const	/**< is the tree to check */
	{
	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			TipData * td = nd->GetTipData();
			if (td != NULL)
				{
				if (td->parWorkingCLA)
					return true;
				if (td->parCachedCLA)
					return true;
				}
			}
		else
			{
			InternalData * id = nd->GetInternalData();
			if (id != NULL)
				{
				if (id->parWorkingCLA)
					return true;
				if (id->parCachedCLA)
					return true;

				if (id->childWorkingCLA)
					return true;
				if (id->childCachedCLA)
					return true;
				}
			}
		}	// preorder loop over all nodes

	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scans tree for cached CLAs. If any are found, returns true. If no CLAs are currently cached, returns false. False is
|	the expected result, because there should only be cached CLAs in the tree if we are in the middle of performing a
|	Metropolis-Hastings move. If a move has been proposed, but not yet accepted or rejected, cached CLAs provide a way
|	to return to the pre-move state after rejection without having to recalculate the likelihood.
*/
bool TreeLikelihood::debugCheckForUncachedCLAs(
  TreeShPtr t)	 /**< is the tree to check */
  const
	{
	bool cached_CLAs_found = false;
	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			TipData * td = nd->GetTipData();
			if (td != NULL)
				{
				if (td->parCachedCLA)
					{
					cached_CLAs_found = true;
					break;
					}
				}
			}
		else
			{
			InternalData * id = nd->GetInternalData();
			if (id != NULL)
				{
				if (id->parCachedCLA)
					{
					cached_CLAs_found = true;
					break;
					}

				if (id->childCachedCLA)
					{
					cached_CLAs_found = true;
					break;
					}
				}
			}
		}	// preorder loop over all nodes
	return cached_CLAs_found;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scans node for CLAs. If any are found, returns true. If no CLAs are found, returns false.
*/
bool TreeLikelihood::debugCheckCLAsRemainInNode(
  TreeNode * nd)   /**< is the node to check */
  const
	{
	bool CLAs_found = false;
	if (nd->IsTip())
		{
		TipData * td = nd->GetTipData();
		if (td != NULL)
			{
			if (td->parWorkingCLA || td->parCachedCLA)
				CLAs_found = true;
			}
		}
	else
		{
		InternalData * id = nd->GetInternalData();
		if (id != NULL)
			{
			if (id->parWorkingCLA || id->parCachedCLA || id->childWorkingCLA || id->childCachedCLA)
				CLAs_found = true;
			}
		}
	return CLAs_found;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Computes the value of the log-likelihood at substitutional saturation (this occurs when all edges have infinite
|   length). The site log-likelihood in this case is simply the sum of the logs of the equilibrium frequency of each tip
|   state: i.e.,
|>
|   log-like for site j = \sum_{i=1}^{ntax} log(\pi_{s_{ij}})
|>
|   Meaningless if no data has been stored, so assumes that `pattern_vect' is not empty.
*/
double TreeLikelihood::calcLogLikeAtSubstitutionSaturation() const
    {
    double lnL = 0.0;
#if DISABLED_UNTIL_WORKING_WITH_PARTITIONING
    //@POL needs to take ambiguity into account!
    PHYCAS_ASSERT(!no_data);
    PHYCAS_ASSERT(!pattern_vect.empty());

    unsigned nss = partition_model->getNumSubsets();
    for (unsigned i = 0; i < nss; ++i)
        {
		// Create a vector containing the log of each state frequency
		unsigned nstates = getNumStates(i);
		const double * stateFreq = &partition_model->subset_model[i]->getStateFreqs()[0]; //PELIGROSO
		double_vect_t logfreq(nstates, 0.0);
		for (unsigned k = 0; k < nstates; ++k)
			logfreq[k] = log(stateFreq[k]);

		for (unsigned j = subset_offset[i]; j < subset_offset[i + 1]; ++j)
			{
			const pattern_vect_t & 	pattern 		= pattern_vect[j];
			const pattern_count_t & count 			= pattern_counts[j];
			double 					site_log_like 	= 0.0;
			pattern_t::const_iterator sit = pattern.begin();
			// skip first element in pattern because it simply indicates the partition subset to which the patterns belongs
			for (++sit; sit != pattern.end(); ++sit)
				{
				state_code_t s = (*sit);
				PHYCAS_ASSERT(s >= (int8_t)0);
				PHYCAS_ASSERT(s < (int8_t)nstates);
				site_log_like += logfreq[s];
				}
			lnL += count*site_log_like;
			}
        }
#endif
	return lnL;
    }

int calcLnLLevel = 0;
/*----------------------------------------------------------------------------------------------------------------------
|	This is the function that needs to be called to recompute the log-likelihood. If `likelihood_root' is not NULL,
|	the calculation will use that node as the root of the likelihood calculation, and it will be assumed that all
|	conditional likelihood arrays (CLAs) are correctly calculated or have been invalidated if they need to be
|	recomputed. On the other hand, if `likelihood_root' is NULL, then all CLAs will be invalidated and recomputed
|	(useful if a parameter has been changed that requires recalculation of all CLAs).
*/
double TreeLikelihood::calcLnL(
  TreeShPtr t)
	{
	if (no_data)
		return 0.0;

//@TEMP force crash to test entry into debugger
//char*p=NULL;
//*p='a';

	// The variable nevals keeps track of the number of times the likelihood has been calculated
	// You can reset this value to 0 using resetNumLikelihoodEvals()
	incrementNumLikelihoodEvals();

	// If likelihood_root has not been specified, set it to the subroot node and invalidate
	// all CLAs. If likelihood_root does already point to a node, assume that the necessary
	// CLA invalidations have already been performed.
	TreeNode * nd = likelihood_root;
	if (nd == NULL)
		{
		// If no likelihood_root has been specified, invalidate the entire tree to be safe
		nd = storeAllCLAs(t);

		// The subroot node will be the new likelihood_root
		likelihood_root = nd;
		}

	//if (0)
	//	{
	//	likelihood_root = nd->GetLeftChild()->GetRightSib();
	//	nd = likelihood_root;
	//	}
	//nd->SelectNode();
	//startTreeViewer(t, "in TreeLikelihood::calcLnL");
	//nd->UnselectNode();

	PHYCAS_ASSERT(nd);
	PHYCAS_ASSERT(nd->IsInternal());

	// Calculate log-likelihood using nd as the likelihood root
	double lnL = calcLnLFromNode(*nd, t);

	// 	if (calcLnLLevel == 0)
	// 		{
	// 		calcLnLLevel = 1;
	// 		storeAllCLAs(t);
	// 		double lnLRecalc = calcLnL(t);
	// 		PHYCAS_ASSERT(fabs(lnL-lnLRecalc) < 0.000001);
	// 		calcLnLLevel = 0;
	// 		}

    //startTreeViewer(t, "lnL = %.5f" % lnL);

	return lnL;
	}

void TreeLikelihood::debugSaveCLAs(TreeShPtr t, std::string fn, bool overwrite)
	{
	std::ofstream tmpf;
	if (overwrite)
		tmpf.open(fn.c_str());
	else
		{
		tmpf.open(fn.c_str(), std::ios::out | std::ios::app);
		}

	for (TreeNode * nd = t->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			tmpf << "\n\nTip node number " << nd->GetNodeNumber();
			if (nd->IsTipRoot())
				tmpf << "\n	 parent = None";
			else
				tmpf << "\n	 parent = " << nd->GetParent()->GetNodeNumber();
			tmpf << "\n	 children = ";
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
				tmpf << child->GetNodeNumber() << " ";

			TipData * td = nd->GetTipData();
			if (td->parentalCLAValid())
				{
				CondLikelihoodShPtr cla = td->getParentalCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
				tmpf << "\n	 cla length = " << sz;
				tmpf << "\n	 parental = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
				tmpf << "\n	 parental CLA = None";
				}
			}
		else
			{
			tmpf << "\n\nInternal node number " << nd->GetNodeNumber();

			tmpf << "\n	 parent = " << nd->GetParent()->GetNodeNumber();
			tmpf << "\n	 children = ";
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
				tmpf << child->GetNodeNumber() << " ";

			InternalData * id = nd->GetInternalData();
			if (id->parentalCLAValid())
				{
				CondLikelihoodShPtr cla = id->getParentalCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
				tmpf << "\n	 parental CLA length = " << sz;
				tmpf << "\n	 parental CLA = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
				tmpf << "\n	 parental CLA = None";
				}

			if (id->filialCLAValid())
				{
				CondLikelihoodShPtr cla = id->getChildCondLikePtr();
				unsigned sz = cla->getCLASize();
				LikeFltType * arr = cla->getCLA();
				tmpf << "\n	 filial CLA length = " << sz;
				tmpf << "\n	 filial CLA = ";
				for (unsigned i = 0; i < sz; ++i)
					{
					tmpf << str(boost::format("%20.12f ") % arr[i]);
					}
				}
			else
				{
				tmpf << "\n	 filial CLA = None";
				}
			}
		}

	tmpf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log-likelihood using `focal_node' as the likelihood root (i.e. the root of the likelihood calculation,
|	but not necessarily the root of the tree).
*/
double TreeLikelihood::calcLnLFromNode(
  TreeNode & focal_node,	/**< is the likelihood root (i.e. node around which the likelihood will be computed) */
  TreeShPtr t)				/**< is the tree */
	{
	double lnL;
//printf("in calcLnLFromNode with focal_node=%d\n", focal_node.GetNodeNumber());//

	if (no_data)
		lnL =  0.0;
	else
		{
		PHYCAS_ASSERT(!focal_node.IsTip());

		// valid_functor will return true if the conditional likelihood arrays pointing away from the
		// focal_node are up-to-date, false if they need to be recomputed. _1 is the focal node,
		// and _2 is the "avoid" node, the node neighboring the focal node on the path to the likelihood root
		NodeValidityChecker valid_functor = boost::bind(&TreeLikelihood::isValid, this, _1, _2);

		// iter will visit nodes that need their CLAs updated centripetally (like a postorder traversal
		// but also coming from below the focal node). Each node visited is guaranteed by valid_functor
		// to need its CLA updated.
		effective_postorder_edge_iterator iter(&focal_node, valid_functor);
		effective_postorder_edge_iterator iter_end;
		for (; iter != iter_end; ++iter)
			{
			// first is the focal node, second is the avoid node
			//if (iter->second->GetParent() == iter->first)
			//	startTreeViewer(t, boost::str(boost::format("avoid = %d, nd = %d") % iter->second->GetNodeNumber() % iter->first->GetNodeNumber()));
			refreshCLA(*iter->first, iter->second);
			}

		// We have now brought all neighboring CLAs up-to-date, so we can now call harvestLnL to
		// compute the likelihood
		EdgeEndpoints edge(&focal_node, NULL);
		lnL = harvestLnL(edge, t);
        //debugWalkTreeShowCondLikes(t);
		}
	return lnL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a postorder traversal of the tree, printing out conditional likelihood arrays for every node encountered
|   along the way. Makes most sense to invoke this function only for very small datasets.
*/
void TreeLikelihood::debugWalkTreeShowCondLikes(    //POLCURR
  TreeShPtr t)		/**< is the tree to walk */
	{
    for (TreeNode *nd = t->GetLastPreorder(); nd; nd = nd->GetNextPostorder())
        {
            if (nd->IsTip()) {
                std::cerr << "Tip node:" << nd->GetNodeNumber() << std::endl;
            } else {
                std::cerr << "Internal node:" << nd->GetNodeNumber() << std::endl;
                std::cerr << "-->|";
                const InternalData * internal_data = nd->GetInternalData();
                PHYCAS_ASSERT(internal_data != NULL);
                ConstCondLikelihoodShPtr condlike = internal_data->getValidChildCondLikePtr();
                const LikeFltType * cla = condlike->getCLA();
                unsigned num_subsets = partition_model->getNumSubsets();
                for (unsigned i = 0; i < num_subsets; ++i) {
                    unsigned num_patterns = partition_model->subset_num_patterns[i];
                    unsigned num_states = partition_model->subset_num_states[i];
                    unsigned num_rates = partition_model->subset_num_rates[i];
                    for (unsigned r = 0; r < num_rates; ++r) {
                        for (unsigned p = 0; p < num_patterns; ++p, cla += num_states) {
                            for (unsigned s = 0; s < num_states; ++s) {
                                std::cerr << boost::str(boost::format("%d:%g") % s % cla[s]) << '|';
                            }
                        }
                    }
                }
                std::cerr << std::endl;
            }
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
 |	Allocates the TipData data structure needed to store the data for one tip (the tip corresponding to the supplied
 |	`row' index in the data matrix `mat'). Returns a pointer to the newly-created TipData structure. See documentation
 |	for the TipData structure for more explanation.
 */
TipData * TreeLikelihood::allocateTipData(	//POLBM TreeLikelihood::allocateTipData
  unsigned row)		/**< is the row of the data matrix corresponding to the data for this tip node */
    {
    // The first element of a pattern vector is used to store the index of the partition subset to which the pattern belongs
	// Hence, add 1 to the requested row to get the correct element of the pattern vector
	unsigned index_of_taxon = row + 1;

	unsigned nsubsets = partition_model->getNumSubsets();

	std::vector< std::map<int8_t, int8_t> >		globalToLocal(nsubsets);
	std::map<int8_t, int8_t>::const_iterator	foundElement;
	uint_vect_t									nPartialAmbig(nsubsets, 0);
	state_list_pos_vect_t						local_state_list_pos(nsubsets);
	state_list_vect_t							local_state_codes(nsubsets);
	for (unsigned i = 0; i < nsubsets; ++i)
		{
		// these capacities are larger than they need to be
		local_state_list_pos[i].reserve(state_list_pos[i].size());
		local_state_codes[i].reserve(state_list[i].size());
		}
		// Example:
		// Here is the observed data for taxon k. The partition comprises 2 subsets: the first subset
		// includes sites 1-4 and consists of DNA data. The second subset includes sites 5-9 and is
		// amino acid data.
		//           1  2  3  4  5  6   7   8  9
		//  taxon_k  A  C  T  R  G  V  (EQ) ?  A
		//
		//  DNA        ------ amino acid -----
		//  0 A        0 A   5 E   10 L   15 S
		//  1 C        1 R   6 Q   11 K   16 T
		//  2 G        2 N   7 G   12 M   17 W
		//  3 T        3 D   8 H   13 F   18 Y
		//  4 (ACGT)   4 C   9 I   14 P   19 V
		//  5 (AG)    20 (ARNDCEQGHILKMFPSTWYV)
		//            21 (EQ)
		//
		//  If the only ambiguities in the entire data matrix were those found in taxon k,
		//                  |-------------- subset 1 ------------|
		//                  | -A- -C- -G- -T- -----?------ --R-- |
		//  state_list[0]:  | 1 0 1 1 1 2 1 3 5 -1 0 1 2 3 2 0 2 |
		//                    ^   ^   ^   ^   ^            ^
		//  state_list_pos:   0   2   4   6   8           14
		//           index:   0   1   2   3   4            5

		//                  |------------------------------------------------------------------------------------- subset 2 ----------------------------------------------------------|
		//                  | -A- -R- -N- -D- -C- -E- -Q- -G- -H- -I- -L-- -K-- -M-- -F-- -P-- -S-- -T-- -W-- -Y-- -V-- --------------------------?---------------------------- -(EQ)-
		//  state_list[1]:  | 1 0 1 1 1 2 1 3 1 4 1 5 1 6 1 7 1 8 1 9 1 10 1 11 1 12 1 13 1 14 1 15 1 16 1 17 1 18 1 19 21 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 2 5 6
		//                    ^   ^   ^   ^   ^   ^   ^   ^   ^   ^   ^    ^    ^    ^    ^    ^    ^    ^    ^    ^     ^                                                      ^
		//  state_list_pos:   0   2   4   6   8  10  12  14  16  18  20   22   24   26   28   30   32   34   36   38    40                                                     62
		//           index:   0   1   2   3   4   5   6   7   8   9  10   11   12   13   14   15   16   17   18   19    20                                                     21
		//
		for (pattern_vect_t::const_iterator it = pattern_vect.begin(); it != pattern_vect.end(); ++it)
			{
			// get subset-specific info
			const unsigned 	subset 		= (unsigned)(*it)[0];
			const int8_t 	ns			= partition_model->subset_num_states[subset];
			const int8_t	nsPlusOne	= ns + 1;

			const int8_t globalStateCode = (*it)[index_of_taxon];

			if (globalStateCode < nsPlusOne)
				{
				// no partial ambiguity, but may be gap state
				state_code_t s = (globalStateCode < 0 ? ns : globalStateCode);
				local_state_codes[subset].push_back(s);
				}
			else
				{
				// partial ambiguity
				foundElement = globalToLocal[subset].find(globalStateCode);
				if (foundElement == globalToLocal[subset].end())
					{
					// state code needs to be added to globalToLocal[subset] map
					state_code_t s = nPartialAmbig[subset] + nsPlusOne;
					globalToLocal[subset][globalStateCode] = s;
					local_state_list_pos[subset].push_back(state_list_pos[subset][globalStateCode]);
					local_state_codes[subset].push_back(s);
					nPartialAmbig[subset]++;
					}
				else
					{
					// state code is already in the globalToLocal[subset] map
					local_state_codes[subset].push_back(foundElement->second);
					}
				}
			}
#	if 0
		}
#	endif


	return new TipData( partition_model,
						local_state_list_pos,
						local_state_codes,
						cla_pool);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the InternalData data structure needed to store the conditional likelihood arrays and transition matrices
|	needed for likelihood calcuations. Returns a pointer to a newly-constructed InternalData structure. See the
|	documentation for InternalData for more explanation.
*/
InternalData * TreeLikelihood::allocateInternalData()
	{
	return new InternalData(partition_model,
							cla_pool);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the TipData structure allocated by the allocateTipData() function. The intention is for this
|	function to be given to TreeNode objects (in the form of a boost::function object) whenever allocateTipData() is
|	used to allocate memory for a TipData structure. The TreeNode destructor can then use the boost::function object to
|	delete the memory required by the TipData structure without knowing any details of the TipData structure (i.e.,
|	the file that defines the TreeNode destructor does not need to include a header file that declares the details of
|	the TipData class.
*/
void deallocateTipData(TipData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the InternalData structure allocated by the allocateInternalData() function. The intention is
|	for this function to be given to TreeNode objects (in the form of a boost::function object) whenever
|	allocateInternalData() is used to allocate memory for an InternalData structure. The TreeNode destructor can then
|	use the boost::function object to delete the memory required by the InternalData structure without knowing any
|	details of the InternalData structure (i.e., the file that defines the TreeNode destructor does not need to include
|	a header file that declares the details of the InternalData class.
*/
void deallocateInternalData(InternalData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a basic TipData object for all tip nodes and calls allocateInternalData() for all internal nodes in the
|	supplied Tree. This function does not require a data matrix and does not add observed data to TipData structures.
|	Silently returns if root node already has a TipData structure. This is taken to indicate that prepareForLikelihood
|	or prepareForSimulation was previously called for this tree, in which case all data structures needed for
|	simulation are already present.
*/
void TreeLikelihood::prepareForSimulation(
  TreeShPtr t)			/**< is the tree to decorate */
	{
    // formerly DISABLED_UNTIL_SIMULATION_WORKING_WITH_PARTITIONING
	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	preorder_iterator nd = t->begin();

	// Only proceed if root node does not already have a TipData structure
	if (nd->GetTipData() != NULL)
		return;

	for (; nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			TipData * td =	new TipData(partition_model, cla_pool);
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_vect', `pattern_counts', `pattern_to_sites' and `charIndexToPatternIndex' using uncompressed data
|	stored in `mat'.
*/
void TreeLikelihood::copyDataFromDiscreteMatrix(
  const CharSuperMatrix * mat,				/**< is the data source */
  const std::vector<unsigned> & partition_info)	/**< is a vector of indices storing the partition subset used by each site */

	{
	_num_taxa = mat->getNTax();
	_num_char = mat->getNChar();

	// Currently, can only deal with the first matrix stored in the CharSuperMatrix object. The CharSuperMatrix object
	// will contain multiple matrices if the nexus file contains a mixed datatype data block
	NxsCXXDiscreteMatrix & singleMat = *(mat->GetMatrix(0));

    // The compressDataMatrix function first erases, then builds, both pattern_vect and pattern_counts using the uncompressed data contained in mat
	//@POL commented out: compressDataMatrix(*mat->GetMatrix(0), partition_info);
	compressDataMatrix(singleMat, partition_info);

	unsigned nsubsets = partition_model->getNumSubsets();

	state_list.clear();
	state_list.resize(nsubsets);

	state_list_pos.clear();
	state_list_pos.resize(nsubsets);

	for (unsigned i = 0; i < nsubsets; ++i)
		{
		state_list[i].clear();
		state_list_pos[i].clear();
		if (partition_model->subset_model[i]->isCodonModel())
			{
			// any ambiguity is complete ambiguity for codon models, which makes it simple to construct
			// the state_list and state_list_pos vectors. Here is a reminder of how state_list and
			// state_list_pos are laid out for codon models:
			//
			// states          |   AAA |   AAC |   AAG |   AAT |   ACA | ... |   TTT | ambiguous
			// state_list      | 1  0  | 1  1  | 1  2  | 1  3  | 1  4  | ... | 1  60 | 61  0  1  2  3  4  5 ... 60
			// state_list_pos  | 0     | 2     | 4     | 6     | 8     | ... | 120   | 122
			state_list[i].reserve(184);
			state_list_pos[i].reserve(62);
			for (unsigned k = 0; k < 61; ++k)
				{
				state_list[i].push_back((int8_t)1);
				state_list[i].push_back((int8_t)k);
				state_list_pos[i].push_back((int8_t)(2*k));
				}
			state_list[i].push_back((int8_t)61);
			state_list_pos[i].push_back((int8_t)(2*61));
			for (unsigned k = 0; k < 61; ++k)
				{
				state_list[i].push_back((int8_t)k);
				}
			}
		else	// not codon model
			{
			//@POL commented out: NxsCXXDiscreteMatrix & singleMat = *(mat->GetMatrix(i));
			if (nsubsets == 1)
				{
				// get state_list and state_list_pos from mat
				state_list[i] = singleMat.getStateList();
				state_list_pos[i] = singleMat.getStateListPos();
				}
			else
				{
				//@POL Kansas2010 TODO
				// The user has specified partitioning, so we must build state_list and state_list_pos ourselves
				// Hopefully this functionality will be returned to NCL eventually, making this section unnecessary.

				//@POL assuming DNA sequence data here
				PHYCAS_ASSERT(partition_model->subset_num_states[i] == 4);

				//@POL assuming any ambiguity is complete ambiguity (needs to be revised)
				state_list[i].reserve(14);
				state_list_pos[i].reserve(5);

				state_list_pos[i].push_back((int8_t)0);
				state_list[i].push_back((int8_t)1); // A
				state_list[i].push_back((int8_t)0);

				state_list_pos[i].push_back((int8_t)2);
				state_list[i].push_back((int8_t)1); // C
				state_list[i].push_back((int8_t)1);

				state_list_pos[i].push_back((int8_t)4);
				state_list[i].push_back((int8_t)1); // G
				state_list[i].push_back((int8_t)2);

				state_list_pos[i].push_back((int8_t)6);
				state_list[i].push_back((int8_t)1); // T
				state_list[i].push_back((int8_t)3);

				state_list_pos[i].push_back((int8_t)8);
				state_list[i].push_back((int8_t)5); // ?
				state_list[i].push_back((int8_t)-1);
				state_list[i].push_back((int8_t)0);
				state_list[i].push_back((int8_t)1);
				state_list[i].push_back((int8_t)2);
				state_list[i].push_back((int8_t)3);
				}
			}	// not codon model
		}	// loop over subsets

	// The constant states vector stores information about which sites can potentially be
	// constant. A site can be potentially constant for more than one state if there is
	// ambiguity. For example, if all taxa have ?, then the site is potentially constant
	// for all states.
    buildConstantStatesVector();

	// The relative rate means and probabilities vectors need to be recalculated if the
	// number of rate categories subsequently changes
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves information about the original, uncompressed data to a file named `filename'. All of the information displayed
|	is saved by the TreeLikelihood::compressDataMatrix function.
*/
void TreeLikelihood::debugUncompressedDataInfo(
  std::string filename)	/**< is the name of the file that will hold the info about the uncompressed data matrix */
	{
	//unsigned sz = (unsigned)pattern_vect[0].size();
	std::ofstream outf(filename.c_str());
	outf << "index\tsubset\tpattern\n";

	unsigned total_patterns = (unsigned)pattern_vect.size();

    std::vector< unsigned > subset_indices;
    subset_indices.assign(_num_char, UINT_MAX);

    std::vector< std::string > pattern_strings;
    pattern_strings.resize(_num_char);

	//unsigned j = 0;	//index into constant states vector
	for (unsigned i = 0; i < total_patterns; ++i)
		{
		const uint_vect_t & sites = pattern_to_sites[i];//UINT_LIST
		const int8_vect_t & states = pattern_vect[i];
		unsigned subset_index = *(states.begin());

        // Create string showing pattern
        // start at 1 because first element is index of partition subset
        std::string pattern_string = "";
        for (unsigned k = 1; k < states.size(); ++k)
            {
            int8_t el = states[k];
            pattern_string << getStateStr(subset_index, el);
            }

        // Store string in pattern_strings vector for each site listed in the sites vector
        for (uint_vect_t::const_iterator j = sites.begin(); j != sites.end(); ++j)
            {
            unsigned which_site = *j;
            pattern_strings[which_site] = pattern_string;
            subset_indices[which_site] = subset_index;
            }

		//unsigned num_cs = constant_states[j];
		//outf << boost::str(boost::format("%12d\t%12d\t%12d\t") % i % sub % num_cs);

		//std::copy(states.begin() + 1, states.end(), std::ostream_iterator<int>(outf," "));
		//outf << "\t";

		//double cnt = pattern_counts[i];
		//outf << boost::str(boost::format("%12g\t") % cnt);

		//std::copy(sites.begin(), sites.end(), std::ostream_iterator<unsigned>(outf," "));
		//outf << "\t";

		//outf << std::endl;

		//j += num_cs + 1;
		}

	for (unsigned i = 0; i < _num_char; ++i)
		{
        outf << i << '\t' << subset_indices[i] << '\t' << pattern_strings[i] << '\n';
        }
	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves information about the compressed data to a file named `filename'. All of the information displayed is
|	saved by the TreeLikelihood::compressDataMatrix function.
*/
void TreeLikelihood::debugCompressedDataInfo(
  std::string filename)	/**< is the name of the file that will hold the info about the compressed data matrix */
	{
	// This function should be combined with listPatterns, which has similar functionality
	unsigned sz = (unsigned)pattern_vect[0].size();
	std::string format_str = boost::str(boost::format("%%12s\t%%12s\t%%12s\t%%%ds\t%%12s\t%%s") % (2*(sz - 1)));
	std::ofstream outf(filename.c_str());
	outf << "pcs = potentially constant states\n" << std::endl;
	outf << boost::str(boost::format(format_str) % "index" % "subset" % "pcs" % "pattern " % "count" % "sites") << std::endl;

	unsigned total_patterns = (unsigned)pattern_vect.size();
	unsigned j = 0;	//index into constant states vector
	for (unsigned i = 0; i < total_patterns; ++i)
		{
		const uint_vect_t & sites = pattern_to_sites[i];//UINT_LIST
		const int8_vect_t & states = pattern_vect[i];

		unsigned sub = *(states.begin());
		unsigned num_cs = constant_states[j];
		outf << boost::str(boost::format("%12d\t%12d\t%12d\t") % i % sub % num_cs);

		std::copy(states.begin() + 1, states.end(), std::ostream_iterator<int>(outf," "));
		outf << "\t";

		double cnt = pattern_counts[i];
		outf << boost::str(boost::format("%12g\t") % cnt);

		std::copy(sites.begin(), sites.end(), std::ostream_iterator<unsigned>(outf," "));
		outf << "\t";

		outf << std::endl;

		j += num_cs + 1;
		}
	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies data from `mat' to `pattern_vect' and `pattern_counts'. The `pattern_vect' vector holds the patterns while
|	`pattern_counts' holds the count of the number of sites having each pattern. Additionally, the vectors
|	`pattern_to_sites' and `charIndexToPatternIndex' are built: `pattern_to_sites' allows you to get a list of sites
|	given a specific pattern, and `charIndexToPatternIndex' lets you find the index of a pattern in `pattern_vect' and
|	`pattern_counts' given an original site index.
*/
unsigned TreeLikelihood::compressDataMatrix(
  const NxsCXXDiscreteMatrix    & mat,              /**< is the data source */
  const uint_vect_t             & partition_info)	/**< is a vector of indices storing the partition subset used by each site */
	{
	unsigned 						ntax 		= mat.getNTax();
	unsigned 						nchar 		= mat.getNChar();
	unsigned 						nsubsets 	= partition_model->getNumSubsets();

	//typedef std::vector< pattern_map_t > pattern_map_vect_t;
	pattern_map_vect_t 				pattern_map(nsubsets);
	pattern_to_sites_map_t 			pattern_to_sites_map;

	// If user has not decided to partition the data, then the partition_info vector should be empty
	bool default_partition = (partition_info.size() == 0 ? true : false);

	// If user did decide to partition the data, then the partition_info vector should have a length
	// equal to the number of sites
	if (!default_partition && partition_info.size() != nchar)
		throw XLikelihood(str(boost::format("Partition scheme accounts for %d characters; however, there are %d characters in the data matrix") % partition_info.size() % nchar));

	// Initialize the charIndexToPatternIndex vector, which will allow us to later locate the
	// pattern associated with any given site in the original data matrix
	charIndexToPatternIndex.assign(nchar, UINT_MAX);

    // Create actingWeights vector and copy the integer weights from mat into it
    // If there are no integer weights in mat, copy the floating point weights instead
    // if floating point weights have been defined
	const int_vect_t & iwts = mat.getIntWeightsConst();
	double_vect_t actingWeights(nchar, 1.0);
	if (!iwts.empty())
		{
		PHYCAS_ASSERT(iwts.size() >= nchar);
		for (unsigned j = 0; j < nchar; ++j)
			actingWeights[j] = (double)iwts.at(j);
		}
	else
		{
		const double_vect_t & dwts = mat.getDblWeightsConst();
		if (!dwts.empty())
			{
			actingWeights = dwts;
			PHYCAS_ASSERT(actingWeights.size() >= nchar);	//@POL why > nchar? shouldn't it just be == nchar?
			}
		}

    // Set corresponding actingWeights elements to zero if any characters have been excluded in mat
	const uint_set_t & excl = mat.getExcludedCharIndices();
	for (uint_set_t::const_iterator eIt = excl.begin(); eIt != excl.end(); ++eIt)
		{
		PHYCAS_ASSERT(*eIt < nchar);
		actingWeights[*eIt] = 0.0;
		}
	const double * wts = &(actingWeights[0]);	// PELIGROSO

	int8_vect_t pattern;
	for (unsigned j = 0; j < nchar;)
		{
		unsigned subset = default_partition ? 0 : partition_info[j];
		bool codon_model = partition_model->subset_model[subset]->isCodonModel();

		if (codon_model)
			{
			if (!default_partition)
				{
				// 2nd and 3rd positions should be assigned to same subset if codon model is used
				PHYCAS_ASSERT(partition_info.size() > j + 2);
				PHYCAS_ASSERT(partition_info[j + 1] == subset);
				PHYCAS_ASSERT(partition_info[j + 2] == subset);
				}

			// Suppose nchar=5 and j=3: not enough sites to make a second codon. In this example,
			// j+2 = 5, which equals nchar, so break out of loop over sites
			if (j+2 >= nchar)
				break;

			// here we (arbitrarily) use the max weight of any char in the codon
			pattern_count_t charWt = (wts ? std::max(wts[j], std::max(wts[j+1], wts[j+2])) : 1.0);
			if (charWt > 0.0)
				{
				//std::cerr << std::endl;

				// Build up a vector representing the pattern of state codes at this site
				pattern.clear();

				// The index of the partition subset fills the first slot in the pattern
				pattern.push_back((int8_t)subset);

				for (unsigned i = 0; i < ntax; ++i)
					{
					const int8_t *	row			= mat.getRow(i);
					const int8_t	code1		= row[j];
					const int8_t	code2		= row[j + 1];
					const int8_t	code3		= row[j + 2];
					bool			code1_ok	= (code1 >= 0 && code1 < 4);
					bool			code2_ok	= (code2 >= 0 && code2 < 4);
					bool			code3_ok	= (code3 >= 0 && code3 < 4);
					if (code1_ok && code2_ok && code3_ok)
						{
						//std::cerr << boost::str(boost::format("(%d,%d,%d)") % (unsigned)code1 % (unsigned)code2 % (unsigned)code3);
						const int8_t code = codon_state_codes[16*code1 + 4*code2 + code3]; // CGT = 27 = 16*1 + 4*2 + 3
						if (code > 60)
							throw XLikelihood(str(boost::format("Stop codon encountered for taxon %d at sites %d-%d") % (i+1) % (j+1) % (j+4)));
						pattern.push_back(code);
						}
					else
						{
						//std::cerr << "(???)";
						// if any site within codon is ambiguous, entire codon is treated as completely ambiguous
						pattern.push_back((int8_t)61);
						}
					}

				//std::cerr << boost::str(boost::format("\n%6d -> ") % j);

				storePattern(pattern_map[subset], pattern_to_sites_map, pattern, j, charWt, codon_model);

				//for (std::vector<int8_t>::const_iterator it = pattern.begin(); it != pattern.end(); ++it)
				//	{
				//	std::cerr << boost::str(boost::format("%s|") % codon_labels[(unsigned)(*it)]);
				//	}
				//std::cerr << std::endl;
				}
			j += 3;
			}
		else	// model for site j is not a codon model
			{
			unsigned nStates = getNumStates(subset);

			pattern_count_t charWt = (wts ? wts[j] : 1.0);
			if (charWt > 0.0)
				{
				// Build up a vector representing the pattern of state codes at this site
				pattern.clear();

				// The index of the partition subset fills the first slot in the pattern
				pattern.push_back((int8_t)subset);

				unsigned num_all_missing = 0;
				for (unsigned i = 0; i < ntax; ++i)
					{
					const int8_t * row	= mat.getRow(i);
					const int8_t   code = row[j];

					if (default_partition)
						{
						pattern.push_back(code);
						if ((unsigned)code == nStates)
							++num_all_missing;
						}
					else
						{
						//@POL any ambiguity is treated as ? for partitioned data (for now)
						//@POL if this is relaxed, need to revisit TreeLikelihood::copyDataFromDiscreteMatrix,
						//@POL which currently constructs a state_list without partial ambiguities
						if (code >= 0 && code < (state_code_t)nStates)	//POLBM
							pattern.push_back(code);
						else
							{
							pattern.push_back((state_code_t)nStates);
							++num_all_missing;
							}
						}
					}

				// Do not include the pattern if it contains only completely missing data
				// for all taxa
				if (num_all_missing == ntax)
					{
					all_missing.push_back(j++);
					continue;
					}

				storePattern(pattern_map[subset], pattern_to_sites_map, pattern, j, charWt, codon_model);
				}
			++j;
			}	// not a codon model
		}	// loop over sites

// ***** below here now duplicated by patternMapToVect *****

	// Build subset_offset, pattern_counts, pattern_vect, pattern_to_sites and charIndexToPatternIndex before
	// leaving this function (whereupon pattern_map and pattern_to_sites_map will be destroyed)
	unsigned npatterns = 0;
	std::vector<unsigned> npatterns_vect(nsubsets, 0);
	for (unsigned i = 0; i < nsubsets; ++i)
		{
		unsigned np = (unsigned)pattern_map[i].size();
		npatterns_vect[i] = np;
		npatterns += np;
		}

	subset_offset.clear();
	subset_offset.reserve(nsubsets + 1);

	pattern_counts.clear();
	pattern_counts.reserve(npatterns);

	pattern_vect.clear();
	pattern_vect.reserve(npatterns);

	pattern_to_sites.clear();
	pattern_to_sites.reserve(npatterns);

	unsigned pattern_index = 0;
	unsigned n_inc_chars = 0;
	std::vector<unsigned> nsites_vect(nsubsets, 0);
	for (unsigned i = 0; i < nsubsets; ++i)
		{
		unsigned num_sites_this_subset = 0;
		subset_offset.push_back(pattern_index);
		for (pattern_map_t::iterator mapit = pattern_map[i].begin(); mapit != pattern_map[i].end(); ++mapit)
			{
            // mapit->first holds the pattern in the form of a vector of int8_t values
            pattern_vect.push_back(mapit->first);

            // mapit->second holds the pattern count
            pattern_counts.push_back(mapit->second);
            num_sites_this_subset += (unsigned)mapit->second;

            // get list of sites that had this pattern
            const uint_vect_t & sites = pattern_to_sites_map[mapit->first];//UINT_LIST

            // add this sites list to pattern_to_sites vector
            pattern_to_sites.push_back(sites);

            // For each site index in the sites list, add an element to the map charIndexToPatternIndex
            // Now, charIndexToPatternIndex[i] points to the index in pattern_vect for the pattern found at site i
            for (uint_vect_t::const_iterator sitesIt = sites.begin(); sitesIt != sites.end(); ++sitesIt)//UINT_LIST
                {
                charIndexToPatternIndex[*sitesIt] = pattern_index;
                ++n_inc_chars;
                }

            ++pattern_index;
			}	// loop over patterns in subset i
		nsites_vect[i] = num_sites_this_subset;
		}	// loop over subsets
	partition_model->setNumSitesVect(nsites_vect);
    partition_model->setNumPatternsVect(npatterns_vect);

	PHYCAS_ASSERT(partition_model->getTotalNumPatterns() == pattern_index);
	subset_offset.push_back(pattern_index);
// ***** above here now duplicated by patternMapToVect *****

	// There should no longer be any elements in charIndexToPatternIndex that have the value UINT_MAX
	// If there are, the elements that still have the value UINT_MAX should correspond with indices stored in the all_missing vector
	// or the excl set (excluded characters)
	PHYCAS_ASSERT(!excl.empty() || !all_missing.empty() || (std::find(charIndexToPatternIndex.begin(), charIndexToPatternIndex.end(), UINT_MAX) == charIndexToPatternIndex.end()));

	// pattern_map and pattern_to_sites_map are just temporary containers. The information originally in
	// pattern_map is now in pattern_vect and pattern_counts, and the information originally in
	// pattern_to_sites_map is now in charIndexToPatternIndex and pattern_to_sites.
	pattern_map.clear();
	pattern_to_sites_map.clear();

	return npatterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_vect' and `pattern_counts' using data stored in supplied `pattern_map'.
*/
void TreeLikelihood::patternMapToVect(
  const pattern_map_vect_t & pattern_map_vect,          /**< is a vector of maps (one element per subset) whose keys are patterns and whose values are counts */
  pattern_to_sites_map_t & pattern_to_sites_map)  /**< is a map whos keys are patterns and whose values are lists of site indices */
    {
	unsigned nsubsets = partition_model->getNumSubsets();

	// Build subset_offset, pattern_counts, pattern_vect, pattern_to_sites and charIndexToPatternIndex before
	// leaving this function (whereupon pattern_map and pattern_to_sites_map will be destroyed)
	unsigned npatterns = 0;
	std::vector<unsigned> npatterns_vect(nsubsets, 0);
	for (unsigned i = 0; i < nsubsets; ++i)
		{
		unsigned np = (unsigned)pattern_map_vect[i].size();
		npatterns_vect[i] = np;
		npatterns += np;
		}

	subset_offset.clear();
	subset_offset.reserve(nsubsets + 1);

	pattern_counts.clear();
	pattern_counts.reserve(npatterns);

	pattern_vect.clear();
	pattern_vect.reserve(npatterns);

	pattern_to_sites.clear();
	pattern_to_sites.reserve(npatterns);

	unsigned pattern_index = 0;
	unsigned n_inc_chars = 0;
	std::vector<unsigned> nsites_vect(nsubsets, 0);
	for (unsigned i = 0; i < nsubsets; ++i)
		{
		unsigned num_sites_this_subset = 0;
		subset_offset.push_back(pattern_index);
		for (pattern_map_t::const_iterator mapit = pattern_map_vect[i].begin(); mapit != pattern_map_vect[i].end(); ++mapit)
			{
            // mapit->first holds the pattern in the form of a vector of int8_t values
            pattern_vect.push_back(mapit->first);

            // mapit->second holds the pattern count
            pattern_counts.push_back(mapit->second);
            num_sites_this_subset += (unsigned)mapit->second;

            // get list of sites that had this pattern
            const uint_vect_t & sites = pattern_to_sites_map[mapit->first];//UINT_LIST

            // add this sites list to pattern_to_sites vector
            pattern_to_sites.push_back(sites);

            // For each site index in the sites list, add an element to the map charIndexToPatternIndex
            // Now, charIndexToPatternIndex[i] points to the index in pattern_vect for the pattern found at site i
            for (uint_vect_t::const_iterator sitesIt = sites.begin(); sitesIt != sites.end(); ++sitesIt)//UINT_LIST
                {
                charIndexToPatternIndex[*sitesIt] = pattern_index;
                ++n_inc_chars;
                }

            ++pattern_index;
            }	// loop over patterns in subset i
		nsites_vect[i] = num_sites_this_subset;
		}	// loop over subsets
	partition_model->setNumSitesVect(nsites_vect);
    partition_model->setNumPatternsVect(npatterns_vect);

	PHYCAS_ASSERT(partition_model->getTotalNumPatterns() == pattern_index);
	subset_offset.push_back(pattern_index);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_vect' and `pattern_counts' using data stored in `sim_data'.
*/
void TreeLikelihood::copyDataFromSimData(
  SimDataShPtr sim_data)	/**< is the data source */
	{
    // formerly DISABLED_UNTIL_SIMULATION_WORKING_WITH_PARTITIONING
	// Copy simulated data to pattern_map
	pattern_map_vect_t pattern_map_vect;
    pattern_map_vect.push_back(sim_data->getSimPatternMap()); //POLSIM: 0 is first and only subset, need to generalize
	pattern_to_sites_map_t & pattern_to_sites_map = sim_data->getPatternToSitesMap();
    patternMapToVect(pattern_map_vect, pattern_to_sites_map);

	// Build up counts vector
	//pattern_counts.clear();
	//for (const pattern_map_t::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
	//	{
	//	pattern_counts.push_back(it->second);
	//	}

	//_num_tax = sim_data->getPatternLength();
	//num_patterns = (unsigned)pattern_map.size();

	//model->buildStateList(state_list, state_list_pos);

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `pattern_vect' to the patterns already in `other'. Assumes that `pattern_length' for
|	this SimData object is identical to the `pattern_length' of `other'.
*/
void TreeLikelihood::addDataTo(SimData & other)
	{
	if (pattern_vect.empty())
		return;

	// If other is empty, then it most likely needs to be initialized
	if (other.getTotalCount() == 0)
		{
		PHYCAS_ASSERT(_num_taxa > 0);
		other.resetPatternLength(_num_taxa);
		}
	PHYCAS_ASSERT(_num_taxa == other.getPatternLength());

    // Make sure SimData::sim_pattern_vect is long enough to accommodate all sites
    other.resizePatternVect(_num_char);

    count_vect_t::iterator cit = pattern_counts.begin();
    pattern_vect_t::iterator pit = pattern_vect.begin();
	for (; pit != pattern_vect.end() && cit != pattern_counts.end(); ++pit, ++cit)
		{
		pattern_count_t count = *cit;

		pattern_t & other_pattern = other.getCurrPattern();

        // skip first element of pattern because that is the subset index
		std::copy(pit->begin() + 1, pit->end(), other_pattern.begin());

		//other.insertPattern(*sit, count); //POLOLD
		other.insertPatternOnly(count); //POLOLD
		}
    PHYCAS_ASSERT(pit == pattern_vect.end());
    PHYCAS_ASSERT(cit == pattern_counts.end());
    //PHYCAS_ASSERT(sit == pattern_to_sites.end()); //POLOLD
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateInternalData() to add an InternalData structure to `nd' containing the conditional likelihood arrays
|	needed for likelihood calculations.
*/
void TreeLikelihood::prepareInternalNodeForLikelihood(
  TreeNode * nd)	/**< is the node to decorate */
	{
	if (nd)
		{
    	InternalData * ndID = nd->GetInternalData();
        if (ndID == NULL)
            {
    		TreeNode::InternalDataDeleter cl_deleter = &deallocateInternalData;
	    	InternalData * cl = allocateInternalData();
		    nd->SetInternalData(cl, cl_deleter);
            }
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Decorates an interior node with conditional likelihood and transition probability structures and add to the tree's
|	StoreInternalNode vector.
*/
void TreeLikelihood::addDecoratedInternalNode(
  TreeShPtr t,		/**< is the tree to decorate */
  unsigned num)		/**< is the number to assign to this node */
	{
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
	InternalData * cl = allocateInternalData();
	TreeNode * nd = t->AllocNewNode();
	nd->SetNodeNum(num);
	nd->SetInternalData(cl, cl_deleter);
	t->StoreInternalNode(nd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Decorates a tip node with data and add to the tree's tipStorage vector.
*/
void TreeLikelihood::addOrphanTip(
  TreeShPtr t,		/**< is the tree to decorate */
  unsigned row,		/**< is the row in the data matrix to associate with the added tip */
  std::string name) /**< is the taxon name */
	{
	TreeNode::TipDataDeleter td_deleter = &deallocateTipData;
	TipData * td = allocateTipData(row);
	TreeNode * nd = t->AllocNewNode();
	nd->SetTipData(td, td_deleter);
	nd->SetNodeNum(row);
	nd->SetNodeName(name);
	t->StoreLeafNode(nd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateTipData() for all tip nodes and allocateInternalData() for all internal nodes in the supplied Tree.
|	Assumes each tip node number in the tree equals the appropriate row in the data matrix.
*/
void TreeLikelihood::prepareForLikelihood( //POL_BOOKMARK TreeLikelihood::prepareForLikelihood
  TreeShPtr t)		/**< is the tree to decorate */
	{
	// If no_data is true, it means that calcLnL will always return 0.0 immediately and
	// will thus never need the TipData or InternalData data structures
	if (no_data)
		return;

	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	// Put all existing conditional likelihood arrays already back into storage
	//storeAllCLAs(t);	//@POL why are these two lines commented out? Seems like a good idea to clear the cla stack at this point.
	//cla_pool->clearStack();	//@POL this here only because prepareForLikelihood called in NCatMove::proposeNewState when ncat is increased

	preorder_iterator nd = t->begin();

	if (t->IsRooted())		//skip artificial root node for rooted trees
		++nd;

	for ( ; nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			unsigned row = nd->GetNodeNumber();
			TipData * td = allocateTipData(row);
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	useAsLikelihoodRoot(NULL);

	//@POL	should decorate any nodes in tipStorage and nodeStorage now, but should wait until those
	// are changed to vectors
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns reference to the vector data member `all_missing', which is a list of indices of sites for which all taxa
|   show completely missing data. The first site has index 0 in the `all_missing' vector.
*/
const std::vector<unsigned> & TreeLikelihood::getListOfAllMissingSites() const
    {
    return all_missing;
    }

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|   Returns a list of all sites with patterns comprising only two primary states and no missing data or ambiguities for
|   any taxon.
*/
std::vector<unsigned> TreeLikelihood::findDataBipartitions() const
    {
    // needs to be written
    }
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Builds up the vector data member `constant_states' based on the patterns in `pattern_vect'. The `constant_states'
|   vector has a similar structure to the `state_list' vector. Here is an example `constant_states' vector for 4 sites:
|>
|   +---+---+---+---+---+---+---+---+
|   | 0 | 1 | 0 | 1 | 3 | 2 | 0 | 2 |
|   +---+---+---+---+---+---+---+---+
|     ^   ^       ^       ^
|     |   |       |       |
|     |   |       |       4th site can potentially be constant for either A (state 0) or G (state 2)
|     |   |       3rd site can be potentially constant for T (state 3)
|     |   2nd site can be potentially constant for 1 state: that state (A, state 0) follows in the next cell
|     1st site is definitely variable (hence the 0 meaning no states follow)
|>
|   The function returns the number of potentially constant patterns found (note the return value is NOT the number
|   of potentially constant sites). The modifier "potentially" is needed because of ambiguities: i.e. if every taxon
|	had missing data for a site, then one does not know if the site is constant or variable, but it is a "potentially"
|	constant site.
*/
unsigned TreeLikelihood::buildConstantStatesVector()
	{
    PHYCAS_ASSERT(!pattern_vect.empty());
    constant_states.clear();

    unsigned 	num_potentially_constant = 0;
    int_set_t 	common_states;		// holds running intersection across taxa for this pattern
    int_set_t 	curr_states;		// holds states for taxon under consideration
    int_vect_t 	xset;				// temporarily holds intersection of common_states and curr_states

    for (pattern_vect_t::const_iterator pat = pattern_vect.begin(); pat != pattern_vect.end(); ++pat)
		{
		// get subset-specific info
		//@POL if slow, consider using subset_offset here to avoid setting subset and ns for each pattern
		const unsigned subset = (unsigned)(*pat)[0];
        bool site_potentially_constant = true;
        for (unsigned taxon = 0; taxon < _num_taxa; ++taxon)
            {
            if (!site_potentially_constant)
                break;

            // Reminder of how data members state_list and state_list_pos are laid out (for nucleotide data anyway):
            //                                         ?         N      {CGT}  {ACG}    R,(AG)    Y
            // states          | A | C | G | T |  - A C G T | A C G T | C G T | A G T |  A  G  | C T
            // state_list      1 0 1 1 1 2 1 3 5 -1 0 1 2 3 4 0 1 2 3 3 1 2 3 3 0 2 3 2  0  2  4 1 3
            // state_list_pos  0   2   4   6   8           14        19      23      27       30
			unsigned index_of_taxon = taxon + 1;	// add 1 because first element of pattern is the subset index
            int code = (int)(*pat)[index_of_taxon];
            PHYCAS_ASSERT(code >= 0);   // Mark, why don't ? and - states trigger this assert? Does this have to do with the fact that we have abandoned the Cipres version of NCL? We translate - to ? in the NxsCXXDiscreteMatrix constructor
            PHYCAS_ASSERT(code < (int)state_list[subset].size());
            PHYCAS_ASSERT(code < (int)state_list_pos[subset].size());
            unsigned pos = (unsigned)state_list_pos[subset][code];
            unsigned nstates = (unsigned)state_list[subset][pos];
			++pos;	// skip the number of states so that pos indicates index of first (of perhaps several) possible state(s)

            // Insert all states for the current taxon into the curr_states set
            curr_states.clear();
            for (unsigned x = pos; x < pos + nstates; ++x)
                {
                int c = (int)state_list[subset][x];
                PHYCAS_ASSERT(c >= -1);
                if (site_potentially_constant)
                    curr_states.insert(c);
                }

            if (taxon == 0)
                {
                // For first taxon, let common_states equal curr_states
                common_states.clear();
                common_states.insert(curr_states.begin(), curr_states.end());
                }
            else
                {
                // For subsequent taxa, common_states should contain only those states possessed
                // by all taxa (i.e. it will be the empty set if the site is variable)
                xset.clear();
                std::set_intersection(common_states.begin(), common_states.end(), curr_states.begin(), curr_states.end(), back_inserter(xset));
                if (xset.empty())
                    site_potentially_constant = false;
                else
                    {
                    common_states.clear();
                    common_states.insert(xset.begin(), xset.end());
                    }
                }
            }

        if (site_potentially_constant)
            {
            ++num_potentially_constant;
            constant_states.push_back((unsigned)common_states.size());
            for (std::set<int>::const_iterator it = common_states.begin(); it != common_states.end(); ++it)
                {
                constant_states.push_back(*it);
                }
            }
        else
            {
            constant_states.push_back(0);
            }
        }
    return num_potentially_constant;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Increments the count corresponding to the supplied `pattern' in `pattern_vect' (adding a new map element if `pattern'
|	has not yet been seen) and adds the `pattern_index' to a list of pattern indices corresponding to the supplied
|	`pattern' in `patternToIndex' (adding a newly-created list containing just `pattern_index' if `pattern' has not yet
|	been seen.
*/
void TreeLikelihood::storePattern(
  pattern_map_t & pattern_map,					/**< is the map into which to store pattern */
  pattern_to_sites_map_t & pattern_to_site_map,	/**< is the map into which to store the original site index of this pattern */
  const std::vector<int8_t> & pattern,			/**< is the pattern to store */
  const unsigned int site_index,				/**< is the index of the site in the original data matrix */
  const pattern_count_t weight, 				/**< is the weight of this pattern (normally 1.0) */
  bool codon_model)								/**< if true, site_index and following two sites will be stored in the vector of sites having this pattern; if false, only site_index will be stored */
	{
	// Add the pattern to pattern_map if it has not yet been seen, otherwise increment
	// the count of this pattern if it is already in pattern_map (see item 24, p. 110, in Meyers' Efficient STL)
	pattern_map_t::iterator lowb = pattern_map.lower_bound(pattern);
	if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
		{
		// pattern is already in pattern_map, increment count
		lowb->second += weight;
		//std::cerr << boost::str(boost::format("%.1f -> |") % lowb->second);
		}
	else
		{
		// pattern has not yet been stored in pattern_map
		pattern_map.insert(lowb, pattern_map_t::value_type(pattern, weight));
		//std::cerr << boost::str(boost::format("%.1f -> |") % weight);
		}

	// Add the pattern to pattern_to_site_map if it has not yet been seen, otherwise increment
	// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
	pattern_to_sites_map_t::iterator pToILowB = pattern_to_site_map.lower_bound(pattern);
	if (pToILowB != pattern_to_site_map.end() && !(pattern_to_site_map.key_comp()(pattern, pToILowB->first)))
		{
		// a list of site indices has already been created for this pattern in pattern_to_site_map
		// so all we need to do is append the index site_index to the existing list
		pToILowB->second.push_back(site_index);
		if (codon_model)
			{
			pToILowB->second.push_back(site_index + 1);
			pToILowB->second.push_back(site_index + 2);
			}
		}
	else
		{
		// this pattern has not yet been seen, so need to create a list of site indices whose only
		// element (so far) is site_index and insert this list into the pattern_to_site_map map
		uint_vect_t ilist(1, site_index);	// create a list containing 1 element whose value is site_index //UINT_LIST
		if (codon_model)
			{
			ilist.push_back(site_index + 1);
			ilist.push_back(site_index + 2);
			}
		pattern_to_site_map.insert(pToILowB, pattern_to_sites_map_t::value_type(pattern, ilist));
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a string representation of the supplied `state'. For example, if `state' equals 1 (where model is a
|	standard DNA model), the string returned would be "C". If, however, `state' was 7 (again, standard DNA model), then
|	there is some ambiguity present and the string returned might look something like "{AC}" (the string returned
|	depends of course on the actual meaning of the global state code 7).
*/
std::string TreeLikelihood::getStateStr(
  unsigned i,				/**< is the partition subset */
  state_code_t state) const /**< is the global state code to be converted to a std::string */
	{
	std::string s;
	const state_code_t ns = partition_model->subset_num_states[i];
	const state_code_t nsPlusOne = ns + 1;

	if (state < nsPlusOne)
		{
		// either no ambiguity or complete ambiguity
		s << partition_model->subset_model[i]->lookupStateRepr((int)state);
		}
	else
		{
		// `state' represents partial ambiguity

		// First, find location of the definition of `state' in the global state list
		unsigned pos = state_list_pos[i][(unsigned)state];
		pattern_t::const_iterator it = state_list[i].begin() + pos;

		// Now get the number of basic states composing `state'
		unsigned n = *it++;

		// Walk down global state list converting states into strings
		s << "(";
		for (unsigned k = 0; k < n; ++k)
			{
			s << partition_model->subset_model[i]->lookupStateRepr((int)*it++);
			}
		s << ")";
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets private data member `nevals' to 0.
*/
void TreeLikelihood::resetNumLikelihoodEvals()
	{
	nevals = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor providing access to the current value of the private data member `nevals'.
*/
unsigned TreeLikelihood::getNumLikelihoodEvals() const
	{
	return nevals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor providing access to the current value of the private data member `nevals'.
*/
void TreeLikelihood::incrementNumLikelihoodEvals()
	{
	++nevals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assuming TreeLikelihood::compressDataMatrix has been called, so that `pattern_vect' is up-to-date, returns a string
|	listing all observed patterns and their frequencies.
*/
std::string TreeLikelihood::listPatterns(
  bool show_coded_states)	/**< if true, output global state codes used internally; otherwise, do the translation back to the original state representations */
	{
	string_vect_t		state_repr;
	std::string 		s;
	unsigned 			nsubsets = partition_model->getNumSubsets();

	uint_vect_t						pattern_tally(nsubsets, 0);
	std::vector< pattern_count_t >	sites_tally(nsubsets, 0);

	for (unsigned i = 0; i < nsubsets; ++i)
		{
		for (unsigned j = subset_offset[i]; j < subset_offset[i + 1]; ++j)
			{
			const pattern_t & 	p = pattern_vect[j];
			pattern_count_t 	c = pattern_counts[j];

			// Increment tally of patterns and sites for this subset
			pattern_tally[i] += 1;
			sites_tally[i] += c;

			s << str(boost::format("%6d %6d %6.1f ") % j % i % c);

			if (show_coded_states)
				{
				// start at 1 because first element is index of partition subset
				for (unsigned k = 1; k < p.size(); ++k)
					{
					int8_t el = p[k];
					s << str(boost::format("%d ") % (int)(el));
					}
				}
			else
				{
				// start at 1 because first element is index of partition subset
				for (unsigned k = 1; k < p.size(); ++k)
					{
					int8_t el = p[k];
					s << str(boost::format("%s") % getStateStr(i,el));
					}
				}
			s << '\n';
			}
		}

	s << "\nTotals for each subset of the partition:\n";
	s << boost::str(boost::format("\n%12s %12s %12s") % "subset" % "patterns" % "sites");
	for (unsigned i = 0; i < nsubsets; ++i)
		{
		s << boost::str(boost::format("\n%12d %12d %12d") % i % pattern_tally[i] % sites_tally[i]);
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First erases `pattern_repr', then builds it into a vector of pattern representation strings. When the function
|	returns, each element of `pattern_repr' is a string containing global state codes separated by whitespace. This
|	function is intended to be used for debugging purposes, where it is sometimes helpful to be sure of exactly which
|	data pattern is being operated upon. Thus, no particular effort has been made to make this function efficient.
*/
void TreeLikelihood::buildPatternReprVector(string_vect_t & pattern_repr, TreeShPtr t)
	{
	unsigned nTips	= t->GetNTips();
	unsigned npatterns = partition_model->getTotalNumPatterns();

	pattern_repr.clear();
	pattern_repr.reserve(npatterns);

	// Recreate the data matrix in terms of tip-specific state codes
	state_code_t * * m = NewTwoDArray<state_code_t>(nTips, npatterns);

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			// get node number i
			unsigned i = nd->GetNodeNumber();
			PHYCAS_ASSERT(i >= 0 && i < nTips);

			// get tip-specific state code array
			TipData * tipData = nd->GetTipData();
			PHYCAS_ASSERT(tipData);

			// fill row i of data matrix
			unsigned subset = 0;
			for (unsigned j = 0; j < npatterns; ++j)
				{
				if (j >= subset_offset[subset + 1])
					subset++;
				const state_code_t * tipCodes = tipData->getConstStateCodes(subset);

				// get position vector that allows us to translate tip-specific state codes into global state codes
				const state_list_pos_t & local_statelist_pos = tipData->getConstStateListPos(subset);

				int8_t global_code = tipCodes[j];
				int8_t offset = global_code - ((int8_t)npatterns + 1);
				if (offset >= 0)
					{
					unsigned pos = local_statelist_pos[offset];
					for (unsigned m = 0; m < (unsigned)state_list_pos[subset].size(); ++m)
						{
						if (state_list_pos[subset][m] == pos)
							global_code = (state_code_t)m;
						}
					}
				m[i][j] = global_code;
				}
			}
		}

	for (unsigned j = 0; j < npatterns; ++j)
		{
		std::string s;
		for (unsigned i = 0; i < nTips; ++i)
			{
			s << (int)m[i][j] << " ";
			}
		pattern_repr.push_back(s);
		}

	DeleteTwoDArray<state_code_t>(m);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The node `slider' is moved along the combined edges of `slider' and `other'. If `fraction' is positive, `slider'
|	gains edge length at the expense of `other'. If `fraction' is negative, `other' gains edge length at the expense of
|	`slider'. This function not only changes edge lengths, but also takes care of transferring any univents that lie
|	along the affected path. If `slider' is a sibling of other, then `slider' is somewhat of a misnomer since it does
|	not actually slide; instead, the parent of both `slider' and `other' does the sliding.
*/
void TreeLikelihood::slideNode(
  double fraction,
  TreeNode * slider,
  TreeNode * other)
	{
	bool other_is_child = (other->GetParent() == slider);
	PHYCAS_ASSERT(other_is_child || slider->GetParent() == other->GetParent());
	if (other_is_child)
		{
		if (fraction < 0.0)
			{
			// Fraction is negative, which means that other gains edge length at its base at the
			// expense of slider, which loses edge length at its end
			//						   |
			//						   |
			//	-----------------------*--------------------------------*
			//				 <---------^slider							^other
			//	<----xnew---><------------------ynew-------------------->
			//	<-----------x----------><--------------y---------------->
			//				 <--delta-->
			//
			//PHYCAS_ASSERT(which_case == 2 || which_case == 5 || which_case == 11);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = -x*fraction;
			double xnew	  = x - delta;
			double ynew	  = y + delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			}	// fraction negative
		else
			{
			// Fraction is positive, so slider gains edge length at its end at the expense of other,
			// which loses edge length at its base
			//						   |
			//						   |
			//	-----------------------*--------------------------------*
			//					 slider^-------------->					^other
			//	<------------------xnew---------------><------ynew------>
			//	<----------x----------><--------------y----------------->
			//						   <-----delta---->
			//
			//PHYCAS_ASSERT(which_case == 1 || which_case == 4 || which_case == 10);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = y*fraction;
			double xnew	  = x + delta;
			double ynew	  = y - delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			}	// fraction positive
		}	// other is child
	else
		{
		// other is the sibling, not the child, of slider
		if (fraction < 0.0)
			{
			// Fraction is negative, which means that other gains edge length at its base at the
			// expense of slider, which loses edge length at its base
			//
			//	   +-----------x-----------* slider
			//	---+
			//	   +----------------y-------------* other
			//
			//				<-----delta---->
			//
			//	   +--xnew--* slider
			//	---+
			//	   +----------------------ynew----------------------* other
			//
			//PHYCAS_ASSERT(which_case == 8);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = -x*fraction;
			double xnew	  = x - delta;
			double ynew	  = y + delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			}	// fraction negative
		else
			{
			// Fraction is positive, so slider gains edge length at its base at the expense of other,
			// which loses edge length at its base
			//
			//	   +-----------x-----------* slider
			//	---+
			//	   +----------------y-------------* other
			//
			//							   <-----delta----->
			//
			//	   +------------------xnew------------------* slider
			//	---+
			//	   +---------ynew------* other
			//
			//PHYCAS_ASSERT(which_case == 7);
			double x	  = slider->GetEdgeLen();
			double y	  = other->GetEdgeLen();
			double delta  = y*fraction;
			double xnew	  = x + delta;
			double ynew	  = y - delta;
			slider->SetEdgeLen(xnew);
			other->SetEdgeLen(ynew);
			}	// fraction positive
		}	// other is sibling
	}	// TreeLikelihood::slideNode

}	// namespace phycas

