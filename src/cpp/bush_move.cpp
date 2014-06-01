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

#include <fstream>
#include <boost/format.hpp>

//#include "phycas/force_include.h"
#include "probability_distribution.hpp"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "xlikelihood.hpp"
#include "mcmc_chain_manager.hpp"
#include "bush_move.hpp"
#include "basic_tree.hpp"
#include "tree_manip.hpp"

//#define KEEP_BRANCHES_LEGAL
//#define MAX_LEGAL_BRLEN 285

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets `num_taxa' and `num_nodes_in_fully_resolved_tree' to zero, `edgelen_mean' to 1.0, `is_move' to
|   true, creates the `topo_prior_calculator', and finally calls BushMove::reset to re-initialize all other variables.
*/
BushMove::BushMove()
  : MCMCUpdater()
	{
	is_move = true;
	edgelen_mean = 1.0;

	num_taxa = 0;
	num_nodes_in_fully_resolved_tree = 0;

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does not need to do anything.
*/
BushMove::~BushMove()
    {
    //std::cerr << "BushMove dying..." << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true because this class is capable of computing the tree topology prior.
*/
bool BushMove::computesTopologyPrior() const
    {
    return true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or
|	revert(), whichever is appropriate.
*/
bool BushMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed BushMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	if (num_taxa == 0)
		{
		throw XLikelihood("Must call finalize() before calling update() for BushMove");
		}

    // chain_mgr is a weak_ptr: calling lock() creates a shared_ptr which ensures that the ChainManager
    // object pointed to does not go away before we are done using it
	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

	double prev_ln_like = p->getLastLnLike();
	double prev_ln_prior = jpm->getLogJointPrior();

	proposeNewState();

	curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
	curr_ln_prior = jpm->getLogJointPrior();

	//likelihood->startTreeViewer(tree, str(boost::format("Bush move PROPOSED (%s)") % (add_edge_move_proposed ? "add edge" : "delete edge")));

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;

    if (is_standard_heating)
        {
        prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
	    curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
        PHYCAS_ASSERT(!use_ref_dist);   // Bush move not ready for generalized steppingstone yet
        }
    else
        {
        prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
	    curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
        }

	double ln_accept_ratio = curr_posterior - prev_posterior + ln_hastings + ln_jacobian;

    bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform()) <= ln_accept_ratio);
    if (save_debug_info)
        {
    	if (add_edge_move_proposed)
            {
            debug_info = str(boost::format("Bush: add (%s)") % (accepted ? "accepted" : "rejected"));
            }
        else
            {
            debug_info = str(boost::format("Bush: del (%s)") % (accepted ? "accepted" : "rejected"));
            }
        }

	if (accepted)
		{
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
		curr_ln_like = p->getLastLnLike();
		revert();
		return false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Proposes either an add-edge move (addition of an edge to a polytomy) or a delete-edge move (deletion of an existing
|	edge) with equal probability. If tree is currently fully-resolved, proposes only a delete-edge move. If tree is
|	currently the star tree, proposes only an add-edge move.
*/
void BushMove::proposeNewState()
	{
	unsigned nnodes = tree->GetNNodes();
	const bool fully_resolved_before = (num_nodes_in_fully_resolved_tree == nnodes);
	const bool star_tree_before = (nnodes == num_taxa + 1);
	add_edge_move_proposed = (!fully_resolved_before) && (star_tree_before || rng->Uniform() < 0.5);

	refreshPolytomies();
	num_polytomies = (unsigned)polytomies.size();

	// One internal edge leads to the tip node serving as the root, but that one doesn't count
	const unsigned num_internal_edges_before = tree->GetNInternals() - 1;

	TreeNode * nd = NULL;
	if (add_edge_move_proposed)
		{
		// Choose a polytomy at random to split
		//
		PHYCAS_ASSERT(num_polytomies > 0);
		unsigned i = rng->SampleUInt(num_polytomies);
		nd = polytomies[i];
		proposeAddEdgeMove(nd);

		// Compute the Jacobian by working backwards to obtain the uniform used to sample new_edgelen (done in
		// proposeAddEdgeMove function). Remember that U = 1 - e^{-theta*v} where v is the new edge length,
		// and theta = 1/edgelen_mean is the hazard parameter of the exponential distribution used for
		// sampling new edge lengths. log(1 - U) = -theta*v = -v/edgelen_mean. The Jacobian for this move
		// can be obtained using the Green (2003) approach. The new edge length v that exists after the
		// add-edge move corresponds with the uniform random deviate U before the add-edge move. The
		// Jacobian is the matrix of the derivatives of all post-move quantities with respect to all
		// pre-move quantities. In this case, there is only one post-move quantity (v) and only one pre-move
		// quantity (U), so the Jacobian is a 1x1 matrix consisting only of the derivative of v with respect
		// to U. The quantity v can be written as -edgelen_mean*log(1 - U), so the partial derivative is
        // simply edgelen_mean/(1 - U). The log of the jacobian is thus:
		//
		//      ln_jacobian = log(edgelen_mean) - log(1 - U) = log(edgelen_mean) + v/edgelen_mean
		//
		ln_jacobian = log(edgelen_mean) + new_edgelen/edgelen_mean;

		// The Hastings ratio is the probability of the reverse move divided by probability of the forward move.
		// The add-edge move involves the following 4 steps:
		//   1) choosing a polytomy to break (the probability of this is the inverse of the number of polytomies)
		//   2) choosing the number of spokes to move over to the new node
		//   3) selecting exactly which spokes to move
		//   4) choosing U ~ Uniform(0,1) to determine the new edge length (probability density is 1).
		// The probability of steps 2 and 3 are easier to compute jointly because of some convenient cancellation:
		// For step 2, assuming x spokes out of n total are moved, the probability is (see explanation in
		// documentation for BushMove::computePolytomyDistribution):
		//
		//              {n choose x} + {n choose (n-x)}
		// Pr(step 2) = -------------------------------
		//                       2^n - 2*(n+1)
		//
		// If x is exactly half of n, the numerator will have only the one term (n choose x), whereas the
		// denominator will remain the same. For example, let the tree be the 4-taxon star tree (A,B,C,D),
		// so n = 4, and let 2 nodes (A and C) be moved, then
		//
		//              {4 choose 2}       6
		// Pr(step 2) = ------------- = ------- = 1
		//              2^4 - 2*(4+1)   16 - 10
		//
		// For step 3, the probability is simply 1 over the number of distinct trees in which an n-spoke
		// polytomy has been divided into a node with x+1 spokes and another with n - x + 1 spokes:
		//
		//                              2
		// Pr(step 3) = -------------------------------
		//               {n choose x} + {n choose n-x}
		//
		// Again, if x is exactly half of n, the denominator will have only the one term {n choose x}. For
		// example, continuing the example above with n = 4 and x = 2, and the tree after the add-edge
		// move is (B, D, (A, C))
		//
		//                   2        2     1
		// Pr(step 3) = ---------- = --- = ---
		//              4 choose 2    6     3
		//
		// The 6 in the denominator corresponds with these possibilities, and the 2 in the numerator are
		// indicated by asterisks:
		//
		//  A     C   A     B   A     B   B     A   B     A   C     A
		//   \___/     \___/     \___/     \___/     \___/     \___/
		//   /   \     / * \     /   \     /   \     / * \     /   \  <not a backslash>
		//  B     D   C     D   C     D   C     D   D     C   D     B
		//
		//
		// Fortunately, the hard part of both calculations cancels, and the probability of steps 2 and 3 is:
		//
		//                                2                1
		// Pr(step 2 and step 3) = -------------- = ---------------
		//                           2^n - 2(n+1)   2^(n-1) - n - 1
		//
		// A different way to obtain Pr(step 2 and step 3) is to imagine placing nodes randomly in two piles,
		// pile 1 and pile 2. After the first node is placed, arbitrarily establishing pile 1, there are
		// n - 1 nodes left. For each one let the probability it is placed in pile 1 be 0.5, with the
		// probability of placing the node in pile 2 being 0.5 also. The probability of any one outcome is
		// thus 2^{n-1}, but some outcomes are not acceptable. The n outcomes in which one taxon is separated
		// from all other taxa are not acceptable because they are indistinguishable from the orignal polytomy.
		// In addition, the outcome in which all taxa end up in pile 1 is unacceptable for the same reason.
		// Thus, there are n + 1 unacceptable outcomes, and these must be subtracted from the 2^{n-1} possible
		// outcomes when determining the denominator of the probability.
		//
		// The probability of proposing an add-edge move is thus:
		//
		//                          1                   1             1                     1
		// Pr(add-edge move) = -------------- X ---------------- X ------- = ----------------------------------
		//                     num_polytomies    2^(n-1) - n - 1   (1 - 0)   num_polytomies * [2^(n-1) - n - 1]
		//
		// The delete-edge move is simpler, involving only the choice of the edge that must be deleted from the tree to
		// reverse the add-edge move (probability is inverse of the number of edges in the post-add-edge-move tree).
		//
		// Thus, the Hastings ratio for an add-edge move is:
		//
		//                                                                1
		//                                                  -----------------------------
		//  Pr(delete-edge of proposed add-edge move)       num_internal_edges_before + 1
		//  ----------------------------------------- = -------------------------------------
		//      Pr(proposed add-edge move)                                1
		//                                              -------------------------------------
		//                                               num_polytomies * [2^(n-1) - n - 1]
		//
		//                                                num_polytomies * [2^(n-1) - n - 1]
		//                                            =  -------------------------------------
		//                                                   num_internal_edges_before + 1
		//
		// If proposed state is the fully-resolved tree, or if the reverse move would generate a fully-resolved tree,
		// then the Hastings ratio must account for the fact that one of the moves is only attempted half the time
		// whereas the other move is the only move possible. If the proposal generates a fully-resolved tree,
		// multiply the Hastings ratio as calculated above by 2.0. If the proposal takes us away from a fully-resolved
		// tree, then divide the Hastings ratio as calculated above by 2.0.

		// nspokes is n_p in Lewis, Holder and Holsinger (2005)
		double nspokes = (double)polytomy_size;

		// Compute the log of the Hastings ratio
		ln_hastings = log((double)num_polytomies);
		ln_hastings += log(pow(2.0, nspokes - 1.0) - nspokes - 1.0);
		ln_hastings -= log((double)num_internal_edges_before + 1.0);

		// Now multiply by the value of the quantity labeled gamma_b in the paper
		nnodes = tree->GetNNodes();
		const bool fully_resolved_after = (nnodes == num_nodes_in_fully_resolved_tree);
		if (star_tree_before && !fully_resolved_after)
			ln_hastings -= log(2.0);
		else if (fully_resolved_after && !star_tree_before)
			ln_hastings += log(2.0);
		//@POL the pow will be problematic for large polytomies
		}
	else
		{
		// Choose an internal node at random (but not the only child of the root node)
		// and delete its edge to create a polytomy (or a bigger polytomy if there is
		// already a polytomy)
		//
		unsigned num_internals = nnodes - num_taxa - 1;
		PHYCAS_ASSERT(num_internals > 0);
		PHYCAS_ASSERT(num_internals < num_taxa);
		unsigned i = rng->SampleUInt(num_internals);

		for (nd = tree->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (nd->IsTip() || nd->GetParent()->IsTipRoot())
				continue;

			if (i == 0)
				break;
			else
				--i;
			}
        //nd->SelectNode();
		//likelihood->startTreeViewer(tree, "Delete edge move");
		proposeDeleteEdgeMove(nd);

		// The Jacobian is the matrix of the derivatives of all post-move quantities with respect to all
		// pre-move quantities. In this case, there is only one pre-move quantity (v) and only one post-move
		// quantity (U), so the Jacobian is a 1x1 matrix consisting only of the derivative of U with respect
		// to v. U = 1 - exp{-theta*v}, so the Jacobian for this delete-edge move is theta*exp(-theta*v),
		// where v is the length of the edge being deleted and theta = 1/edgelen_mean is the hazard parameter
		// of the exponential distribution used for sampling new edge lengths. Thus,
		//
		// ln_jacobian = log[theta*exp(-theta*v)]
		//             = log[exp(-v/edgelen_mean)/edgelen_mean]
		//             = -(v/edgelen_mean) - log(edgelen_mean)
		//
		// Note that v = orig_edgelen; orig_edgelen was set in ProposeDeleteEdgeMove.
		ln_jacobian = -orig_edgelen/edgelen_mean - log(edgelen_mean);

		// Hastings ratio is inverse of that for the add-edge move (see extensive notes above), but be
		// careful to use the number of polytomies in tree *after* a delete-edge move and the number of internal
		// edges *before* a delete-edge move in the formula. Here, n is the number of spokes in the polytomy
		// created by the delete-edge move. Both polytomy_size and num_polytomies are correctly computed by
		// ProposeDeleteEdgeMove.
		//
		//  Pr(add-edge move that reverts proposed delete-edge move)         num_internal_edges_before
		//  -------------------------------------------------------- = -----------------------------------------
		//                Pr(proposed delete-edge move)                 num_polytomies_after * [2^(n-1) - n - 1]
		//
		// If the proposed state is the star tree, then the Hastings ratio must account for the fact that
		// the forward move is only attempted half the time whereas the reverse move is the only move possible
		// for the star tree. Thus, if the number of nodes in the tree after the proposal is 1 more than the
		// number of taxa, the Hastings ratio as computed above must be multiplied by 2.0. If the proposal takes
		// us away from the star tree, then the Hastings ratio must be divided by 2.0.
		//
		double nspokes = (double)polytomy_size;

		// Compute the log of the Hastings ratio
		ln_hastings = log((double)num_internal_edges_before);
		ln_hastings -= log((double)num_polytomies);
		ln_hastings -= log(pow(2.0, nspokes - 1.0) - nspokes - 1.0);

		// Now multiply by the value of the quantity labeled gamma_d in the paper
		nnodes = tree->GetNNodes();
		const bool star_tree_after = (nnodes == num_taxa + 1);
		if (fully_resolved_before && !star_tree_after)
			{
			ln_hastings -= log(2.0);
			}
		else if (star_tree_after && !fully_resolved_before)
			{
			ln_hastings += log(2.0);
			}
		//@POL the pow will be problematic for large polytomies
		}

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();
    jpm->skipNextDebugCheck();
    jpm->allEdgeLensModified(tree);
    jpm->topologyModified("tree_topology", tree);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets `num_taxa' to the number of tips in `tree', computes the polytomy distribution, and sets the number of taxa for the
|	`topo_prior_calculator' object to `num_taxa'.
*/
void BushMove::finalize()
	{
	num_taxa = tree->GetNTips();
	num_nodes_in_fully_resolved_tree = 2*(num_taxa - 1);    // ASSUMES UNROOTED TREE

    // 1   2 3   4    1   2   3   4    Fully-resolved rooted trees have 2*num_taxa nodes
    //  \ /   \ /      \   \   \ /     if the root node (fake tip) is included, regardless
    //   5     6        \   \   5      of symmetry.
    //    \   /          \   \ /
    //     \ /            \   6
    //      7              \ /
    //      |               7
    //      8               |
    //                      8

	computePolytomyDistribution(num_taxa);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void BushMove::revert()
	{
	MCMCUpdater::revert();
	if (add_edge_move_proposed)
		{
		// An add-edge move was previously proposed, which means a new internal node (orig_lchild) and its edge were created.
		// orig_lchild's parent is known as orig_par. The children of orig_lchild represent a subset of the original children
		// of orig_par. To reverse the add-edge move, we need only transfer all children of orig_lchild to orig_par, then store
		// orig_lchild for use in some later add-edge move proposal

		// Must discard CLAs before moving children, otherwise orig_lchild will look like a tip node
		// which causes invalidateBothEndsDiscardCache to attempt to access its TipData structure
		likelihood->invalidateBothEndsDiscardCache(orig_lchild); //@POL should not be any cached CLAs here

		// Transfer all children of orig_lchild to orig_par
		while (orig_lchild->GetLeftChild() != NULL)
			{
			TreeNode * s = orig_lchild->GetLeftChild();
			tree_manipulator.DetachSubtree(s);
			tree_manipulator.InsertSubtree(s, orig_par, TreeManip::kOnRight);
			}

		// Now delete orig_lchild
		//
		tree_manipulator.DetachSubtree(orig_lchild);
		tree->StoreInternalNode(orig_lchild);

		likelihood->useAsLikelihoodRoot(orig_par);
		likelihood->restoreFromCacheAwayFromNode(*orig_par);

		//likelihood->startTreeViewer(tree, "Add edge move REVERTED");

		tree->InvalidateNodeCounts();
		}
	else
		{
		// A delete-edge move was proposed, which means an internal node and its edge have been deleted. This deleted
		// node's parent is known as orig_par. The children of the deleted node have been added as children of orig_par,
		// with the first being identified now as orig_lchild and the last in the series now pointed to by orig_rchild.
		// To revert the delete-edge move, create a new child of orig_par having edgelen equal to orig_edgelen, then transfer
		// the nodes from orig_lchild to orig_rchild (following rSib pointers, and including both orig_lchild and
		// orig_rchild) to the new node.
		//
		TreeNode * u = tree->GetNewNode();
		u->SetEdgeLen(orig_edgelen);
		tree_manipulator.InsertSubtree(u, orig_par, TreeManip::kOnRight);

		TreeNode * nd = orig_lchild;
		for (;;)
			{
			TreeNode * s = nd;
			PHYCAS_ASSERT(s != NULL);
			nd = nd->GetRightSib();
			tree_manipulator.SibToChild(u, s, TreeManip::kOnRight);
			if (s == orig_rchild)
				break;
			}

		likelihood->useAsLikelihoodRoot(orig_par);
		likelihood->restoreFromCacheAwayFromNode(*orig_lchild);
		likelihood->restoreFromCacheParentalOnly(orig_lchild);

		//orig_par->UnselectNode();
		//likelihood->startTreeViewer(tree, "Delete edge move REVERTED");

		tree->InvalidateNodeCounts();
		}

	reset();

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();
    jpm->skipNextDebugCheck();
    jpm->allEdgeLensModified(tree);
    jpm->topologyModified("tree_topology", tree);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears the `polytomies' vector, then walks the tree adding nodes that represent polytomies (more than two children)
|	to the `polytomies' vector.
*/
void BushMove::refreshPolytomies()
	{
	//@POL should keep polytomies list up to date rather than building anew
	polytomies.clear();
	for (TreeNode *nd = tree->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTipRoot())
			continue;

		unsigned s = nd->CountChildren();
		if (s > 2)
			{
			polytomies.push_back(nd);
			}
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Determines distribution of x given nspokes, where x is the number spokes assigned to the newly created node in an
|	add-edge move. The number of ways of choosing x spokes to move out of nspokes total is {nspokes \choose x}. We are not
|	interested in the values x = 0, x = 1, x = nspokes - 1, and x = n because these lead either to a non-move (x = 0 and
|	x = nspokes) or to a tree that has an invalid structure (x = 1 and x = nspokes - 1). Thus, the total number of possible
|	trees considered is {nspokes \choose 2} + {nspokes \choose 3} + ... + {n \choose nspokes - 2} =
|	2^nspokes - 2*(nspokes + 1). The 2^nspokes comes from the fact that 2^nspokes is the sum of all binomial coefficients
|	for a sample size of nspokes. The subtracted term 2(nspokes + 1) comes from the fact that the first and last binomial
|	coefficients - {nspokes \choose 0} and {nspokes \choose nspokes} - are always 1 and the second and penultimate binomial
|	coefficients - {nspokes \choose 1} and {nspokes \choose nspokes - 1} - always equal nspokes. Thus, if one wishes to
|	choose randomly from all possible ways of splitting the polytomy into two groups of spokes, select x with probability:
|>
|	                    (nspokes \choose x}
|	  Pr(X = x) = -----------------------------, x = 2, 3, ..., nspokes - 2
|	                2^nspokes - 2*(nspokes + 1)
|>
*/
const VecPolytomyDistr & BushMove::computePolytomyDistribution(unsigned nspokes)
	{
	PHYCAS_ASSERT(nspokes > 2);
	std::pair<PolytomyDistrMap::const_iterator, bool> retval;
	PolytomyDistrMap::const_iterator i = poly_prob.find(nspokes);
	if (i == poly_prob.end())
		{
		// There is no existing probability distribution vector corresponding to nspokes
		// Need to calcuate it and insert into poly_prob map.
		double ln_nfact = cdf.LnGamma((double)(nspokes + 1));
		double denom = exp((double)nspokes * log(2.0)) - 2.0*nspokes - 2.0; //@POL may need to factor if n large
		double ln_denom = log(denom);
		VecPolytomyDistr v;
		for (unsigned x = 2; x <= nspokes - 2; ++x)
			{
			double ln_numer = ln_nfact - cdf.LnGamma((double)(x + 1)) - cdf.LnGamma((double)(nspokes - x + 1));
			double prob_x = exp(ln_numer - ln_denom);
			v.push_back(prob_x);
			}
		retval = poly_prob.insert(PolytomyDistrMap::value_type(nspokes, v));
		PHYCAS_ASSERT(retval.second == true);
		i = retval.first;
		}
	return (i->second);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Split up the polytomy at `u' by creating a new internal node v and a new edge connecting u with v. Node u is saved
|	as `orig_par' and node v is saved as `orig_lchild' in case we need to revert the proposed move.
*/
void BushMove::proposeAddEdgeMove(TreeNode * u)
	{
	PHYCAS_ASSERT(u != NULL);
	polytomy_size = 1 + u->CountChildren();

	const std::vector<double> & prob_n = computePolytomyDistribution(polytomy_size);

	// Select number of spokes to move over to new node
	unsigned x;
	double p = rng->Uniform();
	double cum = 0.0;
	for (unsigned k = 0; k <= polytomy_size - 4; ++k)
		{
		x = k + 2;
		double prob_k_given_n = prob_n[k];
		cum += prob_k_given_n;
		if (p < cum)
			break;
		}
	PHYCAS_ASSERT(x < polytomy_size - 1);

	// Create the new node that will receive the x randomly-chosen spokes
	new_edgelen = edgelen_dist->Sample();
	TreeNode * v = tree->GetNewNode();
	v->SetEdgeLen(new_edgelen);
	v->SetNodeNum(tree->GetNNodes());
	likelihood->prepareInternalNodeForLikelihood(v);
	tree_manipulator.InsertSubtree(v, u, TreeManip::kOnLeft);

	// Save u and v. If revert is necessary, all of orig_lchild's nodes will be returned
	// to orig_par, and orig_lchild will be deleted.
	//
	orig_par = u;
	orig_lchild = v;

	// After the move, either v or u should have x spokes and the other node polytomy_size - x spokes (u and v will
	// each have 1 additional connector spoke).Choose x spokes randomly out of the polytomy_size available.
	// If u->par is included, let u retain the x spokes and move polytomy_size - x spokes to v. Otherwise, move the
	// x spokes to v leaving polytomy_size - x spokes behind.
	//
	std::vector<TreeNode *> uspokes;
	uspokes.push_back(u->GetParent());
	for (TreeNode *uchild = u->GetLeftChild(); uchild != NULL; uchild = uchild->GetRightSib())
		{
		if (uchild != v)
			uspokes.push_back(uchild);
		}
	PHYCAS_ASSERT(uspokes.size() == polytomy_size);

	bool reverse_polarity = false;
	std::vector<TreeNode *> vspokes;
	typedef std::vector<TreeNode *>::iterator::difference_type vec_it_diff;
	for (unsigned k = 0; k < x; ++k)
		{
		unsigned num_u_spokes = (unsigned)uspokes.size();
		PHYCAS_ASSERT(num_u_spokes > 0);
		unsigned j = rng->SampleUInt(num_u_spokes);
		TreeNode * s = uspokes[j];
		if (s == u->GetParent())
			reverse_polarity = true;
		vspokes.push_back(s);
		uspokes.erase(uspokes.begin() + (vec_it_diff) j);
		}
	PHYCAS_ASSERT(uspokes.size() + vspokes.size() == polytomy_size);

	if (reverse_polarity)
		{
		// transfer nodes in uspokes to v
		//
		std::vector<TreeNode *>::iterator s;
		for (s = uspokes.begin(); s != uspokes.end(); ++s)
			{
			tree_manipulator.DetachSubtree(*s);
			tree_manipulator.InsertSubtree(*s, v, TreeManip::kOnRight);
			}
		}
	else
		{
		// transfer nodes in vspokes to v
		//
		std::vector<TreeNode *>::iterator s;
		for (s = vspokes.begin(); s != vspokes.end(); ++s)
			{
			tree_manipulator.DetachSubtree(*s);
			tree_manipulator.InsertSubtree(*s, v, TreeManip::kOnRight);
			}
		}

	//orig_lchild->SelectNode();

    //v->SelectNode();
    //likelihood->startTreeViewer(tree, "proposeAddEdgeMove");

	likelihood->useAsLikelihoodRoot(orig_lchild);
	likelihood->invalidateAwayFromNode(*orig_lchild);
	likelihood->invalidateBothEnds(orig_lchild);	//@POL really just need invalidateParentalOnly function

	tree->InvalidateNodeCounts();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Delete the edge associated with `u' to create a polytomy (or a bigger polytomy if `u->par' was already a polytomy).
|	The supplied node u should not be the only child of the root node.
|>
|	      b       c
|	       \     /
|	        \   /
|	         \ /
|	  a       u		   a   b   c
|	   \     /		    \  |  /
|	    \   /		     \ | /
|	     \ /		      \|/
|         v                v
|	     /			      /
|
|	    Before           After
|>
|	Returns the number of polytomies in the tree after the proposed delete-edge move. The return value will be incorrect if
|	the polytomies vector is not up-to-date.
*/
void BushMove::proposeDeleteEdgeMove(TreeNode * u)
	{
	// Node u will be stored, so get rid of any CLAs
	//@POL should not discard cache here in case move is reverted - right now we're discarding because
	// tree_manipulator.DeleteLeaf stores node u, and we do not want any CLAs to travel to node storage
	// and come back to haunt us later. Better solution would be to keep u around until we know the move
	// will be accepted, at which point the CLAs can definitely be discarded
	likelihood->invalidateBothEndsDiscardCache(u);

    bool u_is_internal = u->IsInternal();

    // Save nd's edge length in case we need to revert
	//
	orig_edgelen = u->GetEdgeLen();

	// This operation should not leave the root node (which is a tip) with more than
	// one child, so check to make sure that the supplied node is not the root nor a
	// child of root
	//
	orig_par = u->GetParent();
	PHYCAS_ASSERT(orig_par != NULL);
	PHYCAS_ASSERT(!orig_par->IsTipRoot());

	num_polytomies = (unsigned)polytomies.size();

	// Compute size of polytomy after the delete-edge move, a quantity that is needed for computing the Hastings ratio.
	// Note that one of v's children (i.e. u) is deleted but this is made up for by considering v->par, which is
	// also a spoke that counts.
	//
	unsigned u_children = u->CountChildren();
	unsigned v_children = orig_par->CountChildren();
	polytomy_size = v_children + u_children;

	bool u_polytomy_before = (u_children > 2);
	bool v_polytomy_before = (v_children > 2);
	if (u_polytomy_before && v_polytomy_before)
		{
		// No. polytomies will decrease by one as a result of this delete-edge move
		//
		--num_polytomies;
		}
	else if (!u_polytomy_before && !v_polytomy_before)
		{
		// No. polytomies will increase by one as a result of this delete-edge move
		//
		++num_polytomies;
		}

	// Make all of u's children left siblings (i.e. children of u->par)
	//
	orig_lchild = u->GetLeftChild();
	while (u->GetLeftChild() != NULL)
		{
		orig_rchild = u->GetLeftChild();
		tree_manipulator.DetachSubtree(orig_rchild);
		tree_manipulator.InsertSubtree(orig_rchild, orig_par, TreeManip::kOnRight);
		}

	tree_manipulator.DeleteLeaf(u, u_is_internal);

	//orig_par->SelectNode();

	likelihood->useAsLikelihoodRoot(orig_par);
	likelihood->invalidateAwayFromNode(*orig_par);

	tree->InvalidateNodeCounts();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Setup `edgelen_dist' data member.
*/
void BushMove::setEdgeLenDistMean(
  double mean)
	{
	PHYCAS_ASSERT(mean > 0.0);
	if (rng)
		{
		edgelen_mean = mean;
		edgelen_dist = ExponentialDistributionShPtr(new phycas::ExponentialDistribution(1.0/edgelen_mean));
		edgelen_dist->SetMeanAndVariance(edgelen_mean, edgelen_mean);
		edgelen_dist->SetLot(rng.get());  //@POL should just be able to pass in the shared_ptr rather than a raw pointer
		}
	else
		{
		throw XLikelihood("Must set the pseudorandom number generator before calling setEdgeLenDistMean");
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void BushMove::accept()
	{
	MCMCUpdater::accept();
	if (add_edge_move_proposed)
		{
		// Keeping added edge, where orig_par was original polytomous node and orig_lchild was the added node
		likelihood->useAsLikelihoodRoot(orig_lchild);
		likelihood->discardCacheAwayFromNode(*orig_lchild);
		likelihood->discardCacheBothEnds(orig_lchild);

		//likelihood->startTreeViewer(tree, "Add edge move ACCEPTED");
		}
	else
		{
		// Keeping edge deletion, orig_par is the new polytomous node
		likelihood->useAsLikelihoodRoot(orig_par);
		likelihood->discardCacheAwayFromNode(*orig_par);

		//likelihood->startTreeViewer(tree, "Delete edge move ACCEPTED");
		}

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of `add_edge_move_proposed' data member, which is true if the last move proposed was an
|	add-edge move and false if last move proposed was a delete-edge move.
*/
bool BushMove::addEdgeMoveProposed() const
	{
	return add_edge_move_proposed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Re-initializes all data members that are recomputed each time ProposeNewState is called.
*/
void BushMove::reset()
	{
	// Workspace used for computing edge length prior
	// Should always be length one.
	//if (one_edgelen.empty())
    //    one_edgelen.push_back(0.0);

	polytomies.clear();

	add_edge_move_proposed	= false;
	orig_edgelen			= 0.0;
	orig_lchild				= NULL;
	orig_rchild				= NULL;
	orig_par				= NULL;
	ln_jacobian				= 0.0;
	ln_hastings				= 0.0;
	new_edgelen				= 0.0;
	polytomy_size			= 0;
	num_polytomies			= 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_hastings', which is the natural log of the Hastings ratio for this move and which is
|	computed in both BushMove::ProposeAddEdgeMove and BushMove::ProposeDeleteEdgeMove.
*/
double BushMove::getLnHastingsRatio() const
	{
	return ln_hastings;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_jacobian', which is the natural log of the Jacobian for this move and which is computed in
|	both BushMove::ProposeAddEdgeMove and BushMove::ProposeDeleteEdgeMove.
*/
double BushMove::getLnJacobian() const
	{
	return ln_jacobian;
	}
