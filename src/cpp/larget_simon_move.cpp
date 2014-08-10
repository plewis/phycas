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

#include "probability_distribution.hpp"
#include "model.hpp"
#include "basic_tree_node.hpp"
#include "tree_likelihood.hpp"
#include "xlikelihood.hpp"
#include "mcmc_chain_manager.hpp"
#include "larget_simon_move.hpp"
#include "basic_tree.hpp"

#include "boost/format.hpp"

#define KEEP_BRANCHES_LEGAL
#define MAX_LEGAL_BRLEN 285

namespace phycas
{
bool verbose = false;

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor sets `lambda' to the default value (0.2), sets `topol_changed' to false, and `m' and `mstar'
|	to 0.0. All other data members are automatically initialized (shared pointers) or are initialized via a call to
|	reset().
*/
LargetSimonMove::LargetSimonMove() : MCMCUpdater()
	{
	is_move 		= true;
	topol_changed	= false;
	lambda			= 0.2;
	m				= 0.0;
	mstar			= 0.0;
	three_edgelens.reserve(3);
	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Virtual destructor.
*/
LargetSimonMove::~LargetSimonMove()
	{
	//std::cerr << "\n>>>>> LargetSimonMove dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the probability of the tree under Mark Holder's tree topology reference distribution.
*/
double LargetSimonMove::recalcWorkingPrior() const
	{
	PHYCAS_ASSERT(topo_prob_calc);
	std::pair<double, double> treeprobs = topo_prob_calc->CalcTopologyLnProb(*tree, true);
	const double ln_ref_topo = treeprobs.first;
	const double ln_ref_edges = treeprobs.second;
	double ln_ref_dist = ln_ref_topo + ln_ref_edges;
	return ln_ref_dist;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Samples a tree from Mark's reference distribution.
*/
double LargetSimonMove::sampleWorkingPrior() const
	{
	PHYCAS_ASSERT(topo_prob_calc);
	topo_prob_calc->SampleTree(tree, rng);
	return 0.0;
	}

static bool warned = false;
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true, overriding the base class (MCMCUpdater) version.
*/
bool LargetSimonMove::isPriorSteward() const
	{
	//@MTH@ Warning future bug!!!
	if (!warned) {
	    std::cerr << "\n\nWARNING: HACKY workaround in LargetSimonMove::isPriorSteward() - always returning true!!!\n\n\n\n";
	    warned = true;
	}
	return use_ref_dist; //bool(topo_prob_calc);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or
|	revert(), whichever is appropriate.
*/
bool LargetSimonMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed LargetSimonMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

    tree->renumberInternalNodes(tree->GetNTips()); //@POL this should be somewhere else
	tree->RecalcAllSplits(tree->GetNTips());

    ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

    // compute the likelihood before proposing a new state
	double prev_ln_like = p->getLastLnLike();

	double prev_ln_prior = jpm->getLogJointPrior();

	TreeNode * prev_likelihood_root = likelihood->getLikelihoodRoot();
	double prev_ln_ref_dist = 0.0;
	if (use_ref_dist)
		{
        prev_ln_ref_dist = recalcWorkingPrior();
		}

	//likelihood->startTreeViewer(tree, boost::str(boost::format("LS move BEFORE: %.5f") % prev_ln_like));

	proposeNewState();

    curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);

	//likelihood->startTreeViewer(tree, boost::str(boost::format("LS move AFTER: %.5f") % curr_ln_like));

	curr_ln_prior = jpm->getLogJointPrior();

	double curr_ln_ref_dist = 0.0;
    if (likelihood->getTreeLengthPrior())
        {
        PHYCAS_ASSERT(!use_ref_dist); // not ready for this yet
        curr_ln_prior = likelihood->getTreeLengthPrior()->GetLnPDF(tree);
        }
    else
        {
        if (use_ref_dist)
            {
            curr_ln_ref_dist = recalcWorkingPrior();
            }
        }

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;

	if (is_standard_heating)
		{
		prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
		curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
		if (use_ref_dist)
			{
			prev_posterior += (1.0 - heating_power)*prev_ln_ref_dist;
			curr_posterior += (1.0 - heating_power)*curr_ln_ref_dist;
			}
		}
	else
		{
		prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
		curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
		}

	double ln_accept_ratio = curr_posterior - prev_posterior + getLnHastingsRatio() + getLnJacobian();

    double lnu = DBL_MAX;
    bool accepted = (ln_accept_ratio >= 0.0);
    if (!accepted)
        {
        double u = rng->Uniform();
        lnu = std::log(u);
        accepted = (lnu <= ln_accept_ratio);
        }

    if (save_debug_info)
        {
    	if (star_tree_proposal)
            {
            debug_info = str(boost::format("LS: %.5f -> %.5f (%s)") % orig_edge_len % orig_node->GetEdgeLen() % (accepted ? "accepted" : "rejected"));
            }
        else
            {
            debug_info = boost::str(boost::format("%s, prev_ln_like = %.5f, getLastLnLike() = %.5f, curr_ln_like = %.5f, topology %s, case = %d, x=%f, y=%f, z=%f, newX=%f, newY=%f, newZ=%f, lnu = %.5f, lnr = %.5f, curr = %.5f, prev = %.5f")
                % (accepted ? "ACCEPT" : "REJECT")
                % prev_ln_like
                % p->getLastLnLike()
                % curr_ln_like
                % (topol_changed ? "changed" : "unchanged")
                % which_case
                % x
                % y
                % z
                % (ndX->GetEdgeLen())
                % (ndY->GetEdgeLen())
                % (ndZ->GetEdgeLen())
                % (lnu == DBL_MAX ? -1.0 : lnu)
                % ln_accept_ratio
                % curr_posterior
                % prev_posterior);

            if (is_standard_heating)
                {
                if (use_ref_dist)
                    {
                    debug_info += boost::str(boost::format("\n  prev_posterior = %g = %g*(%g + %g) + (1.0 - %g)*%g")
                        % prev_posterior
                        % heating_power
                        % prev_ln_like
                        % prev_ln_prior
                        % heating_power
                        % prev_ln_ref_dist
                        );
                    debug_info += boost::str(boost::format("\n  curr_posterior = %g = %g*(%g + %g) + (1.0 - %g)*%g")
                        % curr_posterior
                        % heating_power
                        % curr_ln_like
                        % curr_ln_prior
                        % heating_power
                        % curr_ln_ref_dist
                        );
                    }
                else
                    {
                    debug_info += boost::str(boost::format("\n  prev_posterior = %g = %g*(%g + %g)")
                        % prev_posterior
                        % heating_power
                        % prev_ln_like
                        % prev_ln_prior
                        );
                    debug_info += boost::str(boost::format("\n  curr_posterior = %g = %g*(%g + %g)")
                        % curr_posterior
                        % heating_power
                        % curr_ln_like
                        % curr_ln_prior
                        );
                    }
                }
            else
                {
                debug_info += boost::str(boost::format("\n  prev_posterior = %g = %g*%g + %g")
                    % prev_posterior
                    % heating_power
                    % prev_ln_like
                    % prev_ln_prior
                    );
                debug_info += boost::str(boost::format("\n  curr_posterior = %g = %g*%g + %g")
                    % curr_posterior
                    % heating_power
                    % curr_ln_like
                    % curr_ln_prior
                    );
                }


            if (!prev_likelihood_root)
                debug_info += "\n  prev_likelihood_root = NULL";
            else
                debug_info += boost::str(boost::format("\n  prev_likelihood_root = %g") % prev_likelihood_root->GetNodeNumber());

            TreeNode * tmp_curr_likelihood_root = likelihood->getLikelihoodRoot();
            if (!tmp_curr_likelihood_root)
                debug_info += "\n  curr_likelihood_root = NULL";
            else
                debug_info += boost::str(boost::format("\n  curr_likelihood_root = %g") % tmp_curr_likelihood_root->GetNodeNumber());
            }
        }

    if (accepted)
		{
		p->setLastLnLike(curr_ln_like);
        accept();
		}
	else
		{
		curr_ln_like = p->getLastLnLike();
		revert();

		PHYCAS_ASSERT(!prev_likelihood_root || prev_likelihood_root->IsInternal());
		likelihood->useAsLikelihoodRoot(prev_likelihood_root);
		}

    //POLTMP
    lambda = p->adaptUpdater(lambda, nattempts, accepted);
    //std::cerr << boost::str(boost::format("~~~> log(lambda) = %.5f <~~~") % log(lambda)) << std::endl;

	return accepted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If the tree has only one internal node (i.e. it is the star tree), then `star_tree_proosal' is set to true and the
|	starTreeProposeNewState function is called. If the tree is not the star tree, then `star_tree_proosal' is set to
|	false and the defaultProposeNewState function is called.
*/
void LargetSimonMove::proposeNewState()
	{
	if (tree->GetNInternals() == 1)
		{
		starTreeProposeNewState();
		star_tree_proposal = true;
		}
	else
		{
		defaultProposeNewState();
		star_tree_proposal = false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Chooses a random edge and changes its current length m to a new length m* using the following formula, where `lambda' is
|	a tuning parameter.
|>
|	m* = m*exp(lambda*(r.Uniform() - 0.5))
|>
*/
void LargetSimonMove::starTreeProposeNewState()
	{
	// Choose edge randomly.
	//
	unsigned numEdges = tree->GetNNodes() - 1;
	unsigned k = rng->SampleUInt(numEdges);
	unsigned i = 0;
	//@POL this loop is crying out for the for_each algorithm
	for (orig_node = tree->GetFirstPreorder(); orig_node != NULL; orig_node = orig_node->GetNextPreorder())
		{
		// All nodes have an edge associated with them except for the root
		//
		if (!orig_node->IsTipRoot())
			{
			if (i == k)
				{
				orig_edge_len = orig_node->GetEdgeLen();
				break;
				}
			++i;
			}
		}

	// Modify the edge
	//
	double m		= orig_node->GetEdgeLen();
	double mstar	= m*std::exp(lambda*(rng->Uniform() - 0.5));
	orig_node->SetEdgeLen(mstar);

	// Invalidate CLAs to ensure next likelihood calculation will be correct
	orig_node->SelectNode();
	TreeNode * nd = orig_node->IsTip() ? orig_node->GetParent() : orig_node;
	PHYCAS_ASSERT(nd->IsInternal());
	likelihood->useAsLikelihoodRoot(nd);
	likelihood->invalidateAwayFromNode(*orig_node);
	likelihood->invalidateBothEnds(orig_node);

    ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();
    jpm->allEdgeLensModified(tree);
    //jpm->externalEdgeLensModified("external_edgelen", tree);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Selects an internal node at random from a discrete uniform distribution with the constraint that the returned node
|   is not equal to the subroot (the sole child of the tip node serving as the root).
*/
TreeNode * LargetSimonMove::randomInternalAboveSubroot()
    {
	// Avoiding the "subroot" node (only child of the tip serving as the root), so the number of
	// acceptable nodes is one fewer than the number of internal nodes
	unsigned numAcceptableNodes = tree->GetNInternals() - 1;

	unsigned ypos = rng->SampleUInt(numAcceptableNodes);
	unsigned i = 0;
    TreeNode * nd = tree->GetFirstPreorder();
	for (; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsInternal() && !nd->GetParentConst()->IsTipRoot())
			{
			if (i == ypos)
				break;
			++i;
			}
		}
	PHYCAS_ASSERT(nd->GetLeftChild() != NULL);
	PHYCAS_ASSERT(nd->GetParentConst() != NULL);
	PHYCAS_ASSERT(!nd->GetParent()->IsTipRoot());
    return nd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Selects a child of the supplied `nd' at random from a discrete uniform distribution.
*/
TreeNode * LargetSimonMove::randomChild(
  TreeNode * nd)    /**< is the parent node whose children are candidates */
    {
	unsigned ychildren = nd->CountChildren();
	unsigned which_child = rng->SampleUInt(ychildren);
	unsigned k = 0;
    TreeNode * child = nd->GetLeftChild();
	for (; child != NULL; child = child->GetRightSib())
		{
		if (k == which_child)
			break;
		++k;
		}
	PHYCAS_ASSERT(child != NULL);
    return child;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Randomly chooses a node to serve as node Z (the bottom of the three nodes involved in a Larget-Simon move). The
|   supplied node `middle' is the node serving as Y. In the figure below, the nodes labeled Z are all possible
|   candidates for the return value of this function. The node selected as Z should be the owner of the lowermost edge
|   involved in the move (X owns the uppermost edge and Y owns the middle edge).
|>
|	     X  X  X
|	      \ | /
|          \|/
|	  Z  Z  Y
|	   \ | /
|	    \|/
|	     Z
|	     |
|>
*/
TreeNode * LargetSimonMove::chooseZ(
  TreeNode * middle)    /**< is the middle node (Y) */
    {
    TreeNode * nd = NULL;
	TreeNode * U = middle->GetParent();
	PHYCAS_ASSERT(U != NULL);
	unsigned uchildren = U->CountChildren();
	unsigned which_child = rng->SampleUInt(uchildren);
	if (which_child == 0)
		{
		// Selected "child" is actually U's parent
		nd = U;
		}
	else
		{
		// Selected child is one of U's actual children (but cannot be equal to middle)
		unsigned k = 1;
		for (nd = U->GetLeftChild(); nd != NULL; nd = nd->GetRightSib())
			{
			if (nd == middle)
				continue;
			else
				{
				if (k == which_child)
					break;
				++k;
				}
			}
		PHYCAS_ASSERT(nd != NULL);
		}
    return nd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a local perturbation using a generalization of the algorithm "LOCAL Without a Molecular Clock" described by
|	Larget and Simon (1999. Mol. Biol. Evol. 16(6): 750-759). This version allows polytomies (except the case of the
|	star tree, for which the required three-contiguous-edge segment cannot be identified).
|
|	     X  c  d
|	      \ | /
|          \|/
|	  a  b  Y
|	   \ | /
|	    \|/
|	     u
|	     |
|	     Z <-- may or may not be the tip at which tree is rooted
|
|	Pick a random interior node Y whose parent is not the root (i.e. avoid the subroot node directly connected to the
|	tip node at which the tree is rooted. Let u be the parent of y. Let Z be a randomly-chosen child of u (note that in
|	this case u's parent is considered a "child" of u). In the figure above, we (by chance) chose the parent of u to be
|	Z, but we could have chosen any of u's real children (except Y). a and b are the other "children" of u, not
|	including Y. Let X be a randomly chosen child of Y (and here a "child" is really a child). c and d are the other
|	children of Y.
|
|	      a   b   c   d
|	       \ /     \ /
|	Z ===== u ===== Y ===== X
|
|	The path represented by the double line above is either contracted or expanded by a factor m*, where
|
|	m* = m*exp(lambda*(r.Uniform() - 0.5))
|
|	Then, one of {u, Y} is chosen at random to move. Let's say for illustration that u was chosen. u is moved (along
|	with a and c) to a random point along the main path from Z to X. If this makes the node u cross over node Y, then
|	the equivalent of an NNI rearrangement is effected:
|
|	    X  c  d          X  a  b
|	     \ | /            \ | /    In this case, invalidate CLAs away from u and
|         \|/              \|/     make u the likelihood root
|	 a  b  Y          c  d  u
|	  \ | /   -->      \ | /
|	   \|/              \|/
|	    u                Y
|	    |                |
|	    Z                Z
|
|	If there is no NNI rearrangement, the move will only require adjusting edge lengths. In this case, invalidate CLAs
|	away from Y and make Y the likelihood root.
*/
void LargetSimonMove::defaultProposeNewState()
	{
	double xstar, ystar, zstar;

    ChainManagerShPtr p = chain_mgr.lock();
    PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

	// Make sure all the necessary shared pointers have been set to something meaningful
	PHYCAS_ASSERT(rng);
	PHYCAS_ASSERT(tree);
	PHYCAS_ASSERT(model);
	PHYCAS_ASSERT(likelihood);

	// Begin by resetting all the data members involved with reverting a move
	reset();

	// Select an internal node whose parent is not the root node to serve as ndY,
    // whose branch will form the middle segment in the path of three contiguous
    // segments to be modified.
    ndY = randomInternalAboveSubroot();
	y = ndY->GetEdgeLen();

	// Set ndX equal to a randomly-chosen child of ndY
	ndX = randomChild(ndY);
	x = ndX->GetEdgeLen();

	// Set ndZ randomly to either the parent of ndY or one of ndY's sibs
    ndZ = chooseZ(ndY);
    z = ndZ->GetEdgeLen();

    // Set ndBase to the deepest affected node
	//ndBase = ndZ->GetParent();

    // Set node U to the other node (besides ndY) that could slide.
    // Note that U may or may not be equal to ndZ.
    TreeNode * ndU = ndY->GetParent();

	m = x + y + z;
    expand_contract_factor = exp(lambda*(rng->Uniform() - 0.5));
	mstar = m*expand_contract_factor;

	xstar = x*expand_contract_factor;
	ndX->SetEdgeLen(xstar);

    ystar = y*expand_contract_factor;
	ndY->SetEdgeLen(ystar);

    zstar = z*expand_contract_factor;
	ndZ->SetEdgeLen(zstar);

	double d = rng->Uniform()*mstar;

	// Decide whether to move ndY a distance d from the top, or
    // move ndU a distance d from the bottom
	bool moving_Y = true;
	if (rng->Uniform() < 0.5)
		moving_Y = false;
	bool moving_U = !moving_Y;

    // Determine whether proposed move will change the topology
    topol_changed = false;
    if (moving_Y && d >= xstar + ystar)
        topol_changed = true;
    if (moving_U && d >= ystar + zstar)
        topol_changed = true;

	if (topol_changed)
		{
        jpm->allEdgeLensModified(tree);
        jpm->topologyModified("tree_topology", tree);
        //jpm->externalEdgeLensModified("external_edgelen", tree);
        //jpm->internalEdgeLensModified("internal_edgelen", tree);
        if (moving_Y)
            {
            double f = (d - xstar - ystar)/zstar;
            if (ndU == ndZ)
                {
                which_case = 6;

                //ndX->SelectNode();
                //ndY->SelectNode();
                //ndZ->SelectNode();
                //std::string titlestr = str(boost::format("%s, %s: f = %.8f, d = %.8f, r = %.8f, x = %.8f, y = %.8f, z = %.8f") % (moving_Y ? "Y" : "U") % (topol_changed ? "nni" : "nul") % f % d % expand_contract_factor % x % y % z);
                //likelihood->startTreeViewer(tree, titlestr);

                likelihood->slideNode(-1.0, ndY, ndX);  // move all of Y's edge to base of X

                //titlestr = str(boost::format("Before NNISwapSpecial: x = %.8f, y = %.8f, z = %.8f") % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                //likelihood->startTreeViewer(tree, titlestr);

                swap1 = ndX;
			    swap2 = ndZ->GetParent();
			    tree_manipulator.NNISwapSpecial(swap1); // moves Y below Z
                likelihood->swapInternalDataAndEdgeLen(ndY, ndZ);

                //titlestr = str(boost::format("After NNISwapSpecial: x = %.8f, y = %.8f, z = %.8f") % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                //likelihood->startTreeViewer(tree, titlestr);

                likelihood->slideNode(1.0 - f, ndZ, ndY);

                //titlestr = str(boost::format("After final slideNode: x = %.8f, y = %.8f, z = %.8f") % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                //likelihood->startTreeViewer(tree, titlestr);
                }
            else
                {
                which_case = 3;
                likelihood->slideNode(-1.0, ndY, ndX);  // move all of Y's edge to base of X
    			swap1 = ndX;
	    		swap2 = ndZ;
    			tree_manipulator.NNISwap(swap1, swap2);  // nearest-neighbor interchange
                likelihood->slideNode(f, ndY, ndZ);  // move a fraction f of the basal part of Z's edge onto Y
                }
            }
        else
            {
            double f = (d - zstar - ystar)/xstar;   // positive f means slide ndU up past ndY toward ndX
            if (ndU == ndZ)
                {
                which_case = 12;

                //ndX->SelectNode();
                //ndY->SelectNode();
                //ndZ->SelectNode();
                //std::string titlestr = str(boost::format("%s, %s: f = %.8f, d = %.8f, r = %.8f, x = %.8f, y = %.8f, z = %.8f") % (moving_Y ? "Y" : "U") % (topol_changed ? "nni" : "nul") % f % d % expand_contract_factor % x % y % z);
                //likelihood->startTreeViewer(tree, titlestr);

                likelihood->slideNode(-1.0, ndZ, ndY);  // move all of Z's edge to base of Y
			    swap1 = ndX;
			    swap2 = ndZ->GetParent();
			    tree_manipulator.NNISwapSpecial(swap1); // moves Y below Z
                likelihood->swapInternalDataAndEdgeLen(ndY, ndZ);

                //titlestr = str(boost::format("After NNISwapSpecial: x = %.8f, y = %.8f, z = %.8f") % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                //likelihood->startTreeViewer(tree, titlestr);

                likelihood->slideNode(f, ndY, ndX); // move a fraction f of the basal part of X's edge onto Y

                //titlestr = str(boost::format("After slideNode: x = %.8f, y = %.8f, z = %.8f") % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                //likelihood->startTreeViewer(tree, titlestr);
                }
            else
                {
                which_case = 9;
                likelihood->slideNode(-1.0, ndY, ndZ);  // move all of Y's edge to base of Z
                likelihood->slideNode(f, ndY, ndX);     // move a fraction f of the basal part of X's edge onto Y
    			swap1 = ndX;
	    		swap2 = ndZ;
    			tree_manipulator.NNISwap(swap1, swap2);  // nearest-neighbor interchange
                }
            }
        }
    else    // topol_changed
        {
        jpm->allEdgeLensModified(tree);
        //jpm->externalEdgeLensModified("external_edgelen", tree);
        //jpm->internalEdgeLensModified("internal_edgelen", tree);
        if (moving_Y)
            {
            if (d < xstar)
                {
                // sliding Y up
                if (ndU == ndZ)
                    which_case = 4;
                else
                    which_case = 1;
                double f = (xstar - d)/xstar;   // positive f means slide ndY up
                likelihood->slideNode(f, ndY, ndX); // slide ndY up a fraction f into ndX's edge
                }
            else
                {
                // sliding Y down
                if (ndU == ndZ)
                    which_case = 5;
                else
                    which_case = 2;
                double f = -(d - xstar)/ystar;  // negative f means slide ndY down
                likelihood->slideNode(f, ndY, ndX); // slide ndY down a fraction f into its own edge
                }
            }
        else // moving U
            {
            if (d < zstar)
                {
                // sliding U down
                if (ndU == ndZ)
                    {
                    which_case = 11;
                    double f = -(zstar - d)/zstar;   // negative f means slide ndZ down
                    likelihood->slideNode(f, ndZ, ndY); // slide ndZ down a fraction f into its own edge
                    }
                else
                    {
                    which_case = 8;
                    double f = -(zstar - d)/zstar;  // negative f means slide ndU down
                    likelihood->slideNode(f, ndZ, ndY); // slide ndU "down" a fraction f into ndZ's edge
                    }
                }
            else
                {
                // sliding U up
                if (ndU == ndZ)
                    {
                    which_case = 10;
                    double f = (d - zstar)/ystar;   // positive f means slide ndZ up
                    likelihood->slideNode(f, ndZ, ndY); // slide ndZ up a fraction f into ndY's edge
                    }
                else
                    {
                    which_case = 7;
                    double f = (d - zstar)/ystar;   // positive f means slide ndU up

                    //ndX->SelectNode();
                    //ndY->SelectNode();
                    //ndZ->SelectNode();
                    //std::string titlestr = str(boost::format("Before slideNode: %s, %s: f = %.8f, d = %.8f, r = %.8f, x = %.8f, y = %.8f, z = %.8f") % (moving_Y ? "Y" : "U") % (topol_changed ? "nni" : "nul") % f % d % expand_contract_factor % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                    //likelihood->startTreeViewer(tree, titlestr);

                    likelihood->slideNode(f, ndZ, ndY); // slide ndU up a fraction f into ndY's edge

                    //titlestr = str(boost::format("After slideNode: %s, %s: f = %.8f, d = %.8f, r = %.8f, x = %.8f, y = %.8f, z = %.8f") % (moving_Y ? "Y" : "U") % (topol_changed ? "nni" : "nul") % f % d % expand_contract_factor % ndX->GetEdgeLen() % ndY->GetEdgeLen() % ndZ->GetEdgeLen());
                    //likelihood->startTreeViewer(tree, titlestr);
                    }
                }
            }
        }

	ndX->SelectNode();
	ndY->SelectNode();
	ndZ->SelectNode();
	PHYCAS_ASSERT(ndY->IsInternal());
    if (!likelihood->getNoData())
        {
        likelihood->useAsLikelihoodRoot(ndY);
        likelihood->invalidateAwayFromNode(*ndY);
        likelihood->invalidateBothEnds(ndY);

        //@POL temp
        //likelihood->invalidateBothEnds(ndZ);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reverses move made in proposeNewState. Assumes ndX, ndY, and ndZ are non-NULL, which will be true if proposeNewState
|	was just called.
*/
void LargetSimonMove::revert()
	{
    ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

	MCMCUpdater::revert();
	if (star_tree_proposal)
		{
		orig_node->SetEdgeLen(orig_edge_len);
		TreeNode * nd = orig_node->IsTip() ? orig_node->GetParent() : orig_node;
		PHYCAS_ASSERT(nd->IsInternal());
		likelihood->useAsLikelihoodRoot(nd);
		likelihood->restoreFromCacheAwayFromNode(*orig_node);
		likelihood->restoreFromCacheParentalOnly(orig_node);

		orig_node->UnselectNode();
        jpm->allEdgeLensModified(tree);
        //jpm->externalEdgeLensModified("external_edgelen", tree);
		}
	else
		{
		PHYCAS_ASSERT(ndX != NULL);
		PHYCAS_ASSERT(ndY != NULL);
		PHYCAS_ASSERT(ndZ != NULL);
		PHYCAS_ASSERT(topol_changed ? (swap1 != NULL && swap2 != NULL) : (swap1 == NULL && swap2 == NULL));

		if (topol_changed)
			{
			if (swap2 == ndZ)
				{
				// If swap2 equals ndZ, then swap2 was a child of ndBase and we were able to use
				// the standard NNISwap function to swap the two nodes
				//
				tree_manipulator.NNISwap(swap1, swap2);
				}
			else
				{
				// If swap2 is ndZ's parent, then swap2 is ndBase (i.e. it is the "child" node below the
				// lower of the two adjacent internal nodes involved in the swap) and we had to use the
				// NNISwapSpecial function to perform the rearrangment
				//
				tree_manipulator.NNISwapSpecial(swap1);
                likelihood->swapInternalDataAndEdgeLen(ndY, ndZ);
				}
			}
		ndX->SetEdgeLen(x);
		ndY->SetEdgeLen(y);
		ndZ->SetEdgeLen(z);

		PHYCAS_ASSERT(ndY->IsInternal());
        if (!likelihood->getNoData())
            {
            likelihood->useAsLikelihoodRoot(ndY);
            likelihood->restoreFromCacheAwayFromNode(*ndY);
            likelihood->restoreFromCacheParentalOnly(ndY);
            }

		ndX->UnselectNode();
		ndY->UnselectNode();
		ndZ->UnselectNode();

        jpm->allEdgeLensModified(tree);
        if (topol_changed)
            jpm->topologyModified("tree_topology", tree);
		}

    curr_ln_prior = jpm->getLogJointPrior();

	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void LargetSimonMove::accept()
	{
	MCMCUpdater::accept();
	if (star_tree_proposal)
		{
		TreeNode * nd = orig_node->IsTip() ? orig_node->GetParent() : orig_node;
		PHYCAS_ASSERT(nd->IsInternal());
        if (!likelihood->getNoData())
            {
            likelihood->useAsLikelihoodRoot(nd);
            likelihood->discardCacheAwayFromNode(*orig_node);
            likelihood->discardCacheBothEnds(orig_node);
            }

		orig_node->UnselectNode();
		}
	else
		{
		PHYCAS_ASSERT(ndY->IsInternal());
        if (!likelihood->getNoData())
            {
            likelihood->useAsLikelihoodRoot(ndY);
            likelihood->discardCacheAwayFromNode(*ndY);
            likelihood->discardCacheBothEnds(ndY);
            }

		ndX->UnselectNode();
		ndY->UnselectNode();
		ndZ->UnselectNode();
		}

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Forgets information saved to enable reverting a proposed move.
*/
void LargetSimonMove::reset()
	{
	topol_changed   = false;
	swap1		    = NULL;
	swap2		    = NULL;
	ndX			    = NULL;
	ndY			    = NULL;
	ndZ			    = NULL;
	//ndBase		    = NULL;
	x		        = 0.0;
	y		        = 0.0;
	z		        = 0.0;
    which_case      = 0;

    expand_contract_factor = 1.0;

	// three_edgelens should have 3 elements, used for computing the edge length prior for default proposal in update
	three_edgelens.resize(3, 0.0);

	// these are related to the star tree exception
	orig_edge_len	= 0.0;
	orig_node		= NULL;

	// one_edgelen should have 1 element, used for computing the edge length prior for star tree proposal in update
	one_edgelen.resize(1, 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides a way of setting the value for the private data member `lambda'.
*/
void LargetSimonMove::setLambda(double x)
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the private data member `lambda'.
*/
double LargetSimonMove::getLambda() const
	{
	return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is (`mstar'/`m')^3, where `mstar' is
|	the new length of the modified three-edge segment and `m' is the length before the move is proposed. If the tree is
|	the star tree (only one internal node), then the exponent is 1 rather than 3 because only one edge is involved.
*/
double LargetSimonMove::getLnHastingsRatio() const
	{
	if (star_tree_proposal)
		return std::log(orig_node->GetEdgeLen()) - std::log(orig_edge_len);
	else
		return 3.0*(std::log(mstar) - std::log(m));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double LargetSimonMove::getLnJacobian() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of private variable `which_case', which is the index of whichever of the eight possible
|	cases was used last.
*/
unsigned LargetSimonMove::getWhichCase() const
	{
	return which_case;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the private data member `topol_changed'.
*/
bool LargetSimonMove::topologyChanged() const
	{
	return topol_changed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'min_lambda', which is the tuning parameter used for exploring the posterior
|   distribution in this move.
*/
void LargetSimonMove::setPosteriorTuningParam(
  double x) /* is the new value for `min_lambda' */
	{
	min_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member `max_lambda', which is the tuning parameter used for exploring the prior
|   distribution in this move.
*/
void LargetSimonMove::setPriorTuningParam(
  double x) /* is the new value for `max_lambda' */
	{
	max_lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'lambda', which is the tuning parameter for this move, based on a boldness value
|	that ranges from 0 (least bold) to 100 (most bold). Simple linear interpolation is used (i.e. a boldness of 50
|   results in `lambda' halfway between `min_lambda' and `max_lambda').
*/
void LargetSimonMove::setBoldness(
  double x) /* is the new boldness value */
	{
	boldness = x;
	if (boldness < 0.0)
		boldness = 0.0;
	else if (boldness > 100.0)
		boldness = 100.0;

    // compute lambda from boldness value
	lambda = min_lambda + (max_lambda - min_lambda)*boldness/100.0;
	}

}	// namespace phycas

