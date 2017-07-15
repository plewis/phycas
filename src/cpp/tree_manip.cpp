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

#include <numeric>
#include "tree_manip.hpp"
#include "basic_lot.hpp"
#include "probability_distribution.hpp"
#include "tree_length_distribution.hpp"
#include "basic_tree.hpp"
#include "tree_iterators.hpp"

//temporary
#include <fstream>

using namespace phycas;

#include "basic_lot.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|   Selects an internal node at random from a discrete uniform distribution with the constraint that the returned node
|   is not equal to the subroot (the sole child of the tip node serving as the root).
*/
TreeNode * randomInternalAboveSubroot(Tree & tree, Lot & rng)
    {
	// Avoiding the "subroot" node (only child of the tip serving as the root), so the number of
	// acceptable nodes is one fewer than the number of internal nodes
	unsigned numAcceptableNodes = tree.GetNInternals() - 1;

	unsigned ypos = rng.SampleUInt(numAcceptableNodes);
	unsigned i = 0;
    TreeNode * nd = tree.GetFirstPreorder();
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
|	The default constructor does not set the data member `tree'. You should set `tree' using the SetTree() member
|	function before using any of the other member functions, most of which assume that `tree' points to a Tree object.
*/
TreeManip::TreeManip()
  	{
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor requires a single argument, which is a shared pointer to the tree to be manipulated.
*/
TreeManip::TreeManip(
  TreeShPtr t)	/**< Is the tree to be manipulated */
  	{
	PHYCAS_ASSERT(t);
	tree = t;
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	The destructor need not do anything, but it can be useful to uncomment the output statement to see when the
|	destructor is called.
*/
TreeManip::~TreeManip()
	{
	//std::cerr << "\n\n>>>>> TreeManip is dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' data member to the supplied shared pointer to Tree.
*/
void TreeManip::setTree(TreeShPtr t)
  	{
	tree = t;
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds new leaf node to a randomly-selected internal vertex. Does not strip conditional likelihoods, so care should be
|   taken with this function if the likelihood will be calculated on the resulting tree.
*/
void TreeManip::addLeafToRandomNode(unsigned leaf_node_number, LotShPtr rng)
    {
    // POL_BOOKMARK 11-July-2017

    //std::cerr << "\n~~> begin TreeManip::addLeafToRandomNode" << std::endl;

    //tree->debugMode(true);
	//tree->DebugCheckTree(false, false, 0);

    // Choose an internal node at random
	tree->InvalidateNodeCounts();
    unsigned num_nodes = tree->GetNNodes();
    PHYCAS_ASSERT(num_nodes > 0);
    unsigned num_internals = tree->GetNInternals();
    //unsigned num_tips = num_nodes - num_internals;
    //std::cerr << boost::str(boost::format("~~> num_nodes = %d, num_internals = %d, num_tips = %d, leaf_node_number = %d") % num_nodes % num_internals % num_tips % leaf_node_number) << std::endl;

    PHYCAS_ASSERT(num_internals > 0);
    TreeNode * nd = tree->GetFirstPreorder();
    unsigned i = rng->SampleUInt(num_internals);    // returns value from [0, 1, ..., num_internals-1]
    //std::cerr << boost::str(boost::format("~~> adding to node = %d (num_internals = %d)") % i % num_internals) << std::endl;
    for (; nd != NULL; nd = nd->GetNextPreorder())
        {
        if (!nd->IsTip())
            {
            if (i == 0)
                break;
            else
                --i;
            }
        }
    //std::cerr << boost::str(boost::format("~~> nd->nodeNum = %d") % nd->nodeNum) << std::endl;

    // Add a new leaf node to nd to create a polytomy (or a bigger polytomy if nd is already a polytomy).
    TreeNode * new_node = tree->GetNewNode();
    new_node->nodeNum = leaf_node_number;

	InsertSubtree(new_node, nd, TreeManip::kOnRight);
	tree->InvalidateNodeCounts();
    tree->preorderDirty = true;
	tree->RefreshPreorder();

	//tree->DebugCheckTree(false, false, 0);
    //tree->debugMode(false);

    //std::cerr << "~~> end TreeManip::addLeafToRandomNode" << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Adds new leaf node to a randomly-selected edge. Does not strip conditional likelihoods, so care should be taken with
|   this function if the likelihood will be calculated on the resulting tree.
*/
void TreeManip::addLeafToRandomEdge(unsigned leaf_node_number, LotShPtr rng)
    {
    // POL_BOOKMARK 11-July-2017

    //std::cerr << "\n~~> begin TreeManip::addLeafToRandomEdge" << std::endl;

    tree->debugMode(true);
	tree->DebugCheckTree(false, false, 0);

    // Choose a node at random
	tree->InvalidateNodeCounts();
    unsigned num_nodes = tree->GetNNodes();
    PHYCAS_ASSERT(num_nodes > 0);
    unsigned num_internals = tree->GetNInternals();
    //unsigned num_tips = num_nodes - num_internals;
    //std::cerr << boost::str(boost::format("~~> num_nodes = %d, num_internals = %d, num_tips = %d, leaf_node_number = %d") % num_nodes % num_internals % num_tips % leaf_node_number) << std::endl;

    TreeNode * nd = tree->GetFirstPreorder();
    unsigned i = 1 + rng->SampleUInt(num_nodes-1);

    //std::cerr << "~~> adding new node to edge of " << i << (i == 1 ? "st" : (i == 2 ? "nd" : (i == 3 ? "rd" : "th"))) << " node in preorder sequence" << std::endl;

    for (; nd != NULL; nd = nd->GetNextPreorder())
        {
        if (i == 0)
            break;
        else
            --i;
        }

    PHYCAS_ASSERT(nd->GetParent());

    //std::cerr << "~~> nd->nodeNum = " << nd->nodeNum << std::endl;

    // Create new internal node
    TreeNode * new_anc  = tree->GetNewNode();
    new_anc->nodeNum   = num_nodes;

    // Create new leaf node
    TreeNode * new_leaf = tree->GetNewNode();
    new_leaf->nodeNum   = leaf_node_number;

    // Add new leaf to new internal
    new_anc->lChild     = new_leaf;
    new_leaf->par       = new_anc;

    // Detach nd from tree
	TreeNode * nd_parent = nd->GetParent();
    DetachSubtree(nd);

    // Make nd the right sibling of the new leaf
    new_leaf->rSib = nd;
    nd->par        = new_anc;

    // Finally, add new internal to nd's former parent
	InsertSubtree(new_anc, nd_parent, TreeManip::kOnRight);

    tree->preorderDirty = true;
	tree->RefreshPreorder();
	tree->InvalidateNodeCounts();

    num_nodes = tree->GetNNodes();
    num_internals = tree->GetNInternals();


	//tree->DebugCheckTree(false, false, 0);
    //tree->debugMode(false);

    //std::cerr << "~~> end TreeManip::addLeafToRandomEdge\n" << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Labels leaves of tree randomly.
*/
void TreeManip::assignLeafNumbersRandomly(LotShPtr rng)
    {
    // POL_BOOKMARK 11-July-2017
    //std::cerr << "\n~~> begin TreeManip::assignLeafNumbersRandomly\n" << std::endl;

	tree->InvalidateNodeCounts();
    unsigned num_nodes = tree->GetNNodes();
    PHYCAS_ASSERT(num_nodes > 0);
    unsigned num_internals = tree->GetNInternals();
    unsigned num_tips = num_nodes - num_internals;
    //std::cerr << boost::str(boost::format("~~> num_nodes = %d, num_internals = %d, num_tips = %d") % num_nodes % num_internals % num_tips) << std::endl;

    // Create pool of available tip node numbers
    std::set<unsigned> pool;
    for (unsigned i = 0; i < num_tips; ++i)
        pool.insert(i);

    TreeNode * nd = tree->GetFirstPreorder();
    for (; nd != NULL; nd = nd->GetNextPreorder())
        {
        if (nd->IsTip())
            {
            unsigned i = rng->SampleUInt((unsigned)pool.size());
            std::set<unsigned>::iterator it = pool.begin();
            for (unsigned j = 0; j < i; ++j)
                it++;
            nd->nodeNum = *it;
            pool.erase(it);
            }
        }
    //std::cerr << "~~> end TreeManip::assignLeafNumbersRandomly\n" << std::endl;
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	Deletes a randomly-selected edge from the tree. Does not strip conditional likelihoods before storing internal
|   nodes, so care should be taken with this function if the likelihood will be calculated on the resulting tree.
*/
void TreeManip::deleteRandomInternalEdge(LotShPtr rng)
    {
    // Choose an internal node at random (but not the only child of the root node) and delete
    // its edge to create a polytomy (or a bigger polytomy if there is already a polytomy)
    unsigned num_internals = tree->GetNInternals();
    PHYCAS_ASSERT(num_internals > 0);
    unsigned i = rng->SampleUInt(num_internals - 1);

    TreeNode * nd = tree->GetFirstPreorder();
    for (; nd != NULL; nd = nd->GetNextPreorder())
        {
        bool is_tip = nd->IsTip();
        bool is_subroot = (is_tip ? false : nd->GetParent()->IsTipRoot());
        if (!is_tip && !is_subroot)
            {
            if (i == 0)
                break;
            else
                --i;
            }
        }

	// This operation should not leave the root node (which is a tip) with more than one child,
	// so check to make sure that nd is not the root nor a child of root
	TreeNode * parent = nd->GetParent();
	PHYCAS_ASSERT(parent != NULL);
	PHYCAS_ASSERT(!parent->IsTipRoot());

	// Make all of nd's children children of parent instead
	while (nd->GetLeftChild() != NULL)
		{
		TreeNode * rchild = nd->GetLeftChild();
		DetachSubtree(rchild);
		InsertSubtree(rchild, parent, TreeManip::kOnRight);
		}

	DeleteLeaf(nd, true);
	tree->InvalidateNodeCounts();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns all edge lengths in the tree using independent draws from the ProbabilityDistribution object pointed to by
|	`d'.
*/
void TreeManip::setRandomEdgeLens(ProbDistShPtr d)
	{
	//@POL would like to do something like that shown below, but doesn't compile with d wrapped in var, and not using
	// var wrapper results in Sample() only getting called once (all edge lengths set to the same value)
	//std::for_each(tree->begin(), tree->end(), boost::lambda::bind(&TreeNode::SetEdgeLen, boost::lambda::_1, boost::lambda::var(d)->Sample()));
	for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
		{
		nd->SetEdgeLen(d->Sample());
		}
    tree->hasEdgeLens = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns all internal edge lengths in the tree using independent draws from the ProbabilityDistribution object
|   pointed to by `internal_dist', and all external edge lengths using draws from the ProbabilityDistribution
|   object pointed to by `external_dist'.
*/
void TreeManip::setRandomInternalExternalEdgeLens(
  ProbDistShPtr internal_dist,  /**< is the distribution to use for choosing internal edge lengths */
  ProbDistShPtr external_dist)  /**< is the distribution to use for choosing external edge lengths */
	{
	for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
		{
		if (nd->IsTip())
		    nd->SetEdgeLen(external_dist->Sample());
		else
		    nd->SetEdgeLen(internal_dist->Sample());
		}
    tree->hasEdgeLens = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Begins with left child of parent of `this' and calls GetRightSib() until the left sibling of `this' is located.
*/
//TreeNode * TreeNode::FindLeftSib()
//	{
//	TreeNode * nd = this->GetParent();
//
//	// If this has no parent, then there's no way it can have a left sibling
//	if (nd == NULL)
//		return NULL;
//
//	nd = nd->GetLeftChild();
//
//	// Parent of this should have at least one child
//	PHYCAS_ASSERT(nd != NULL);
//
//	// If left child of this's parent equals this, then this is an only child and has no siblings
//	if (nd == this)
//		return NULL;
//
//	TreeNode * leftsib = NULL;
//	while (nd != this)
//		{
//		if (nd == NULL)
//			{
//			throw XPhylogeny("pointer inconsistency in FindLeftSib");
//			}
//		leftsib = nd;
//		nd = nd->GetRightSib();
//		}
//	return leftsib;
//	}

/*----------------------------------------------------------------------------------------------------------------------
|	Begins with left child of `this'. If left child is NULL, returns NULL, otherwise, returns the rightmost sibling of
|	the left child of `this'
*/
//TreeNode * TreeNode::FindRightmostChild()	/**< is the parent node whose children will be searched */
//	{
//	TreeNode * curr = GetLeftChild();
//	TreeNode * rightmost = NULL;
//	while (curr != NULL)
//		{
//		rightmost = curr;
//		curr = curr->GetRightSib();
//		}
//	return rightmost;
//	}

/*----------------------------------------------------------------------------------------------------------------------
|	Zig-zags up the right side of the clade starting with the node `start' until reaching a tip. This tip node is the
|	last node in the clade according to the preorder sequence. Zig-zagging is necessary because there is no pointer to
|	the right-most child of a node, so we must use FindRightmostChild to move up the clade. Assumes `start' is non-NULL.
*/
TreeNode * TreeManip::FindLastPreorderInClade(
  TreeNode * start)	/**< is the deepest node belonging to the clade in question */
	{
	PHYCAS_ASSERT(start != NULL);
	TreeNode * curr = start;
	TreeNode * rChild = FindRightmostChild(curr);
	while (rChild != NULL)
		{
		curr = rChild;
		rChild = FindRightmostChild(curr);
		}
	return curr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts `subtree' into the edge owned by `edge_nd', adding an internal node.
*/
void TreeManip::InsertSubtreeIntoEdge(
  TreeNode * subtree, /**< */
  TreeNode * edge_nd) /**< */
    {
    // get new_node from internal node storage
    TreeNode * new_node = tree->GetNewNode();

    // add new_node to edge_nd->par
	//if (Tree::gDebugOutput)
    //	{
	//	std::cerr << "*** InsertSubtreeIntoEdge:  add new_node to edge_nd->par" << std::endl; //temporary
	//	std::cerr << "new_node" << new_node->oneLineDebugReport() << std::endl;
	//	std::cerr << "edge_nd" << edge_nd->oneLineDebugReport() << std::endl;
	//	}
	PHYCAS_ASSERT(edge_nd->par);
	InsertSubtree(new_node, edge_nd->par, TreeManip::kOnRight);
	//tree->DebugCheckTree(false, false, 2);    // POL_BOOKMARK 14-July-2017

    // add subtree to new_node
    InsertSubtree(subtree, new_node, TreeManip::kOnRight);
	//if (Tree::gDebugOutput)
    //	{
	//	std::cerr << "*** InsertSubtreeIntoEdge:  add subtree to new_node" << std::endl; //temporary
	//	tree->DebugCheckTree(true, true, 2);
	//	}
    // make edge_nd a child of new_node
    SibToChild(new_node, edge_nd,  TreeManip::kOnRight);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts node 's' as child of 'u' keeping the subtree rooted at `s' intact. Assumes both `s' and `u' are non-NULL. If
|	`targetSib' is specified, `s' will become an immediate sibling of `targetSib' (right or left sibling determined by
|	`m'). If `targetSib' is not specified, `s' will be added as the rightmost or leftmost child of `u' (again, depending
|	on the value of `m').
*/
void TreeManip::InsertSubtree(
  TreeNode * s,			/**< is the node to be added as a child of u */
  TreeNode * u,			/**< is the new parent of s */
  InsertMode m,			/**< if kOnRight, s will be added to the right of targetSib or, if targetSib is not specified, s will become rightmost child of u (vice versa if kOnLeft specified) */
  TreeNode * targetSib)	/**< if non-NULL, s will become an immediate sibling of this node (right or left depending on value of `m') */
	{
	PHYCAS_ASSERT(u != NULL);
	PHYCAS_ASSERT(s != NULL);
	PHYCAS_ASSERT(targetSib == NULL || (targetSib->par == u));

	TreeNode * slast				= FindLastPreorderInClade(s);
	TreeNode * u_lChild				= u->lChild;
	TreeNode * u_rChild				= FindRightmostChild(u);
	TreeNode * ulast				= FindLastPreorderInClade(u);
	TreeNode * ulast_nextPreorder	= ulast->nextPreorder;
	TreeNode * u_nextPreorder		= u->nextPreorder;

	// Identify possible reductions in complexity
	//
	if (targetSib != NULL)
		{
		if ((m == TreeManip::kOnLeft) && (targetSib == u_lChild))
			targetSib = NULL;
		if ((m == TreeManip::kOnRight) && (targetSib == u_rChild))
			targetSib = NULL;
		if ((m == TreeManip::kOnLeft) && (targetSib != NULL))
			{
			targetSib = FindLeftSib(targetSib);
			m = TreeManip::kOnRight;
			}
		}

	// Make s a child of u
	//
	if (targetSib == NULL)
		{
		if (m == TreeManip::kOnRight)
			{
            s->ptr = (void *)u;    //temporary!
			s->rSib				= NULL;
			s->par				= u;
			s->prevPreorder 	= ulast;
			slast->nextPreorder = ulast_nextPreorder;
			ulast->nextPreorder	= s;

			if (ulast_nextPreorder == NULL)
				tree->lastPreorder = slast;
			else
				ulast_nextPreorder->prevPreorder = slast;

			if (u_lChild == NULL)
				u->lChild		= s;
			else
				u_rChild->rSib	= s;
			}
		else
			{
			s->rSib					= u_lChild;
			s->par					= u;
			s->prevPreorder 		= u;
			slast->nextPreorder		= u_nextPreorder;

			if (u_lChild != NULL)
				u_lChild->prevPreorder	= slast;
			u->lChild				= s;
			u->nextPreorder			= s;
			}
		}
	else
		{
		// Save pointers to relevant nodes
		//
		TreeNode * targetSib_rSib				= targetSib->rSib;
		TreeNode * targetSiblast				= FindLastPreorderInClade(targetSib);
		TreeNode * targetSiblast_nextPreorder	= targetSiblast->nextPreorder;

		// Can assume that targetSib is not rightmost child of u and also that m == TreeManip::kOnRight
		//
		s->rSib						= targetSib_rSib;
		s->par						= u;
		s->prevPreorder 			= targetSiblast;
		slast->nextPreorder 		= targetSiblast_nextPreorder;
		targetSiblast->nextPreorder	= s;
		targetSib->rSib				= s;
		}

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a star tree having `ntips' tips. Each edge has a length drawn from the supplied ProbabilityDistribution. It
|	is assumed that the random number generator has already been set for `edge_len_dist'.
*/
void TreeManip::starTree(
  unsigned		ntips,			/**< Number of tip nodes (tree will have one more tip than this if rooted) */
  ProbDistShPtr	edge_len_dist)	/**< Probability distribution from which to draw the edge lengths */
	{
	PHYCAS_ASSERT(ntips > 0);
	//PHYCAS_ASSERT(edge_len_dist);
	PHYCAS_ASSERT(tree);

	tree->Clear();
    tree->DeleteAllStoredNodes();
    //tree.reset(new Tree());
	tree->hasEdgeLens = true;

	// Create the root node
	//
	TreeNode * rootNd = tree->GetNewNode();
	rootNd->SetEdgeLen(0.0);
    unsigned first_tip_number;
    if (tree->IsRooted())
        {
        first_tip_number = 0;
    	rootNd->SetNodeNum(ntips);
        rootNd->SetObservable(false);
        }
    else
        {
        first_tip_number = 1;
    	rootNd->SetNodeNum(0);
        rootNd->SetObservable(true);
        }
	tree->firstPreorder = rootNd;

	// Create the hub node
	//
    const double hubLen = (edge_len_dist ? edge_len_dist->Sample() : 1.0);
	TreeNode * hub = tree->GetNewNode();
	hub->SetEdgeLen(hubLen);
	hub->SetNodeNum(ntips);
	InsertSubtree(hub, rootNd, TreeManip::kOnLeft);

	for (unsigned i = first_tip_number; i < ntips; ++i)
		{
		const double ndEdgeLen = (edge_len_dist ? edge_len_dist->Sample() : 1.0);
		TreeNode * nd = tree->GetNewNode();
		nd->SetEdgeLen(ndEdgeLen);
		nd->SetNodeNum(i);
        nd->SetObservable(true);
		InsertSubtree(nd, hub, TreeManip::kOnLeft);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If the `yule' parameter is true, and the supplied probability distribution is ExponentialDist(lambda), this function
|	builds a Yule tree with speciation rate `lambda'. If the `yule' parameter is false, then the supplied probability
|	distribution is used to generate independent edge lengths and the resulting tree is not expected to be ultrametric.
|	If you want the same random number generator to be used for both edge lengths and branching order, call the setLot()
|	function of the probability distribution object, passing the same Lot object used in the call to
|	TreeManip::randomTree. The resulting tree is unrooted, so the supplied `ntips' includes the tip serving as the root
|	node. A Yule tree is more appropriately interpreted as a rooted tree, but thus far Phycas only deals with unrooted
|	trees so we are using this function as simply a way to generate a random topology.
|
|	Yule trees: Let Sk be the sojourn time between the kth and (k+1)st birth. Let Wk be	the time of the kth birth. Let
|	X(t) be the number of lineages existing at time t. Each lineage is independent and splits with probability `lambda'
|	dt in each infinitesimal time interval.
|>
|	------------------------------------------ W4 = time of 4th birth = S0 + S1 + S2 + S3
|	 \             \        /     |
|	  \             \      /      |
|	   \             \    /       S3    X(t) = 3
|	    \             \  /        |
|	     \             \/         |
|	------------------------------------------ W3 = time of 3rd birth = S0 + S1 + S2
|	      \            /          |
|	       \          /           |
|	        \        /            |
|	         \      /             S2    X(t) = 2
|	          \    /              |
|	           \  /               |
|	            \/                |
|	------------------------------------------ W2 = time of 2nd birth = S0 + S1
|	            /                 |
|	           /                  |
|	          /                   S1    X(t) = 1
|	         /                    |
|	        /                     |
|	------------------------------------------ W1 = time of 1st birth = S0 = 0.0
|>
|	The cumulative distribution function of the sojourn time S1 given speciation rate `lambda' is
|>
|	Pr(S1 < t) = 1 - Pr{no speciation in lineage | lambda, t} = 1 - exp{-lambda t}
|>
|	Thus, S1 is exponentially distributed with rate parameter `lambda'. The cumulative distribution function of the
|	sojourn time S2 is
|>
|	Pr(S2 < t) = 1 - Pr{no speciation in either lineage | lambda, t}
|	           = 1 - Pr{no speciation in lineage 1} Pr{no speciation in lineage 2}
|	           = 1 - exp{-lambda t} exp{-lambda t}
|	           = 1 - exp{-2 lambda t}
|>
|	Thus, S2 is exponentially distributed with hazard parameter 2*`lambda'. In general, for j > 0, sojourn time Sj is
|	exponentially distributed with hazard parameter j*`lambda', where j is the number of extant lineages.
|
|	While finding the MLE of `lambda' is not relevant here, it might be useful to calling functions to obtain the
|	MLE of `lambda' in order to supply that value to this function (to create a Yule tree that is similar to a certain
|	known tree, for example). If one is given T, the total tree length, and the number of non-root tip nodes n, the MLE
|	of `lambda' is simply n/T. Here is the derivation for a tree with 4 non-root tips as an example:
|>
|	Pr(S1) = lambda*exp(-lambda*S1)
|	Pr(S2) = 2*lambda*exp(-2*lambda*S2)
|	Pr(S3) = 3*lambda*exp(-3*lambda*S3)
|	Pr(S4) = 4*lambda*exp(-4*lambda*S4)
|>
|	The likelihood of `lambda' is thus proportional to the product of these four expressions:
|>
|	L = (4!)*(lambda^4)*exp(-lambda*T)
|>
|	The log-likelihood is thus
|>
|	lnL = n ln(lambda) - lambda*T
|>
|	The derivative of lnL with respect to lambda is n/lambda - T, the maximum of which is at n/T.
*/
void TreeManip::randomTree(
  unsigned ntips,				/**< is the number of tip nodes in the final tree (includes the tip node serving as the root) */
  LotShPtr	rng,				/**< is the random number generator used to determine the branching order (but not edge lengths) */
  ProbDistShPtr	edge_len_dist,	/**< is the probability distribution used to generate new edge lengths (should be exponential(lambda) if `yule' is true, where lambda is the instantaneous speciation rate) */
  bool yule)					/**< if true, a Yule tree will be generated; otherwise, edge lengths will be independent draws from `edge_len_dist' */
	{
	//@POL should assert that edge_len_dist is of type ExponentialDist if yule is true?
	PHYCAS_ASSERT(edge_len_dist);
	PHYCAS_ASSERT(tree);

    // Empty the tree of everything and simply return if ntips is 0
	tree->Clear();
    tree->DeleteAllStoredNodes();
    //tree.reset(new Tree());
    if (ntips == 0)
        return;

	tree->hasEdgeLens = true;

	unsigned i, j, k;
	double new_edgelen;
	unsigned nextTipNodeNum = 2; // 0 is reserved for root node, and 1 is given to the first tip created
	unsigned nextInternalNodeNum = ntips;

	// This vector keeps all nodes representing growing tips handy
	//
	std::vector<TreeNode *> shootTips;

	// Create the root node
	//
	TreeNode * rootNd = tree->GetNewNode();
	rootNd->SetEdgeLen(0.0);
	rootNd->SetNodeNum(0);
	tree->firstPreorder = rootNd;

	// Create the first leaf node (note: root is a tip but not a leaf, shoot tip is synonymous with leaf)
	//
	TreeNode * nd = tree->GetNewNode();
	nd->SetEdgeLen(yule ? 0.0 : edge_len_dist->Sample());
	nd->SetNodeNum(1);
	shootTips.push_back(nd);

	// Connect the root node to the leaf node
	//
	rootNd->lChild = nd;
	rootNd->rSib = NULL;
	rootNd->par = NULL;

	nd->lChild = NULL;
	nd->rSib = NULL;
	nd->par = rootNd;

	for (i = 2; i < ntips; ++i)
		{
		unsigned num_leaves = i - 1;

		// Choose one of the current leaves to speciate
		//
		j = rng->SampleUInt(num_leaves);
		nd = shootTips[j];

		// Create an internal node
		//
		TreeNode * parent = tree->GetNewNode();
		parent->SetEdgeLen(yule ? 0.0 : edge_len_dist->Sample()); // was DBL_MAX, can that be right?
		parent->SetNodeNum(nextInternalNodeNum++);

		// parent becomes nd, nd can then remain a leaf
		//
		parent->lChild = nd;
		parent->rSib = nd->rSib;
		parent->par = nd->par;

		// Check for rSib pointer directed at nd, redirect to parent
		//
		TreeNode * lSib = FindLeftSib(nd);
		if (lSib != NULL)
			lSib->rSib = parent;

		// Check for lChild pointer directed at nd, redirect to parent
		//
		PHYCAS_ASSERT(nd->par != NULL);
		if (nd->par->lChild == nd)
			nd->par->lChild = parent;

		nd->lChild = NULL;
		nd->rSib = NULL;
		nd->par = parent;

		if (yule)
			{
			parent->SetEdgeLen(nd->GetEdgeLen());
			nd->SetEdgeLen(0.0);
			}

		// Create a new tip representing the sister taxon
		//
		TreeNode * sister = tree->GetNewNode();
		sister->SetEdgeLen(yule ? 0.0 : edge_len_dist->Sample());
		sister->SetNodeNum(nextTipNodeNum++);
		nd->rSib = sister;

		sister->lChild = NULL;
		sister->rSib = NULL;
		sister->par = parent;

		// Choose a new edge length if building a Yule tree
		//
		if (yule)
			{
			new_edgelen = edge_len_dist->Sample()/num_leaves;

			// Add the new edge length to all existing tips
			//
			for (k = 0; k < num_leaves; ++k)
				{
				if (k == j)
					parent->SetEdgeLen(new_edgelen + parent->GetEdgeLen());
				else
					shootTips[k]->SetEdgeLen(new_edgelen + shootTips[k]->GetEdgeLen());
				}
			}

		// Add sister to the vector of shoot tips
		//
		shootTips.push_back(sister);
		}

	// Choose a final edge length if building a Yule tree
	//
	if (yule)
		{
		PHYCAS_ASSERT(i == ntips);
		new_edgelen = edge_len_dist->Sample()/ntips;

		// Add the new edge length to all existing tips
		//
		for (k = 0; k < ntips - 1; ++k)
			{
			shootTips[k]->SetEdgeLen(new_edgelen + shootTips[k]->GetEdgeLen());
			}
		}

	tree->RefreshPreorder();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs an equiprobable tree topology and sets all edge lengths using draws from the supplied `edge_len_dist'.
*/
void TreeManip::equiprobTree(
  unsigned ntips,				    /**< is the number of tip nodes in the final tree (includes the tip node serving as the root) */
  LotShPtr	rng,				    /**< is the random number generator used to determine the branching order (but not edge lengths) */
  ProbDistShPtr	internal_edge_dist,	/**< is the probability distribution used to generate new internal edge lengths  */
  ProbDistShPtr	external_edge_dist)	/**< is the probability distribution used to generate new external edge lengths  */
	{
	PHYCAS_ASSERT(tree);

    // Empty the tree of everything and simply return if ntips is 0
	tree->Clear();
    tree->DeleteAllStoredNodes();
    //tree.reset(new Tree());
    if (ntips == 0)
        return;

    unsigned tip_number = 0;
    unsigned internal_number = ntips;

    // Create the root node (first tip)
	TreeNode * rootNd = tree->GetNewNode();
    rootNd->SetNodeNum(tip_number++);
    rootNd->SetObservable();
	tree->firstPreorder = rootNd;

	// Create the subroot node (which is second tip)
	TreeNode * nd = tree->GetNewNode();
    nd->SetNodeNum(tip_number++);
    nd->SetObservable();
    nd->SetEdgeLen(1.0);

    // This vector makes it easy to randomly choose next node to bisect
    //typedef std::vector<TreeNode *> NodeVector;
	//NodeVector receptive;
    std::vector<TreeNode *> receptive;
    receptive.reserve(2*ntips - 3);
    receptive.push_back(nd);

	// Connect the root node to the subroot node
	rootNd->lChild = nd;
	rootNd->rSib = NULL;
	rootNd->par = NULL;
	nd->lChild = NULL;
	nd->rSib = NULL;
	nd->par = rootNd;

	for (unsigned i = 2; i < ntips; ++i)
		{
        unsigned num_receptive_nodes = 2*i - 3;
        PHYCAS_ASSERT(num_receptive_nodes == (unsigned)receptive.size());

		// Choose one of the receptive nodes (i.e. any node except rootNd)
		unsigned j = rng->SampleUInt(num_receptive_nodes);
		nd = receptive[j];

		// Create an internal node on nd's edge
		TreeNode * parent = tree->GetNewNode();
        parent->SetNodeNum(internal_number++);
        parent->SetEdgeLen(1.0);
        receptive.push_back(parent);
		parent->lChild = nd;
		parent->rSib = nd->rSib;
		parent->par = nd->par;
		TreeNode * lSib = FindLeftSib(nd);
		if (lSib != NULL)
			lSib->rSib = parent;
		PHYCAS_ASSERT(nd->par != NULL);
		if (nd->par->lChild == nd)
			nd->par->lChild = parent;
		nd->rSib = NULL;
		nd->par = parent;

        // Create a new tip representing nd's sister
		TreeNode * sister = tree->GetNewNode();
        sister->SetNodeNum(tip_number++);
        sister->SetObservable();
        sister->SetEdgeLen(1.0);
        receptive.push_back(sister);
		nd->rSib = sister;
        sister->par = parent;
		sister->lChild = NULL;
		sister->rSib = NULL;
		}

	tree->RefreshPreorder();

    // Set edge lengths
	bool set_all_edge_lengths_to_one = false;
    if (internal_edge_dist && !external_edge_dist)
        external_edge_dist = internal_edge_dist;
    else if (external_edge_dist && !internal_edge_dist)
        internal_edge_dist = external_edge_dist;
	else if (!internal_edge_dist && !external_edge_dist)
		set_all_edge_lengths_to_one = true;

    //double tree_length = 0.0;
    for (std::vector<TreeNode *>::iterator it = receptive.begin(); it != receptive.end(); ++it)
        {
        TreeNode * nd = (*it);
        double new_edgelen = 1.0;
        if (nd->IsInternal() && !set_all_edge_lengths_to_one)
            {
			new_edgelen = internal_edge_dist->Sample();
            }
        else if (!set_all_edge_lengths_to_one)
            {
			new_edgelen = external_edge_dist->Sample();
            }
        //tree_length += new_edgelen;
        nd->SetEdgeLen(new_edgelen);
        }
	tree->hasEdgeLens = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets edge lengths of `tree' using supplied TreeLengthDistribution.
*/
void TreeManip::setRandomEdgeLensFromTreeLengthDist(
  TreeLengthDistributionShPtr tree_length_dist)   /**< is the tree length distribution to be used to generate edge lengths */
    {
	PHYCAS_ASSERT(tree);
	PHYCAS_ASSERT(tree_length_dist);
	PHYCAS_ASSERT(!tree->IsRooted());

    unsigned num_external_edges = tree->GetNTips();
    unsigned num_internal_edges = tree->GetNInternals() - 1;    // one internal node represents a tip

    // Get edge lengths from tree length distribution: external edges are first, followed by internals
    std::vector<double> edge_lengths = tree_length_dist->Sample(num_external_edges, num_internal_edges);

    //temporary!
    //std::ofstream tmpf("debug-edgelens.txt", std::ios::out | std::ios::app);
    //std::copy(edge_lengths.begin(), edge_lengths.end(), std::ostream_iterator<double>(tmpf, "\t"));
    //tmpf << std::endl;
    //tmpf.close();

    //temporary!
    //std::cerr << "|~~~~> TreeManip::setRandomEdgeLensFromTreeLengthDist <~~~~~|" << std::endl;
    //std::cerr << boost::str(boost::format("|~~~~> num_external_edges = %d <~~~~~|") % num_external_edges) << std::endl;
    //std::cerr << boost::str(boost::format("|~~~~> num_internal_edges = %d <~~~~~|") % num_internal_edges) << std::endl;
    //std::copy(edge_lengths.begin(), edge_lengths.end(), std::ostream_iterator<double>(std::cerr, "|"));
    //std::cerr << "\n|~~~~~> tree length (from vect) = " << std::accumulate(edge_lengths.begin(), edge_lengths.end(), 0.0) << std::endl;

    unsigned ntips = 0;
    unsigned ninternals = 0;
    preorder_iterator nd = tree->begin();
    for (++nd; nd != tree->end(); ++nd)
        {
        if (nd->IsSubroot() || nd->IsTip())
            {
            //std::cerr << boost::str(boost::format("|~~~~> assigning edge length %g to tip node <~~~~~|") % edge_lengths[ntips]) << std::endl;
            nd->SetEdgeLen(edge_lengths[ntips]);
            ntips++;
            }
        else
            {
            //std::cerr << boost::str(boost::format("|~~~~> assigning edge length %g to internal node <~~~~~|") % edge_lengths[num_external_edges + ninternals]) << std::endl;
            nd->SetEdgeLen(edge_lengths[num_external_edges + ninternals]);
            ninternals++;
            }
        }
	tree->hasEdgeLens = true;

    //if (num_internal_edges != ninternals)
    //    {
    //    std::cerr << "num_external_edges = " << num_external_edges << std::endl;
    //    std::cerr << "ntips              = " << ntips << std::endl;
    //
    //    std::cerr << "num_internal_edges = " << num_internal_edges << std::endl;
    //    std::cerr << "ninternals         = " << ninternals << std::endl;
    //    }

	PHYCAS_ASSERT(num_external_edges == ntips);
	PHYCAS_ASSERT(num_internal_edges == ninternals);

    //temporary!
    //std::cerr << "|~~~~~> tree length (from tree) = " << tree->EdgeLenSum() << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs an equiprobable tree topology and sets all edge lengths to 1.0.
*/
void TreeManip::buildEquiprobTree(
  unsigned ntips,                               /**< is the number of tip nodes in the final tree (includes the tip node serving as the root) */
  LotShPtr	rng)                                /**< is the random number generator used to determine the branching order */
	{
	PHYCAS_ASSERT(tree);

    // Empty the tree of everything and simply return if ntips is 0
	tree->Clear();
    tree->DeleteAllStoredNodes();
    //tree.reset(new Tree());
    if (ntips == 0)
        return;

    unsigned tip_number = 0;
    unsigned internal_number = ntips;

    // Create the root node (first tip)
	TreeNode * rootNd = tree->GetNewNode();
    rootNd->SetNodeNum(tip_number++);
    rootNd->SetObservable();
	tree->firstPreorder = rootNd;

	// Create the subroot node (which is the only child of root tip)
	TreeNode * nd = tree->GetNewNode();
    nd->SetNodeNum(tip_number++);
    nd->SetObservable();
    nd->SetEdgeLen(1.0);

    // This vector makes it easy to randomly choose next node to bisect
    std::vector<TreeNode *> receptive;
    receptive.reserve(2*ntips - 3);
    receptive.push_back(nd);

	// Connect the root node to the subroot node
	rootNd->lChild = nd;
	rootNd->rSib = NULL;
	rootNd->par = NULL;
	nd->lChild = NULL;
	nd->rSib = NULL;
	nd->par = rootNd;

	for (unsigned i = 2; i < ntips; ++i)
		{
        unsigned num_receptive_nodes = 2*i - 3;
        PHYCAS_ASSERT(num_receptive_nodes == (unsigned)receptive.size());

		// Choose one of the receptive nodes (i.e. any node except rootNd)
		unsigned j = rng->SampleUInt(num_receptive_nodes);
		nd = receptive[j];

		// Create an internal node on nd's edge
		TreeNode * parent = tree->GetNewNode();
        parent->SetNodeNum(internal_number++);
        parent->SetEdgeLen(1.0);
        receptive.push_back(parent);
		parent->lChild = nd;
		parent->rSib = nd->rSib;
		parent->par = nd->par;
		TreeNode * lSib = FindLeftSib(nd);
		if (lSib != NULL)
			lSib->rSib = parent;
		PHYCAS_ASSERT(nd->par != NULL);
		if (nd->par->lChild == nd)
			nd->par->lChild = parent;
		nd->rSib = NULL;
		nd->par = parent;

        // Create a new tip representing nd's sister
		TreeNode * sister = tree->GetNewNode();
        sister->SetNodeNum(tip_number++);
        sister->SetObservable();
        sister->SetEdgeLen(1.0);
        receptive.push_back(sister);
		nd->rSib = sister;
        sister->par = parent;
		sister->lChild = NULL;
		sister->rSib = NULL;
		}

	tree->RefreshPreorder();
	tree->hasEdgeLens = false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a specialized NNI swap for the following situation (where swap2 is the deepest node in the affected
|	subtree):
|>
|	      swap1   B
|	         \   /
|	          \ /
|	     A     V
|	      \   /
|	       \ /
|	        U
|	       /
|	      /
|	   swap2
|>
|	This move effectively moves V down (carrying B along with it) until V is between U and swap2 (you can also think of
|	it as moving U up - with A attached - until it is between V and swap1). If polytomies are allowed, V carries not
|	only B but any other siblings of B (except swap1) as it slides along the path from swap1 to swap2. Note: because
|	swap2 must be the parent of U, there is no need to actually supply it to this function.
*/
void TreeManip::NNISwapSpecial(
  TreeNode * swap1)	/**< is the node whose parent is slid down one node, carrying all of its children except swap1 */
	{
	PHYCAS_ASSERT(swap1 != NULL);

	TreeNode * V = swap1->GetParent();
	PHYCAS_ASSERT(V != NULL);

	TreeNode * U = V->GetParent();
	PHYCAS_ASSERT(U != NULL);

	TreeNode * R = FindRightmostChild(U);
	if (R == V)
		R = swap1;

	TreeNode * L = U->GetLeftChild();
	if (L == V)
		L = swap1;

	while (U->GetLeftChild() != V)
		LChildToLSib(U, swap1);
	while (V->GetRightSib() != NULL)
		RChildToRSib(U, swap1);
	while (V->GetLeftChild() != L)
		LChildToLSib(V, V);
	while (R->GetRightSib() != NULL)
		RChildToRSib(V, V);

	// This rearrangement invalidates the TreeID and node counts
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps the two nodes specified by swap1 and swap2. Assumes tree, tree->root, swap1 and swap2 are all non-NULL.
|	Polytomies are ok at u or v or both, but assumes that the tree is NOT rooted at either swap1 or swap2.
*/
void TreeManip::NNISwap(
  TreeNode * swap1,		/**< */
  TreeNode * swap2)		/**< */
	{
	PHYCAS_ASSERT(tree != NULL);
	PHYCAS_ASSERT(tree->GetFirstPreorder() != NULL);
	PHYCAS_ASSERT(tree->GetFirstPreorder() != swap1);
	PHYCAS_ASSERT(tree->GetFirstPreorder() != swap2);
	PHYCAS_ASSERT(swap1 != NULL);
	PHYCAS_ASSERT(swap2 != NULL);
	PHYCAS_ASSERT(swap1->GetParent() != NULL);
	PHYCAS_ASSERT(swap2->GetParent() != NULL);

	// Let node v be the higher of the two parent nodes, with u being the lower of the two parent nodes,
	// and let x be v's child and y be u's child:
	//
	//           \ /   /
	//            x   /
	//             \ /
	//        \ /   v
	//         y   /
	//          \ /
	//           u
	//          /
	//
	// Assume swap1 is the higher of the two swap nodes
	//
	TreeNode * x = swap1;
	TreeNode * v = swap1->GetParent();
	TreeNode * y = swap2;
	TreeNode * u = swap2->GetParent();

	if (u->GetParent() == v)
		{
		// swap2 was actually the higher of the two swap nodes
		//
		y = swap1;
		u = swap1->GetParent();
		x = swap2;
		v = swap2->GetParent();
		}

	PHYCAS_ASSERT(v->GetParent() == u);

	// +------------------------------------------------------------------+
	// | First order of business is to save the new navigational pointer  |
	// | targets before we go changing things                             |
	// +------------------------------------------------------------------+

	TreeNode * xlast 					= FindLastPreorderInClade(x);
	TreeNode * ylast 					= FindLastPreorderInClade(y);

	TreeNode * xlastnext 				= xlast->nextPreorder;
	TreeNode * ylastnext 				= ylast->nextPreorder;

	TreeNode * u_lChild					= u->lChild;
	TreeNode * v_lChild					= v->lChild;

	TreeNode * x_rSib					= x->rSib;
	TreeNode * y_rSib					= y->rSib;

	TreeNode * x_par					= x->par;
	TreeNode * y_par					= y->par;

	TreeNode * x_prevPreorder			= x->prevPreorder;
	TreeNode * y_prevPreorder			= y->prevPreorder;

	TreeNode * xlast_nextPreorder		= xlast->nextPreorder;
	TreeNode * ylast_nextPreorder		= ylast->nextPreorder;

	TreeNode * xlsib 					= FindLeftSib(x);
	TreeNode * ylsib 					= FindLeftSib(y);

	TreeNode * xlsiblast				= (xlsib == NULL ? NULL : FindLastPreorderInClade(xlsib));
	TreeNode * ylsiblast				= (ylsib == NULL ? NULL : FindLastPreorderInClade(ylsib));

	// +--------------------------------------------------------------------------+
	// | Now that we are done making pointer reassigments, we can make the switch |
	// +--------------------------------------------------------------------------+

	u->lChild				= (u_lChild == y ? x : u_lChild);
	v->lChild				= (v_lChild == x ? y : v_lChild);

	x->rSib					= y_rSib;
	y->rSib					= x_rSib;

	x->par					= y_par;
	y->par					= x_par;

	x->prevPreorder			= (y_prevPreorder == xlast ? ylast : y_prevPreorder);
	y->prevPreorder			= x_prevPreorder;

	xlast->nextPreorder		= ylast_nextPreorder;
	ylast->nextPreorder		= (xlast_nextPreorder == y ? x : xlast_nextPreorder);

	u->nextPreorder			= (u_lChild == y ? x : u_lChild);
	v->nextPreorder			= (v_lChild == x ? y : v_lChild);

	if (xlsib != NULL)
		{
		xlsib->rSib = y;
		xlsiblast->nextPreorder = y;
		}
	if (ylsib != NULL)
		{
		ylsib->rSib = x;
		if (ylsiblast != xlast)
			ylsiblast->nextPreorder = x;
		}

	if (xlastnext == NULL)
		tree->lastPreorder = ylast;
	else if (xlastnext != y)
		xlastnext->prevPreorder	= ylast;

	if (ylastnext == NULL)
		tree->lastPreorder = xlast;
	else if (ylastnext != x)
		ylastnext->prevPreorder	= xlast;

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes leftmost child of node `u' and makes it the immediate left sibling of `w'. Used in performing Larget-Simon
|	moves when polytomies are present. Does not necessarily leave tree in a consistent state (i.e. tree may contain
|	nodes of order 2 afterwards).
|>
|	Before:                After:
|
|	      c       w              c   a   w
|	       \     /		          \  |  /
|	        \   /		           \ | /
|	         \ /		            \|/
|	  a   b   v			     b       v
|	   \  |  /			      \     /
|	    \ | /			       \   /
|	     \|/			        \ /
|	      u				         u
|	     /				        /
|>
*/
void TreeManip::LChildToLSib(
  TreeNode * u,	/* this node's leftmost child will be removed */
  TreeNode * w)	/* the removed node will become this node's left sibling */
	{
	PHYCAS_ASSERT(u				!= NULL);
	PHYCAS_ASSERT(u->lChild		!= NULL);	// must have a node to mode
	PHYCAS_ASSERT(w				!= NULL);
	PHYCAS_ASSERT(w->par			!= NULL);
	PHYCAS_ASSERT(w->par->lChild	!= NULL);

	// Make copies of the important pointers
	//
	TreeNode * v					= w->par;
	TreeNode * a					= u->lChild;
	TreeNode * a_last				= FindLastPreorderInClade(a);
	TreeNode * a_rSib				= a->rSib;
	TreeNode * a_prevPre			= a->prevPreorder;
	TreeNode * w_prevPre			= w->prevPreorder;
	TreeNode * w_lSib				= FindLeftSib(w);
	TreeNode * w_lSibLast			= (w_lSib == NULL ? NULL : w_prevPre);
	//TreeNode * v_lChild			= v->lChild;

	// Make the switch
	//
	a->rSib							= w;
	a->par							= v;
	a->prevPreorder					= w_prevPre;
	a_last->nextPreorder			= w;
	a_rSib->prevPreorder			= a_prevPre;
	u->lChild						= a_rSib;
	u->nextPreorder					= a_rSib;
	w->prevPreorder					= a_last;
	if (w_lSib == NULL)
		{
		v->lChild					= a;
		v->nextPreorder				= a;
		if (u == w)
			u->lChild				= a_rSib;
		}
	else
		{
		w_lSib->rSib				= a;
		w_lSibLast->nextPreorder	= a;
		}

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes rightmost child of node `u' and makes it the immediate right sibling of w. Used in performing Larget-Simon
|	moves when polytomies are present. Does not necessarily leave tree in a consistent state (i.e. tree may contain
|	nodes of order 2 afterwards). Assumes `w' either equals `u' or is above `u' in the tree as it is currently rooted.
|>
|	Before:                After:
|
|	w       c                w   a   c
|	 \     /		          \  |  /
|	  \   /		               \ | /
|	   \ /		                \|/
|	    v	b  a                 v       b
|	     \  |  /		          \     /
|	      \ | /			           \   /
|	       \|/			            \ /
|	        u	                     u
|	       /				        /
|>
*/
void TreeManip::RChildToRSib(
  TreeNode * u,	/**< this node's rightmost child will be removed */
  TreeNode * w)	/**< the removed node will become this node's right sibling */
	{
	PHYCAS_ASSERT(u				!= NULL);
	PHYCAS_ASSERT(u->lChild		!= NULL);
	PHYCAS_ASSERT(w				!= NULL);
	PHYCAS_ASSERT(w->par			!= NULL);
	PHYCAS_ASSERT(w->par->lChild	!= NULL);

	// Make copies of the important pointers
	//
	TreeNode * w_par					= w->par;
	TreeNode * a						= FindRightmostChild(u);
	//	TreeNode * a_par				= a->GetParent();
	TreeNode * w_last					= FindLastPreorderInClade(w);
	TreeNode * w_last_nextPre			= w_last->nextPreorder;
	TreeNode * a_last					= FindLastPreorderInClade(a);
	TreeNode * a_last_nextPre			= a_last->nextPreorder;
	TreeNode * w_rSib					= w->rSib;
	TreeNode * a_lSib					= FindLeftSib(a);
	TreeNode * a_lSib_last				= FindLastPreorderInClade(a_lSib);

	// Make the switch
	//
	w->rSib								= a;
	a->par								= w_par;
	a->rSib								= w_rSib;
	a_lSib->rSib						= NULL;

	if (u == w)
		{
		a_lSib_last->nextPreorder			= a;
		a->prevPreorder						= a_lSib_last;
		}
	else
		{
		w_last->nextPreorder				= a;
		a->prevPreorder						= w_last;

		if (w_last_nextPre != a)
			{
			w_last_nextPre->prevPreorder	= a_last;
			a_last->nextPreorder			= w_last_nextPre;
			a_lSib_last->nextPreorder		= a_last_nextPre;
			}

		if (a_last_nextPre == NULL && a_lSib_last != w_last)
			tree->lastPreorder				= a_lSib_last;
		else if (w_last_nextPre != a)
			a_last_nextPre->prevPreorder	= a_lSib_last;

		if (w_rSib != NULL)
			w_rSib->prevPreorder			= a_last;
		}

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches node 's' from tree, keeping intact subtree rooted at `s'. Assumes `s' is non-NULL and is neither the root
|	node nor the subroot node (i.e. only child of the root node). Sets all pointers of `s' to NULL except `lChild' and
|	`nextPreorder'.
*/
void TreeManip::DetachSubtree(
  TreeNode * s)	/**< is the root of the subtree to be detached */
	{
	PHYCAS_ASSERT(s != NULL);
	PHYCAS_ASSERT(!s->IsTipRoot());
	//PHYCAS_ASSERT(!s->GetParent()->IsTipRoot());

	// Save pointers to relevant nodes
	//
	TreeNode * s_par				= s->par;
	TreeNode * s_lSib				= FindLeftSib(s);
	TreeNode * s_rSib				= s->rSib;
	TreeNode * s_prevPreorder		= s->prevPreorder;
	TreeNode * slast				= FindLastPreorderInClade(s);
	TreeNode * slast_nextPreorder	= slast->nextPreorder;

	// Completely detach s and seal up the wound
	//
	s->par = NULL;

	s->rSib = NULL;
	if (s_lSib == NULL)
		s_par->lChild = s_rSib;
	else
		s_lSib->rSib = s_rSib;

	s->prevPreorder = NULL;
	if (s_prevPreorder == NULL)
		tree->firstPreorder = slast_nextPreorder;
	else
		s_prevPreorder->nextPreorder = slast_nextPreorder;

	slast->nextPreorder = NULL;
	if (slast_nextPreorder == NULL)
		tree->lastPreorder = s_prevPreorder;
	else
		slast_nextPreorder->prevPreorder = s_prevPreorder;

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes `s' (a sibling of `u') and makes it a child of `u'. If `m' equals `TreeManip::kOnRight', `s' will be added
|	as the immediate right sibling of `targetSib' or, if `targetSib' is not specified, `s' will become the rightmost
|	child of `u'. If `m' is `TreeManip::kOnLeft', `s' will be added as the immediate left sibling of `targetSib' or,
|	if `targetSib' is not specified, `s' will become the leftmost child of `u'. Assumes that `u', `u->par',`s' and
|	`s->par' are all non-NULL and that `u->par' equals `s->par'. Also assumes that targetSib is a child of u if
|	targetSib is specified.
|>
|	Before:                After:
|
|	a       b                s   a   b
|	 \     /		          \  |  /
|	  \   /		               \ | /
|	   \ /		                \|/
|	    u	c   s                u       c
|	     \  |  /		          \     /
|	      \ | /			           \   /
|	       \|/			            \ /
|	        v	                     v
|	       /				        /
|>
*/
void TreeManip::SibToChild(
  TreeNode *u,			/* this node's sibling will be removed and attached as the leftmost child of this node */
  TreeNode *s,			/* the sibling of u to be removed */
  InsertMode ,			/* if kOnRight, s will be added to the right of targetSib or, if targetSib is not specified, s will become rightmost child of u (vice versa if kOnLeft specified) */
  TreeNode *targetSib)	/* if specified, s will become and immediate sibling of this node (which should be a child of u) */
	{
	PHYCAS_ASSERT(u				!= NULL);
	PHYCAS_ASSERT(u->par			!= NULL);
	PHYCAS_ASSERT(s->par			!= NULL);
	PHYCAS_ASSERT(u->par == s->par);
	PHYCAS_ASSERT(targetSib == NULL || (targetSib->par == u));

	DetachSubtree(s);
	InsertSubtree(s, u, TreeManip::kOnRight, targetSib);

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes leaf node `u' and its associated edge.
*/
void TreeManip::DeleteLeaf(
  TreeNode * u,	            /**< this node's edge will be removed */
  bool store_as_internal)   /**< if true, store u in internal node storage; otherwise, store it in leaf node storage */
	{
	// Note which nodes have pointers potentially aimed at u (will need to be redirected)
	TreeNode * u_par		= u->par;
	TreeNode * u_rSib	    = u->rSib;
	TreeNode * u_lSib	    = FindLeftSib(u);
	TreeNode * u_prevPre    = u->prevPreorder;
	TreeNode * u_nextPre    = u->nextPreorder;

	// Redirect pointers aimed at u
	if (u_nextPre == NULL)
		tree->lastPreorder = u_prevPre;
	else
		u_nextPre->prevPreorder = u_prevPre;

	if (u_lSib == NULL)
		u_par->lChild = u_rSib;
	else
		u_lSib->rSib = u_rSib;

	PHYCAS_ASSERT(u_prevPre != NULL);
	u_prevPre->nextPreorder = u_nextPre;

	// Detach and then delete u
	PHYCAS_ASSERT(u->lChild == NULL);
	u->rSib			= NULL;
	u->par			= NULL;
	u->prevPreorder = NULL;
	u->nextPreorder = NULL;

    if (store_as_internal)
	    tree->StoreInternalNode(u);
    else
	    tree->StoreLeafNode(u);

	if (u_par->CountChildren() == 1)
		{
		double edge_len = u_par->GetEdgeLen();
		TreeNode * lchild = u_par->GetLeftChild();
		TreeNode * u_par_lSib	= FindLeftSib(u_par);
		if (u_par_lSib)
			u_par_lSib->rSib = lchild;
		else
			u_par->par->lChild = u_par->lChild;

		lchild->SetEdgeLen(edge_len + lchild->GetEdgeLen());
		lchild->rSib = u_par->rSib;
		lchild->par = u_par->par;
		lchild->prevPreorder = u_par->prevPreorder;
		u_par->prevPreorder->nextPreorder = lchild;

        // Detach and then delete u_par
	    u_par->lChild       = NULL;
	    u_par->rSib			= NULL;
	    u_par->par			= NULL;
	    u_par->prevPreorder = NULL;
	    u_par->nextPreorder = NULL;
		tree->StoreInternalNode(u_par);
		}

	// This rearrangement invalidates the TreeID and node counts
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Builds a star tree with the number of taxa implied by the first element of `split_vect', then applies each aplit
|   from `split_vect' in turn until finished. Assumes all splits are mutually compatible, that none are trivial splits,
|   and that each split is represented by a string of '-' and '*' characters.
*/
void TreeManip::buildTreeFromSplitVector(
    const std::vector<std::string> & split_vect,	/**> is the vector of string representations of splits specifying the edges in the tree */
    ProbDistShPtr edge_len_dist)	                /**< is the probability distribution from which to draw the edge lengths */
	{
	//std::ofstream tmpf("idbuild.txt");

    // Bail out now if user has supplied an empty split_vect
	if (split_vect.empty())
		{
		//tmpf << "split vector empty" << std::endl;
        throw XPhylogeny("could not build tree from empty split vector");
		}

    // Create a vector of split objects from the supplied vector of split representations
    // Determine number of taxa from first element of split_vect.
    std::vector<Split> splits;
    splits.resize(split_vect.size());
    std::vector<std::string>::const_iterator sit = split_vect.begin();
    unsigned ntips = (unsigned)sit->length();
    bool is_rooted = tree->IsRooted();

    unsigned i = 0;
    for (; sit != split_vect.end(); ++sit, ++i)
        {
        // Create a Split object from the pattern representation
        Split & s = splits[i];
        s.CreateFromPattern(*sit);
        PHYCAS_ASSERT(s.GetNTaxa() == ntips);
        PHYCAS_ASSERT(s.CountOnBits() < ntips - 1); // no trivial splits allowed
        PHYCAS_ASSERT(s.CountOnBits() > 1); // no trivial splits allowed
        if (!is_rooted && s.IsBitSet(0))
            s.InvertSplit();    // ensure all splits are oriented correctly if tree is unrooted
        }

    // Sort the split objects from smallest to largest so that more inclusive splits will
    // follow the less inclusive splits
    std::sort(splits.begin(), splits.end());

    // Build a star tree to begin with
    starTree(ntips, edge_len_dist);

    TreeNode * root = tree->GetFirstPreorder();
    TreeNode * subroot = root->GetNextPreorder();
    PHYCAS_ASSERT(root != NULL);
    PHYCAS_ASSERT(subroot != NULL);

    // Loop over all splits
    for (std::vector<Split>::iterator splitit = splits.begin(); splitit != splits.end(); ++splitit)
        {
        Split & s = *splitit;

        // Create a new node to hold the taxa in the current split
       	TreeNode * newNd = tree->GetNewNode();
	    double newNdLen = edge_len_dist->Sample();
	    newNd->SetEdgeLen(newNdLen);

        tree->RecalcAllSplits(ntips);

        if (subroot->GetSplit().SubsumedIn(s))
            continue;

        // Identify all nodes that should be detached, but don't detach them just yet
        std::vector<TreeNode *> prune_list;
        for (TreeNode * nd = subroot->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
            {
            Split & ss = nd->GetSplit();
            if (ss.SubsumedIn(s))
                {
                prune_list.push_back(nd);
                nd = tree->FindLastPreorderInClade(nd);
                }
            }

        PHYCAS_ASSERT(!prune_list.empty());

        // Now detach all nodes in prune_list and add each to newNd

        for (std::vector<TreeNode *>::iterator nit = prune_list.begin(); nit != prune_list.end(); ++nit)
            {
            DetachSubtree(*nit);
            InsertSubtree(*nit, newNd, TreeManip::kOnRight);
            }

        // Add newNd to subroot (only descendant of tip serving as the root)
        InsertSubtree(newNd, subroot, TreeManip::kOnRight);

        tree->DebugCheckTree(false, false, 2);
        }
	}

// *********************************************************************************************************************
// *********************************************************************************************************************
// ************************************ Phycas not yet ready for functions below here **********************************
// *********************************************************************************************************************
// *********************************************************************************************************************

#if 0

//#define DEBUG_BUILD_TREE_FROM_ID
#if defined(DEBUG_BUILD_TREE_FROM_ID)
#	include <iostream>
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Expects splits in tree_id to be ordered from smallest (least inclusive) to largest. All nodes in tree have branch
|	lengths equal to 1.0 after this function. Use SplitManager::SetBrlensFromSplits function to assign branch lengths
|	based on splits in SplitManager. No TreeID is made here, so if SplitManager::SetBrlensFromSplits is not called
|	(it creates a TreeID), you should call CreateID explicitly after building the tree. Assumes tree_id contains at
|	least some splits.
|	Note that there is a nearly identical method called TreeManip::SimpleBuildTreeFromID that assumes that all taxa in
|	the taxa manager are in the tree being built.
*/
void TreeManip::BuildTreeFromID(
  NxsTaxaManager *taxaMgr,	/**< is the taxa manager (knows which taxa are currently active) */
  const TreeID& tree_id,	/**< s the TreeID object specifying the splits to be used in constructing the tree */
  unsigned root_at)			/**< is the index of the tip node serving as the root */	//POL added 22-Oct-2004
	{
	PHYCAS_ASSERT(taxaMgr != NULL);
	unsigned nlvs = taxaMgr->GetNumActive();
	unsigned maxlvs = taxaMgr->GetNumTaxa();
	unsigned nextInternalNodeNum = maxlvs;
	const double commonDefEdgeLen = 1.0;

	SplitSet const &ss = tree_id.GetSplitSet();
	PHYCAS_ASSERT(!ss.empty());

#if defined(DEBUG_BUILD_TREE_FROM_ID)
	std::string tmps;
	std::ofstream tmpf("idbuild.txt");
	tmpf << "tree_id comprises these " << (unsigned)ss.size() << " splits:" << std::endl;
	for (SplitSet::const_iterator ssi = ss.begin(); ssi != ss.end(); ++ssi)
		{
		const Split &s = (*ssi);
		tmps.clear();
		s.CreateAndAppendPatternRepresentation(&tmps);
		tmpf << tmps << "\n";
		}
	tmpf << std::endl;
#endif

	typedef list<TreeNode *> NodeList;
	NodeList nodes;
	const NxsIndexSet &active_taxa = taxaMgr->GetActiveSet();
	NxsIndexSet::const_iterator iter = active_taxa.begin();
	for (; iter != active_taxa.end(); ++iter)
		{
		unsigned leaf = *iter;
		if (leaf != root_at)
			{
			TreeNode *leafNode = tree->CreateTreeNode(leaf, commonDefEdgeLen, true);
#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << "Adding node number " << leaf << " to node list" << std::endl;
#endif
			nodes.push_back(leafNode);
			}
		}

	NodeList::iterator tmp;
	for (SplitSet::const_iterator ssi = ss.begin(); ssi != ss.end(); ++ssi)
		{
		const Split &s = (*ssi);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "\n" << std::setw(6) << "s";
		tmpf << " --> ";
		tmps.clear();
		s.CreateAndAppendPatternRepresentation(&tmps);
		tmpf << tmps << std::endl;
#endif

		TreeNode *newnd = tree->CreateTreeNode(nextInternalNodeNum++, commonDefEdgeLen, false);
		NodeList::iterator it = nodes.begin();
		for (; it != nodes.end();)
			{
			TreeNode *child = (*it);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << std::setw(6) << child->GetNodeNumber() << " --> ";
			tmps.clear();
			child->split.CreateAndAppendPatternRepresentation(&tmps);
			tmpf << tmps;
#endif

			bool subsumed = child->split.SubsumedIn(s);
			if (subsumed)
				{
				InsertSubtree(child, newnd, TreeManip::kOnRight);
				newnd->split.CombineWith(child->split);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
				tmpf << " ==> subsumed in s, deleting node";
#endif

				tmp = it++;
				nodes.erase(tmp);
				}
			else
				++it;

#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << std::endl;
#endif
			}

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "\nAdding node number " << newnd->GetNodeNumber() << " to node list" << std::endl;
#endif

		nodes.push_back(newnd);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "Node list now: ";
		for (tmp = nodes.begin(); tmp != nodes.end(); ++tmp)
			{
			tmpf << (*tmp)->GetNodeNumber() << " ";
			}
		tmpf << "\n" << std::endl;
#endif
		}

#if defined(DEBUG_BUILD_TREE_FROM_ID)
	tmpf << "\n\nAll remaining nodes ";
	for (NodeList::iterator it = nodes.begin(); it != nodes.end(); ++it)
		{
		tmpf << "[" << (*it)->GetNodeNumber() << "]";
		}
	tmpf << std::endl;
	//tmpf << " will be added to root (taxon 0):\n  ";
	tmpf << " will be added to root (taxon " << root_at << "):\n  "; //POL 22-Oct-2004

	tmpf.close();
#endif

	tmp = nodes.begin();
	TreeNode *lastNode = *tmp;
	//TreeNode *rootNode = tree->CreateTreeNode(0,0.0, false);
	TreeNode *rootNode = tree->CreateTreeNode(root_at, 0.0, false);	//POL 22-Oct-2004
	InsertSubtree(lastNode, rootNode, TreeManip::kOnRight);

	tree->firstPreorder = rootNode;
	tree->TraverseTree();

	//if (Tree::gDebugOutput)
    //	{
	//	std::cerr << "GetNLeaves() = " << tree->GetNLeaves() << std::endl;
	//	std::cerr << "nlvs = " << nlvs << std::endl;
	//	}
	PHYCAS_ASSERT(tree->GetNLeaves() == nlvs);
	PHYCAS_ASSERT(rootNode == tree->GetFirstPreorder());

	tree->RefreshFullID();
	tree->DebugCheckTreeStructure();
	}

#endif
