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

#include <string>
#include <cstdio>
#include <cctype>
#include "basic_tree.hpp"
#include "tree_manip.hpp"
#include "xphylogeny.hpp"

#include <boost/format.hpp>

// these were at the top of basic_tree.inl
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

//#define TESTING_TOWARD_NODE_ITERATOR
#if defined(TESTING_TOWARD_NODE_ITERATOR)
#	include "edge_iterators.hpp"
#endif
using namespace phycas;
//bool Tree::gDebugOutput = true;

#define INTERNAL_NODENUM_OFFSET 1000

TreeNode * Tree::FindNodeBySplit(const Split & s)
    {
    std::vector<TreeNode *> leafVec;
    leafVec.resize(nTips);
    preorder_iterator ndIt = begin();
    for (; ndIt != end(); ++ndIt)
        {
        TreeNode & nd = (*ndIt);
        if (nd.split == s)
            return &nd;
        }
    return 0L;
    }

/// Get vector of tips sorted by the node number
std::vector<TreeNode *> Tree::GetTips()
    {
    std::vector<TreeNode *> leafVec;
    leafVec.resize(nTips);
    preorder_iterator ndIt = begin();
    for (; ndIt != end(); ++ndIt)
        {
        TreeNode & nd = (*ndIt);
        if (nd.IsTip())
            {
            const unsigned index = nd.GetNodeNumber();
            if (index >= leafVec.size())
                leafVec.resize(index + 1);
            leafVec[index] = &nd;
            }
        }
    return leafVec;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets firstPreorder to NULL and calls Clear() to initialize data members.
*/
Tree::Tree()
  : firstPreorder(NULL)
	{
	//std::cerr << "Tree::Tree()" << std::endl;
	isRooted = false;
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Clear() to store all nodes and delete all other allocated memory, then erases `internalNodeStorage' and
|   `tipStorage' to eliminate the memory allocated for TreeNode objects.
*/
Tree::~Tree()
	{
	//std::cerr << "Tree::~Tree()" << std::endl;//tmp
	Clear();
    DeleteAllStoredNodes();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all nodes stored in `tipStorage' and `internalNodeStorage'.
*/
void Tree::DeleteAllStoredNodes()
    {
	while (!internalNodeStorage.empty())
	    {
	    //std::cerr << "Popping internal node" << std::endl;
	    TreeNode * nd = internalNodeStorage.top();
		internalNodeStorage.pop();
		delete nd;
		nd = NULL;
		}
	while (!tipStorage.empty())
		{
	    //std::cerr << "Popping tip node" << std::endl;
	    TreeNode * nd = tipStorage.top();
		tipStorage.pop();
		delete nd;
		nd = NULL;
		}
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Clears each node in the tree and sets `firstPreorder' to NULL, thus forgetting the tree structure. This is useful
|   to prevent the tree destructor from deleting nodes that are owned by other objects. Identical to Clear except for
|   handling of `firstPreorder'.
*/
void Tree::Forget()
	{
    Clear(true);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns Tree object to just-constructed state. If preorder pointers are valid, walks tree in postorder fashion
|	storing nodes as they are visited in `internalNodeStorage' or `tipStorage' (nodes are cleared but not stored if
|   `forget' is true).
*/
void Tree::Clear(
  bool forget)  /**< if true, nodes will be cleared but not stored; if false, nodes will be stored */
	{
	if (firstPreorder != NULL)
		{
		// Walk tree in postorder fashion, storing nodes as we go
		for (TreeNode *nd = GetLastPreorder(); nd != NULL;)
			{
			TreeNode * next = nd->GetNextPostorder();
            if (!forget)
                {
                if (nd->IsTip())
                    StoreLeafNode(nd);
                else
                    StoreInternalNode(nd);
                }
			nd->Clear();
			nd = next;
			}

		firstPreorder = NULL;
		}

	lastPreorder		= NULL;
	nInternals			= 0;
	nTips				= 0;
	preorderDirty		= true;
	hasEdgeLens			= false;
	nodeCountsValid		= false;
	treeid_valid		= false;
	numbers_from_names	= false;
    debugOutput         = false;
   	}

/*----------------------------------------------------------------------------------------------------------------------
|	Push node onto end of `internalNodeStorage' vector (nodes are stored rather than deleted to save having to reallocate them
|	later.
*/
TreeNode * Tree::PopLeafNode()
	{
    PHYCAS_ASSERT(!tipStorage.empty());
	TreeNode * nd = tipStorage.top();
	tipStorage.pop();
	nd->par = nd->lChild = nd->rSib = nd->nextPreorder = nd->prevPreorder = NULL;
    return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Pop node off end of `internalNodeStorage' vector
*/
TreeNode * Tree::PopInternalNode()
	{
    PHYCAS_ASSERT(!internalNodeStorage.empty());
	TreeNode * nd = internalNodeStorage.top();
	internalNodeStorage.pop();
	nd->par = nd->lChild = nd->rSib = nd->nextPreorder = nd->prevPreorder = NULL;
    return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Push node onto end of `internalNodeStorage' vector (nodes are stored rather than deleted to save having to reallocate them
|	later.
*/
void Tree::StoreInternalNode(TreeNode * u)
	{
	internalNodeStorage.push(u);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Push node onto end of `tipStorage' vector.
*/
void Tree::StoreLeafNode(TreeNode * u)
	{
	tipStorage.push(u);
    //std::cerr << boost::str(boost::format("tipStorage now has %d elements") % (unsigned)(tipStorage.size())) << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Should only be called from within RerootAtThisTip() or RerootAtThisInternal. Makes node `m' (the mover) the
|	rightmost child of node `t' (the target). Assumes both `m' and `t' are non-NULL. Because the purpose of this
|	function is to move nodes below the new root node to the appropriate point above the new root node, this method
|	also assumes that `t' is always a descendant (perhaps remote) of `m'.
*/
void Tree::RerootHelper(TreeNode * m, TreeNode * t)
	{
	PHYCAS_ASSERT(m != NULL);
	PHYCAS_ASSERT(t != NULL);

	// Save nodes to which m attaches
	TreeNode * m_lChild	= m->lChild;
	TreeNode * m_rSib	= m->rSib;
	TreeNode * m_par		= m->par;

	// Starting with t, walk down tree to identify x, the child of m that is on the path from m to t
	TreeNode * x = t;
	while (x->par != m)
		{
		x = x->par;
		PHYCAS_ASSERT(x != NULL);
		}
	TreeNode * x_rSib = x->rSib;

	// identify x_lSib, the immediate left sibling of x (will be NULL if x is lChild of m)
	TreeNode * x_lSib = NULL;
	if (x != m_lChild)
		{
		x_lSib = m_lChild;
		while (x_lSib->rSib != x)
			{
			x_lSib = x_lSib->rSib;
			PHYCAS_ASSERT(x_lSib != NULL);
			}
		}

	// identify m_lSib, the immediate left sibling of m (will be NULL if m is root node or is lChild of its parent)
	TreeNode *m_lSib = NULL;
	if (m_par != NULL && m != m_par->lChild)
		{
		m_lSib = m_par->lChild;
		while (m_lSib->rSib != m)
			{
			m_lSib = m_lSib->rSib;
			PHYCAS_ASSERT(m_lSib != NULL);
			}
		}

	// Put x where m is now
	preorderDirty = true;
	if (m_par == NULL)
		{
		// m is the root node
		PHYCAS_ASSERT(m_rSib == NULL);
		PHYCAS_ASSERT(m_lSib == NULL);
		x->rSib = NULL;
		x->par = NULL;
		if (x == m_lChild)
			m->lChild = x_rSib;
		else
			x_lSib->rSib = x_rSib;
		}
	else if (m == m_par->lChild)
		{
		// m is leftmost child of its parent
		x->rSib = m_rSib;
		x->par = m_par;
		m->rSib = NULL;
		m->par = NULL;
		m_par->lChild = x;
		if (x == m_lChild)
			m->lChild = x_rSib;
		else
			x_lSib->rSib = x_rSib;
		}
	else
		{
		// m is not leftmost child of its parent
		m_lSib->rSib = x;
		x->rSib = m_rSib;
		x->par = m_par;
		m->rSib = NULL;
		m->par = NULL;
		if (x == m_lChild)
			m->lChild = x_rSib;
		else
			x_lSib->rSib = x_rSib;
		}

	// Make m the new rightmost child of t
	m->par = t;
	if (t->lChild == NULL)
		t->lChild = m;
	else
		{
		// Find rightmost child of t
		m_lSib = t->lChild;
		while (m_lSib->rSib != NULL)
			m_lSib = m_lSib->rSib;

		// Make rightmost child of t the left sib of m
		m_lSib->rSib = m;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reroots the tree at the node specified by the TreeNode pointer `nd', which is assumed to be non-NULL and an internal
|	node.
*/
void Tree::RerootAtThisInternal(TreeNode * nd)
	{
	//POL: could be merged with RerootAtThisTip, but note that the two use different root test functions: IsTipRoot vs. IsInternalRoot
	PHYCAS_ASSERT(!nd->IsTip());
	PHYCAS_ASSERT(nd != NULL);
	TreeNode * t = nd;
	TreeNode * m = nd->par;
	while (!nd->IsInternalRoot())
		{
		//std::cerr << "In RerootAt: t = " << t->GetNodeName() << ", m = " << m->GetNodeName() << std::endl; //POL-debug
		PHYCAS_ASSERT(nd->HasParent());

		// Begin by swapping the mover's edge length with nd's edge length
		if (HasEdgeLens())
			{
			double tmp_edgelen = m->edgeLen;
			m->edgeLen = nd->edgeLen;
			nd->edgeLen = tmp_edgelen;
			}

		// Make the "mover" node m (along with all of its children except nd)
		// the rightmost child of the "target" node t
		RerootHelper(m, t);

		// The next target node is always the previous mover, and the next mover node is always nd's parent
		t = m;
		m = nd->par;
		}

	firstPreorder = nd;
	nd->prevPreorder = NULL;

	//std::cerr << "\n\nWalking tree just before leaving RerootAt...\n"; //POL-debug
	//std::cerr << DebugWalkTree(true, 2) << std::endl; //POL-debug
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reroots the tree at the node specified by the TreeNode pointer `nd', which is assumed to be non-NULL and a tip node.
*/
void Tree::RerootAtThisTip(TreeNode * nd)
	{
	//POL: could be merged with RerootAtThisInternal, but note that the two use different root test functions: IsTipRoot vs. IsInternalRoot
	PHYCAS_ASSERT(nd->IsTip());
	PHYCAS_ASSERT(nd != NULL);
	TreeNode * t = nd;
	TreeNode * m = nd->par;
	while (!nd->IsTipRoot())
		{
		//std::cerr << "In RerootAt: t = " << t->GetNodeName() << ", m = " << m->GetNodeName() << std::endl; //POL-debug
		PHYCAS_ASSERT(nd->HasParent());

		// Begin by swapping the mover's edge length with nd's edge length
		if (HasEdgeLens())
			{
			double tmp_edgelen = m->edgeLen;
			m->edgeLen = nd->edgeLen;
			nd->edgeLen = tmp_edgelen;
			}

		// Make the "mover" node m (along with all of its children except nd)
		// the rightmost child of the "target" node t
		RerootHelper(m, t);

		// The next target node is always the previous mover, and the next mover node is always nd's parent
		t = m;
		m = nd->par;
		}

	firstPreorder = nd;
	nd->prevPreorder = NULL;

	//std::cerr << "\n\nWalking tree just before leaving RerootAt...\n"; //POL-debug
	//std::cerr << DebugWalkTree(true, 2) << std::endl; //POL-debug
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Searches for tip (degree = 1) node whose `nodeNum' data member equals `num', then reroots the tree at that node. If
|	a node by that number cannot be found, throws an XPhylogeny exception.
*/
void Tree::RerootAtTip(unsigned num)
	{
	//std::cerr << "*** in RerootAtTip ***" << std::endl;

	TreeNode *nd = FindTipNode(num);
	if (nd == NULL)
		{
		workspace.clear();
		workspace << boost::str(boost::format("there is no tip node having number %d") % num);
		throw XPhylogeny(workspace);
		//throw XPhylogeny(append_unsigned(workspace, num));
		}
	else
		RerootAtThisTip(nd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes preorder pointers (`Tree::firstPreorder', `Tree::lastPreorder', `TreeNode::nextPreorder',
|	`TreeNode::prevPreorder') by starting with the root node (`firstPreorder') and walking through the tree by visiting
|	`lChild', `rSib' and `par' pointers. Sets `preorderDirty' to false when done.
*/
void Tree::RefreshPreorder(TreeNode *nd) const
	{
	if (!preorderDirty)
		return;

	if (nd == NULL)
		nd = firstPreorder;
	if (nd == NULL)
		return; // tolerate calling this on a empty tree

	if (nd->lChild == NULL && nd->rSib == NULL)
		{
		// nd has no children and no siblings, so nextPreorder is the right sibling of
		// the first ancestral node that has a right sibling.
		TreeNode *anc = nd->par;
		while ((anc != NULL) && (anc->rSib == NULL))
			anc = anc->par;
		if (anc == NULL)
			{
			// We descended all the way to the root without finding an ancestor with
			// a right sibling, so nd must be the upper-right-most node in the tree
			lastPreorder = nd;
			nd->nextPreorder = NULL;
			}
		else
			{
			// We found an ancestor with a right sibling without having to go all
			// the way to the root of the tree
			nd->nextPreorder = anc->rSib;
			anc->rSib->prevPreorder = nd;
			}
		}
	else if (nd->lChild == NULL && nd->rSib != NULL)
		{
		// nd has no children (it is a tip), but does have a sibling on its right
		nd->nextPreorder = nd->rSib;
		nd->rSib->prevPreorder = nd;
		}
	else if (nd->lChild != NULL && nd->rSib == NULL)
		{
		// nd has children (it is an internal node) but no siblings on its right
		nd->nextPreorder = nd->lChild;
		nd->lChild->prevPreorder = nd;
		}
	else
		{
		// nd has both children and siblings on its right
		nd->nextPreorder = nd->lChild;
		nd->lChild->prevPreorder = nd;
		}

	if (nd->lChild != NULL)
		{
		RefreshPreorder(nd->lChild);
		}

	if (nd->rSib != NULL)
		{
		RefreshPreorder(nd->rSib);
		}

	if (nd == firstPreorder)
		preorderDirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets internal node numbers starting with the supplied value `start_at'. Visits nodes in preorder sequence, so the
|   subroot node (the only descendant of the tip node serving as the root) is assigned the number `start_at'.
*/
void Tree::renumberInternalNodes(
  unsigned start_at)    /**< is the number to use for the first internal node encountered in the preorder sequence */
	{
    unsigned num = start_at;
    for (preorder_iterator nd = begin(); nd != end(); ++nd)
        {
        if (nd->IsInternal())
            nd->SetNodeNum(num++);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all node names to the empty string "". This forces MakeNewick, for example, to use node numbers as the names of
|   tip nodes.
*/
void Tree::stripNodeNames()
	{
    for (preorder_iterator nd = begin(); nd != end(); ++nd)
        {
        nd->SetNodeName("");
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	The supplied `name_vector' should be a vector of observable tip names. Assumes that observable nodes in the tree
|	already have names, and the goal is to number them appropriately (i.e. set their `nodeNum' to the position of their
|	name in the supplied `name_vector'. If an observable node is encountered that has no name, or if the number of
|	observable nodes is not equal to the length of `name_vector', an XPhylogeny exception is thrown. Depends on the
|	data member `observable' having been set correctly for all nodes in the tree.
*/
void Tree::RectifyNumbers(std::vector<std::string> name_vector)
	{
	unsigned name_vector_len = (unsigned)name_vector.size();
	bool ok = (GetNObservables() == name_vector_len);
	if (!ok)
		throw XPhylogeny("number of names must equal number of observable tips");

	if (preorderDirty)
		RefreshPreorder();

	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsObservable())
			{
			std::vector<std::string>::iterator it = std::find(name_vector.begin(), name_vector.end(), nd->GetNodeName());
			if (it == name_vector.end())
				{
				std::string s = "node name (";
				s << nd->GetNodeName();
				s << ") not found in supplied list of names";
				throw XPhylogeny(s);
				}
			else
				{
				unsigned n = (unsigned)(it - name_vector.begin());
				nd->SetNodeNum(n);
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The supplied `name_vector' should be a vector of tip names. Assumes that the tip nodes have numbers but that their
|	names should be set according to the values supplied in `name_vector'; that is, a node whose number is k will have
|	its name set to `name_vector'[k]. Any tip having the default number TreeNode::nodeNumInitValue will cause an
|	XPhylogeny exception to be thrown. An XPhylogeny exception will also be thrown if the number of tips in the tree
|	does not match the length of `name_vector'. For rooted trees, no attempt is made to look up the root node in the
|   `name_vector' because IsObservable() returns false for the root tip in a rooted tree.
*/
void Tree::RectifyNames(std::vector<std::string> name_vector)
	{
	//std::cerr << "*** in RectifyNames ***" << std::endl;

	unsigned name_vector_len = (unsigned)name_vector.size();
	bool ok = (GetNObservables() == name_vector_len);
	if (!ok)
		throw XPhylogeny(str(boost::format("number of names (%d) not equal to the number of tips (%d)") % name_vector_len % GetNObservables()));

	if (preorderDirty)
		RefreshPreorder();

	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsObservable())
			{
			unsigned n = nd->GetNodeNumber();
			if (n < name_vector_len)
				{
                //std::cerr << boost::str(boost::format("*** setting name for node %d to \"%s\"") % n % name_vector[n]) << std::endl;
				nd->SetNodeName(name_vector[n]);
				}
			else
				{
				std::string s = "node number (";
				s << n;
				s << ") out of range";
				throw XPhylogeny(s);
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `adhocTreeNodePtr' data member to 0 for all nodes in the tree.
*/
void Tree::clearAdHocTreeNodePointers()
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		nd->adhocTreeNodePtr = 0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to false for all nodes in the tree.
*/
void Tree::UnselectAllNodes()
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		nd->UnselectNode();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to true for all nodes in the tree.
*/
void Tree::SelectAllNodes()
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		nd->SelectNode();
		}
	}

#if 0 // including this function results in this VC error: fatal error C1055: compiler limit : out of keys
/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to TreeNode representing the deepest point of intersection between the lineage leading to `tip1' and
|	the lineage leading to `tip2'. Both `tip1' and `tip2' are assumed to be tip nodes, but one or the other could
|	currently be serving as the root. If either `tip1' or `tip2' is serving as the root, then the node returned is the
|	immediate descendant of this node. Otherwise, the chain of parents of both tips is walked until the other tip shows
|	up as a member of the split of the other.
*/
TreeNode * Tree::FindMRCA(unsigned tip1, unsigned tip2)
	{
	PHYCAS_ASSERT(0);	//POL not yet written
	return NULL;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Zig-zags up the right side of the clade starting with the node `start' until reaching a tip. This tip node is the
|	last node in the clade according to the preorder sequence. Zig-zagging is necessary because there is no pointer to
|	the right-most child of a node, so we must use FindRightmostChild to move up the clade. Assumes `start' is non-NULL.
*/
TreeNode * Tree::FindLastPreorderInClade(
  TreeNode * start)	/**< is the deepest node belonging to the clade in question */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
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
|	Begins with left child of parent of `start' and calls GetRightSib() until the left sibling of `start' is located.
|	Assumes `start' is non-NULL.
*/
TreeNode * Tree::FindLeftSib(
  TreeNode * start)	/**< is the node whose left sibling is sought */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(start != NULL);
	TreeNode * nd = start->GetParent();

	// If start has no parent, then there's no way it can have a left sibling
	if (nd == NULL)
		return NULL;

	nd = nd->GetLeftChild();

	// Parent of start should have at least one child
	PHYCAS_ASSERT(nd != NULL);

	// If left child of start's parent equals start, then start is an only child and has no siblings
	if (nd == start)
		return NULL;

	TreeNode * leftsib = NULL;
	while (nd != start)
		{
		if (nd == NULL)
			{
			throw XPhylogeny("pointer inconsistency in FindLeftSib");
			}
		leftsib = nd;
		nd = nd->GetRightSib();
		}
	return leftsib;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Begins with left child of `start'. If left child is NULL, returns NULL, otherwise, returns the rightmost sibling of
|	the left child of `start'. Assumes `start' is non-NULL.
*/
TreeNode * Tree::FindRightmostChild(
  TreeNode * start)	/**< is the parent node whose children will be searched */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(start != NULL);
	TreeNode * curr = start->GetLeftChild();
	TreeNode * rightmost = NULL;
	while (curr != NULL)
		{
		rightmost = curr;
		curr = curr->GetRightSib();
		}
	return rightmost;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts node 's' as child of 'u' keeping the subtree rooted at `s' intact. Assumes both `s' and `u' are non-NULL. If
|	`targetSib' is specified, `s' will become an immediate sibling of `targetSib' (right or left sibling determined by
|	`on_right'). If `targetSib' is not specified, `s' will be added as the rightmost or leftmost child of `u' (again,
|	depending on the value of `on_right').
*/
void Tree::InsertSubtree(
  TreeNode * s,			/**< is the node to be added as a child of u */
  TreeNode * u,			/**< is the new parent of s */
  bool on_right,		/**< if true, `s' will be added to the right of `targetSib' or, if` targetSib' is not specified, `s' will become rightmost child of `u' (vice versa if false specified) */
  TreeNode * targetSib)	/**< if non-NULL, `s' will become an immediate sibling of this node (right or left depending on value of `on_right') */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
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
		if ((!on_right) && (targetSib == u_lChild))
			targetSib = NULL;
		if ((on_right) && (targetSib == u_rChild))
			targetSib = NULL;
		if ((!on_right) && (targetSib != NULL))
			{
			targetSib = FindLeftSib(targetSib);
			on_right = true;
			}
		}

	// Make s a child of u
	//
	if (targetSib == NULL)
		{
		if (on_right)
			{
			s->rSib				= NULL;
			s->par				= u;
			s->prevPreorder 	= ulast;
			slast->nextPreorder = ulast_nextPreorder;
			ulast->nextPreorder	= s;

			if (ulast_nextPreorder == NULL)
				this->lastPreorder = slast;
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

		// Can assume that targetSib is not rightmost child of u and also that on_right == true
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
	this->InvalidateTreeID();
	this->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches node 's' from tree, keeping intact subtree rooted at `s'. Assumes `s' is non-NULL and is neither the root
|	node nor the subroot node (i.e. only child of the root node). Sets all pointers of `s' to NULL except `lChild' and
|	`nextPreorder'.
*/
void Tree::DetachSubtree(
  TreeNode * s)	/**< is the root of the subtree to be detached */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(s != NULL);
	PHYCAS_ASSERT(!s->IsTipRoot());
	//PHYCAS_ASSERT(!s->GetParent()->IsTipRoot());
	//if (s->GetParent()->IsTipRoot())
	//	{
	//	std::cerr << "oops" << std::endl;
	//	}

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
		this->firstPreorder = slast_nextPreorder;
	else
		s_prevPreorder->nextPreorder = slast_nextPreorder;

	slast->nextPreorder = NULL;
	if (slast_nextPreorder == NULL)
		this->lastPreorder = s_prevPreorder;
	else
		slast_nextPreorder->prevPreorder = s_prevPreorder;

	// This rearrangement invalidates the TreeID and node counts
	//
	this->InvalidateID();
	this->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Traverses the tree in postorder fashion looking for a node having only one child. If such a node is found, detaches
|	the child and uses it to replace the parent node, thus eliminating the degree-2 node. Returns the number of such
|	degree-2 nodes found (in case caller wishes to warn user if this number is greater than 1, which probably signifies
|	an error).
*/
unsigned Tree::deroot()
	{
	unsigned number_of_incisions = 0;
	TreeNode * rooted_at = GetFirstPreorder();
	PHYCAS_ASSERT(rooted_at != NULL);
	TreeNode * start_at = rooted_at->GetNextPreorder();
	for (TreeNode * nd = start_at; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsInternal() && nd->CountChildren() == 1)
			{
			TreeNode * nd_rSib   = nd->rSib;
			TreeNode * nd_parent = nd->GetParent();
			TreeNode * nd_child  = nd->GetLeftChild();

			PHYCAS_ASSERT(nd_parent != NULL);
			PHYCAS_ASSERT(nd_child != NULL);

            if (nd_parent == rooted_at)
                {
                // Redundant node is the subroot node (only child of root node)
                TreeManip::InsertMode insert_where = TreeManip::kOnLeft;
                TreeManip tm(shared_from_this());
                tm.DetachSubtree(nd_child);
                TreeNode * target_sib = nd;
                tm.InsertSubtree(nd_child, nd_parent, insert_where, target_sib);
                tm.DetachSubtree(nd);
                StoreInternalNode(nd);
                }
            else
                {
                TreeManip::InsertMode insert_where = TreeManip::kOnLeft;
                TreeNode * target_sib = nd_rSib;
                if (nd_rSib == NULL)
                    {
                    // nd is rightmost child, insert nd_child to right of nd's left sibling
                    insert_where = TreeManip::kOnRight;
                    target_sib = nd_parent->GetLeftChild();
                    while (target_sib->rSib != nd)
                        {
                        target_sib = target_sib->GetRightSib();
                        }
                    }

                TreeManip tm(shared_from_this());
                tm.DetachSubtree(nd_child);
                tm.DetachSubtree(nd);
                StoreInternalNode(nd);
                tm.InsertSubtree(nd_child, nd_parent, insert_where, target_sib);
                }

			nd = nd_parent;
			++number_of_incisions;
			}
		}

	DebugCheckTree(false, true, 2);
	return number_of_incisions;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a postorder traversal and, for each internal node visited, sorts clades from smallest to largest. If the
|	`right' argument is true, then the largest clades will be on the right after ladderization. If the `right' argument
|	is false, then the leftmost child will be the largest clade.
*/
void Tree::Ladderize(
  bool right)	/**> if true, largest clades will be on right; if false, largest clades will be on left */
	{
	for (TreeNode * nd = GetLastPreorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		bool istip = nd->IsTip();
		TreeNode * par = nd->GetParent();
		TreeNode * rsib = nd->GetRightSib();

		if (par)
			{
			if (istip)
				{
				if (rsib)
					par->tmp += 1.0;
				else
					par->tmp = 1.0;
				}
			else
				{
				if (rsib)
					par->tmp += nd->tmp;
				else
					par->tmp = nd->tmp;
				}
			}

		// If node is an internal node, then reorder its children from smallest (left) to
		// largest (right) based on the tmp data member of each child
		if (!istip)
			{
			typedef std::multimap< unsigned, TreeNode * > ChildMap;
			typedef std::pair< unsigned, TreeNode * > ChildMapKeyValuePair;
			ChildMap childmap;
			TreeNode * child = nd->GetLeftChild();
			while (child != NULL)
				{
				unsigned clade_size = (unsigned)child->tmp;
				childmap.insert(ChildMapKeyValuePair(clade_size, child));
				TreeNode * next_child = child->GetRightSib();
				DetachSubtree(child);
				child = next_child;
				}

			// The childmap is ordered from smallest to largest already, so we need only
			// to march down the map adding nodes as we go to the right of the last one
			if (right)
				{
				for (ChildMap::iterator childiter = childmap.begin(); childiter != childmap.end(); ++childiter)
					{
					InsertSubtree(childiter->second, nd, true);
					}
				}
			else
				{
				for (ChildMap::reverse_iterator childiter = childmap.rbegin(); childiter != childmap.rend(); ++childiter)
					{
					InsertSubtree(childiter->second, nd, true);
					}
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Follows preorder pointers until a tip (degree = 1) node is encountered with `nodeNum' equal to `num'. Returns a
|	pointer to this node, or 0 if a tip node numbered `num' cannot be not found.
*/
TreeNode *Tree::FindTipNode(unsigned num)
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip() && num == nd->GetNodeNumber())
			return nd;
		}

	return 0;
	}

#if defined(TESTING_TOWARD_NODE_ITERATOR)
bool isValid(const TreeNode * refNd, const TreeNode * neighborCloserToEffectiveRoot)
	{
	return refNd->IsTip();
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes and returns the longest path from the subroot node to any tip node in the tree, including the tip serving
|   as the root.
*/
double Tree::calcTotalHeight()
	{
    if (!hasEdgeLens)
        throw XPhylogeny("asked to calculate total height of a tree with no edge lengths");

	//@POL: this function could be made more general (i.e. right now it depends on the root being a tip)
	PHYCAS_ASSERT(GetFirstPreorder() != NULL);
	PHYCAS_ASSERT(GetFirstPreorder()->GetNextPreorder() != NULL);

	// Create a map (node number -> height) that will keep track of the maximum height above each internal node
	typedef std::map<unsigned, double> HeightMap;
	HeightMap hmap;

	// Walk through nodes in postorder fashion, skipping only the root tip node and its only child, the subroot node.
    // For each tip node, record the tip's edge length in hmap[parent] if hmap[parent] does not yet exist or
    // if the edge length is greater than the value already stored in hmap[parent]. Do the same for internal
    // nodes, except the "edge length" for an internal node is its edge length plus its hmap value.
    TreeNode * nd = GetLastPreorder();
    PHYCAS_ASSERT(nd->GetParent() != NULL);
	for (; nd->GetParent()->IsInternal(); nd = nd->GetNextPostorder())
		{
        // Note: the only node whose parent is not an internal node is the subroot node
        unsigned ndnum = nd->GetNodeNumber();
        TreeNode * par = nd->GetParent();
        unsigned par_ndnum = par->GetNodeNumber();

        double max_above = nd->GetEdgeLen();
        if (nd->IsInternal())
            max_above += hmap[ndnum];

   		// If parent's hmap entry does not yet exist, add it and make it equal to max_above. If it does exist,
        // make it equal to max_above if max_above is greater than the current value
        // (see item 24, p. 110, in Meyers' Efficient STL for more info on the technique used here)
		HeightMap::iterator lowb = hmap.lower_bound(par_ndnum);
		if (lowb != hmap.end() && !(hmap.key_comp()(par_ndnum, lowb->first)))
			{
			// height entry exists already for this node
            if (lowb->second < max_above)
                lowb->second = max_above;
			}
		else
			{
			// height entry does not yet exist for this node
			hmap.insert(lowb, HeightMap::value_type(par_ndnum, max_above));
			}
		}

    // Maximum tree height is now equal to the maximum height above the subroot
    PHYCAS_ASSERT(nd == GetFirstPreorder()->GetNextPreorder()); // check to make sure that nd is the subroot
	unsigned subroot_ndnum = nd->GetNodeNumber();
	double subroot_height = nd->GetEdgeLen();
	double height_above_subroot = hmap[subroot_ndnum];

    // hmap[subroot_ndnum] should hold the height of the tree above the subroot node. The only exception
    // would be if the root tip edge (which is actually subroot's edge) was itself longer than hmap[subroot_ndnum]
    // in which case the tip root's edge length should be returned
    return (subroot_height > height_above_subroot ? subroot_height : height_above_subroot);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Walks tree using preorder pointers, returning a string documenting the journey. For this tree description:
|	"(c:0.4,d:0.5,(a:0.1, b:0.2):0.3);", the returned string would be: "[1] -> c (0) -> d (1) -> [0] -> a (2) -> b (3)"
|	assuming that the `preorder' parameter was true. Note that node names are output if they were supplied, the internal
|	node numbers are inside square brackets, and the tip node numbers are inside parentheses. The tip nodes are, by
|	default, numbered according to the preorder sequence (and the internal nodes are numbered in postorder sequence).
|	If a vector of tip names was supplied using RectifyNames(), then the tip node numbers will be equal to the 0-offset
|	index of the tip into that vector. This function is intended to be used for debugging purposes and is thus not
|	designed to be particularly fast. Verbosity levels: 0 means minimal (names for named nodes, numbers otherwise);
|	1 means add node numbers and indicate using brackets or parentheses whether node is an [internal], (tip) or {root}
|	node; 2 means show nearly everything about each node (each node requires a paragraph to explain).
*/
std::string Tree::DebugWalkTree(bool preorder, unsigned verbosity)
	{
	std::string s;
#	if defined(TESTING_TOWARD_NODE_ITERATOR)
	NodeValidityChecker validFunctor = isValid;
	for (TreeNode * z = GetFirstPreorder(); z != NULL; z = z->GetNextPreorder())
		{
		if (z->IsTip())
			continue;
		s << "\nStarting with " << z->nodeName << "...\n";
		effective_postorder_edge_iterator i(z, validFunctor);
		effective_postorder_edge_iterator e;
		for (; i != e; ++i)
#	else
	TreeNode *nd = (preorder ? GetFirstPreorder() : GetLastPreorder());
	for (; nd != NULL; nd = (preorder ? nd->GetNextPreorder() : nd->GetNextPostorder()))
#	endif
		{
#		if defined(TESTING_TOWARD_NODE_ITERATOR)
			TreeNode * nd = i->first;
#		endif
		if (verbosity > 2)
			{
			s << "\n";
			nd->AppendNodeInfo(s);
			}
		else
			{
			if (verbosity < 2)
				s << nd->briefDebugReport(verbosity);
			else    // verbosity == 2
				s << nd->oneLineDebugReport();

			// Output an arrow if not the last node in the series
			bool show_arrow = false;
			if (preorder)
				show_arrow = (nd != lastPreorder);
			else
				show_arrow = (nd != firstPreorder);
			if (show_arrow)
				s << " -> ";
			}
		}
#	if defined(TESTING_TOWARD_NODE_ITERATOR)
	}
#	endif

	return s;
	}

static int nodeIndex(TreeNode *p)
	{
	return (p == NULL) ? -1 : p->GetNodeNumber();
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Recursive component of listTree (below).
*/
void Tree::recDebugListTree(TreeNode *p, int nindent)
	{
	TreeNode *q;

	if (p != NULL)
		{
		unsigned nodeNum = p->GetNodeNumber();
		const char * nodeName = p->GetNodeName().c_str();
		if (p->NoChildren())
			fprintf(stderr, "%*cNode %d \"%s\" ", nindent, ' ', nodeNum, nodeName);
		else
			fprintf(stderr, "%*cNode %d ", nindent, ' ', nodeNum);

		fprintf(stderr, "(");
		fprintf(stderr, "a=%d l=%d r=%d", nodeIndex(p->par), nodeIndex(p->lChild), nodeIndex(p->rSib));
#		if SHOW_TRAVSEQ
			fprintf(stderr, " next=%d prev=%d", nodeIndex(p->nextPreorder), nodeIndex(p->prevPreorder));
#		endif

		fprintf(stderr, "): %g\n", p->edgeLen);
//		if (abortRequested())
//			return;

		q = p->lChild;
		while (q != NULL)
			{
			recDebugListTree(q, nindent + 1);
//			if (abortRequested())
//				return;
			q = q->rSib;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Shows tree topology and pointers (for debugging).
*/
void Tree::debugListTree()
	{
	recDebugListTree(GetRoot(), 2);
//	fprintf(stderr, "  (specialNode=%d", nodeDex(tree->specialNode));
//#	if SHOW_TRAVSEQ
//		fprintf(stderr, " lastPreorder=%d", nodeDex(tree->lastPreorder));
//#	endif
//	myputs(")\n");
	}

void Tree::debugMode(bool turn_on)
	{
	//Tree::gDebugOutput = turn_on;
	debugOutput = turn_on;
	//std::cerr << "@@@@@ debugMode(" << (turn_on ? "True" : "False") << ") @@@@@" << std::endl;
	}

bool Tree::DebugCheckTree(bool allowDegTwo, bool checkDataPointers, int verbosity) const
	{
	/**temp*/
	//if (!Tree::gDebugOutput)
	//	return true;
    if (!debugOutput)
        return true;
	const TreeNode * root = firstPreorder;
	PHYCAS_ASSERT(root);
	PHYCAS_ASSERT(root->par == NULL);
	PHYCAS_ASSERT(root->lChild != NULL);
	const TreeNode * currNd = root;
	const TreeNode * prevNd = NULL;
	const bool shouldHaveData = root->tipData != NULL;
	std::stack<const TreeNode *> next_pre_order_stack;
	std::stack<int> level_stack;
	int curr_level = 0;
	int next_level = 0;
	unsigned countedNNodes = 0;
	unsigned countedNLeaves = 1;

#   if !defined(NDEBUG)
	const unsigned expectedNLeaves = nTips; // +tipStorage.size()
	const unsigned expectedNNodes = expectedNLeaves + nInternals; // + internalNodeStorage.size()
#   endif
	std::string indent;
	while (currNd)
		{
		if (verbosity)
			{
			for (int i = 0; i < curr_level; ++i)
				std::cerr << "  ";
			if (verbosity == 1)
				std::cerr << currNd->briefDebugReport();
			else if (verbosity == 2)
				std::cerr << currNd->oneLineDebugReport();
			else
				{
				std::string s;
				currNd->AppendNodeInfo(s);
				std::cerr << s;
				}
			std::cerr << std::endl;
			}
		const TreeNode *nextNd = NULL;
		countedNNodes++;
		if (currNd->lChild)
			{
			PHYCAS_ASSERT(currNd->lChild->par == currNd);
			nextNd = currNd->lChild;
			next_level = curr_level + 1;
			if (!preorderDirty)
				{
				PHYCAS_ASSERT(currNd->nextPreorder == currNd->lChild);
				}
			if (currNd->par)
				{
				if ((!allowDegTwo) && currNd->par)
					{
					PHYCAS_ASSERT(currNd->lChild->rSib);
					}
				if (checkDataPointers)
					{
					if (shouldHaveData)
						{
						PHYCAS_ASSERT(currNd->internalData);
						}
					else
						{
						PHYCAS_ASSERT(currNd->internalData == NULL);
						}
					}
				}
			}
		else
			{
			countedNLeaves++;
			if (currNd->rSib == NULL)
				{
				if (next_pre_order_stack.empty())
					{
					nextNd = NULL;
					if (!preorderDirty)
                        {
						PHYCAS_ASSERT(currNd->nextPreorder == NULL);
                        }
					}
				else
					{
					nextNd = next_pre_order_stack.top();
					next_level = level_stack.top();
					if (!preorderDirty)
                        {
						PHYCAS_ASSERT(currNd->nextPreorder == next_pre_order_stack.top());
                        }
					next_pre_order_stack.pop();
					level_stack.pop();
					}
				}
			if (checkDataPointers)
				{
				if (shouldHaveData)
					{
					PHYCAS_ASSERT(currNd->tipData);
					}
				else
					{
					PHYCAS_ASSERT(currNd->tipData == NULL);
					}
				}
			}
		if (currNd->rSib)
			{
			PHYCAS_ASSERT(currNd->rSib->par == currNd->par);
			if (currNd->lChild)
				{
				next_pre_order_stack.push(currNd->rSib);
				level_stack.push(curr_level);
				}
			else
				{
				nextNd = currNd->rSib;
				next_level = curr_level;
				if (!preorderDirty)
                    {
					PHYCAS_ASSERT(currNd->rSib == currNd->nextPreorder);
                    }
				}
			}
#if 0
		PHYCAS_ASSERT(&*currNd->tree == this);
#endif
		if (!preorderDirty)
			{
			PHYCAS_ASSERT(currNd->prevPreorder == prevNd);
			}
		prevNd = currNd;
		currNd = nextNd;
		curr_level = next_level;
		if (nodeCountsValid)
			{
			PHYCAS_ASSERT(countedNNodes <= nTips + nInternals);
			}
		}
	PHYCAS_ASSERT(level_stack.empty());
	if (nodeCountsValid)
		{
		PHYCAS_ASSERT(countedNLeaves == expectedNLeaves);
		PHYCAS_ASSERT(countedNNodes == expectedNNodes);
		}
	if (!preorderDirty)
		{
		PHYCAS_ASSERT(lastPreorder == prevNd);
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles each element of which is one edge length. Edges are traversed in preorder fashion,
|	starting with the leftmost child of the root node (root node itself has no edge length).
*/
std::vector<double> Tree::EdgeLens()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	if (preorderDirty)
		RefreshPreorder();

	unsigned nnodes = GetNNodes();
	std::vector<double> v(nnodes - 1, 0.0);

	TreeNode * nd = GetFirstPreorder();
	unsigned i = 0;
	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		double brlen = nd->GetEdgeLen();
		v[i++] = brlen;
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles each element of which is one internal edge length. Edges are traversed in preorder
|   fashion, starting with the leftmost child of the root node (root node itself has no edge length).
*/
std::vector<double> Tree::GetInternalEdgeLens()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	if (preorderDirty)
		RefreshPreorder();

	unsigned nnodes = GetNInternals();
	std::vector<double> v(nnodes, 0.0);

	TreeNode * nd = GetFirstPreorder();
	unsigned i = 0;
	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
        if (nd->IsInternal() && !nd->IsSubroot())
            {
            double brlen = nd->GetEdgeLen();
            v[i++] = brlen;
            }
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles each element of which is one internal edge length. Edges are traversed in preorder
|   fashion, starting with the leftmost child of the root node (root node itself has no edge length).
*/
std::vector<double> Tree::GetExternalEdgeLens()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	if (preorderDirty)
		RefreshPreorder();

	unsigned nnodes = GetNTips();
	std::vector<double> v(nnodes - 1, 0.0);

	TreeNode * nd = GetFirstPreorder();
	unsigned i = 0;
	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
        if (nd->IsTip() || nd->IsSubroot())
            {
            double brlen = nd->GetEdgeLen();
            v[i++] = brlen;
            }
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Perform a preorder traversal of the tree and replace the existing edge lengths with those in the supplied vector
|	`new_edgelens'.
*/
void Tree::replaceEdgeLens(
  std::vector<double> new_edgelens)	/**< vector of edge lengths to be used to replace existing edge lengths in preorder traversal */
	{
	PHYCAS_ASSERT(GetNNodes() == new_edgelens.size() + 1);
	unsigned i = 0;
	for (preorder_iterator nd = begin(); nd != end(); ++nd)
		{
		if (!nd->IsAnyRoot())
			{
			double v = new_edgelens[i++];
			//std::cerr << boost::str(boost::format("@@@@@@@@@@ %.5f %.5f %s") % nd->GetEdgeLen() % v % nd->GetNodeName()) << std::endl;
			nd->SetEdgeLen(v);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing a key showing the bipartitions associated with each edge in the tree. Edges are
|	traversed in preorder fashion, starting with the leftmost child of the root node (root node itself is skipped).
*/
std::string Tree::KeyToEdges()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	//unsigned nnodes = GetNNodes();
	unsigned ntaxa = GetNObservables();
	RecalcAllSplits(ntaxa);
	if (preorderDirty)
		RefreshPreorder();

	std::string key_to_edges = "edge\tbipartition\n";

	TreeNode * nd = GetFirstPreorder();
	unsigned i = 1;
	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder(), ++i)
		{
		const Split & s = nd->GetSplitConst();
		key_to_edges += boost::str(boost::format("%d\t%s\n") % i % s.CreatePatternRepresentation());
		}

	return key_to_edges;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the sum of all edge lengths in the tree.  Uses preorder pointers to walk through tree, calling
|	RefreshPreorder() if these pointers have never been set or have been invalidated.
*/
double Tree::EdgeLenSum()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	//std::cerr << "\nEntering EdgeLenSum()" << std::endl; //POL-debug

	// todo: use algorithm accumulate for this, which entails writing node iterator
	//       this would make this function inlineable, I think
	if (preorderDirty)
		RefreshPreorder();

	TreeNode *nd = GetFirstPreorder();
	double sum = 0.0;

	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		//std::cerr << "nd->GetNodeNumber() = " << nd->GetNodeNumber() << ", sum = " << sum << std::endl; //POL-debug
		sum += nd->GetEdgeLen();
		}
	//std::cerr << "sum = " << sum << std::endl; //POL-debug

	//std::cerr << "\nLeaving EdgeLenSum()" << std::endl; //POL-debug
	return sum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of all internal edge lengths in the tree. Uses preorder pointers to walk through tree, calling
|	RefreshPreorder() if these pointers have never been set or have been invalidated.
*/
double Tree::internalEdgeLenSum()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	if (preorderDirty)
		RefreshPreorder();

	double sum = 0.0;
	TreeNode *nd = GetFirstPreorder();
	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsInternal())
			sum += nd->GetEdgeLen();
		}
	return sum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of all internal edge lengths in the tree. Uses preorder pointers to walk through tree, calling
|	RefreshPreorder() if these pointers have never been set or have been invalidated.
*/
double Tree::externalEdgeLenSum()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	if (preorderDirty)
		RefreshPreorder();

	double sum = 0.0;
	TreeNode *nd = GetFirstPreorder();
	for (nd = nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			sum += nd->GetEdgeLen();
		}
	return sum;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Sets `nTips' and `nInternals' to their correct values. Uses preorder pointers to walk through tree, calling
|	RefreshPreorder() if these pointers have never been set or have been invalidated.
*/
void Tree::RefreshNodeCounts()
	{
	nTips		= 0;
	nInternals	= 0;
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			++nTips;
		else
			++nInternals;
		}

	nodeCountsValid = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a tree from a string containing a "Newick" style (nested parentheses) tree description. If `newick' is
|	not a valid tree description, throws an XPhylogeny exception. Assumes that if any tip node names in the tree
|	description are convertible to integers, then tip node numbers should be set according to the converted node names.
|	Thus, it is considered an error to mix names that are not convertible to integers with names that are convertible,
|	and an XPhylogeny exception will be thrown if such a mixture is detected, or if a node number ends up being
|	duplicated or out of range.
*/
void Tree::BuildFromString(
  const std::string & newick,   /**< is the tree description */
  bool zero_based_tips)         /**< is true if tip node numbers used in the tree description start at zero */
	{
	PHYCAS_ASSERT(newick.length() > 0);

	// Assume tip node names represent the index of the taxon represented by the tip in the data matrix
	// and thus that tip node numbers should be set to the tip node name after converting it to an integer.
    // If zero_based_tips is false, then tip node numbers are set to the converted node name minus 1. This
	// saves having to call RectifyNumbers() for the common case of tree descriptions in which taxon names have
	// been replaced by their indices.
	//
	bool tip_numbers_equal_names = true;
	std::set<unsigned> numbers_used;

	// Stow away any existing nodes, reset nTips and nInternals to 0, and reset data members
	//
	Clear();

	// Save first tip node and reroot at this node before returning
	TreeNode *first_tip = NULL;

	try
		{
		// At the end, we will check to make sure that if *any* edge lengths were specified then
		// *all* edge lengths were specified
		unsigned nEdgeLengths = 0;

		// Create a root node
		TreeNode *nd = GetNewNode();
		firstPreorder = nd;

		// Some flags to keep track of what we did last
		enum {
			Prev_Tok_LParen		= 0x01,	/**< previous token was a left parenthesis ('(') */
			Prev_Tok_RParen		= 0x02,	/**< previous token was a right parenthesis (')') */
			Prev_Tok_Colon		= 0x04,	/**< previous token was a colon (':') */
			Prev_Tok_Comma		= 0x08,	/**< previous token was a comma (',') */
			Prev_Tok_Name		= 0x10,	/**< previous token was a node name (e.g. '2', 'P._articulata') */
			Prev_Tok_EdgeLen	= 0x20	/**< previous token was an edge length (e.g. '0.1', '1.7e-3') */
			};
		unsigned previous = Prev_Tok_LParen;

		// Some useful flag combinations
		unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
		unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
		unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
		unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
		unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

		// loop through the characters in newick, building up tree as we go
		std::string::const_iterator i = newick.begin();
		for (; i != newick.end(); ++i)
			{
			char ch = *i;
			if (iswspace(ch))
				continue;
			switch(ch)
				{
				case ';':
					//throw XPhylogeny("premature semicolon");
					break;

				case ')':
					// If nd is bottommost node, expecting left paren or semicolon, but not right paren
					if (nd->NoParent())
						throw XPhylogeny("too many right parentheses");
					// Expect right paren only after an edge length, a node name, or another right paren
					if (!(previous & RParen_Valid))
						throw XPhylogeny("unexpected right parenthesis");
					// Right paren means we are leaving a node for its parent, so set node number now if internal node
					if (nd->IsInternal() && nd->NumberNotYetAssigned())
						{
						nd->nodeNum = INTERNAL_NODENUM_OFFSET + nInternals++;
						}
					// Go down a level
					nd = nd->GetParent();
					if (nd->GetLeftChild()->GetRightSib() == NULL)
						throw XPhylogeny("internal node has only one child");
					previous = Prev_Tok_RParen;
					break;

				case ':':
					// Expect colon only after a node name or another right paren
					if (!(previous & Colon_Valid))
						throw XPhylogeny("unexpected colon");
					previous = Prev_Tok_Colon;
					break;

				case ',':
					// Expect comma only after an edge length, a node name, or a right paren
					if (nd->NoParent() || !(previous & Comma_Valid))
						throw XPhylogeny("unexpected comma");
					// Comma means we are leaving a node to create its sibling, so set node number now if internal
					if (nd->IsInternal() && nd->NumberNotYetAssigned())
						{
//						nd->nodeNum = nInternals++;
						nd->nodeNum = INTERNAL_NODENUM_OFFSET + nInternals++;
						}
					// Create the sibling
					nd->rSib = GetNewNode();
					nd->rSib->par = nd->par;
					nd = nd->rSib;
					previous = Prev_Tok_Comma;
					break;

				case '(':
					// Expect left paren only after a comma or another left paren
					if (!(previous & LParen_Valid))
						throw XPhylogeny("unexpected left parenthesis");
					// Create new_node above and to the left of the current node
					PHYCAS_ASSERT(nd->lChild == NULL);
					nd->lChild = GetNewNode();
					nd->lChild->par = nd;
					nd = nd->lChild;
					previous = Prev_Tok_LParen;
					break;

				case '\'':
					// Encountered an apostrophe, which always indicates the start of a
					// node name (but note that node names do not have to be quoted)

					// Expect node name only after a left paren (child's name), a comma (sib's name)
					// or a right paren (parent's name)
					if (!(previous & Name_Valid))
						{
						throw XPhylogeny("1> unexpected node name");
						}

					// Get the rest of the name
					nd->nodeName.clear();
					for (++i; i != newick.end(); ++i)
						{
						ch = *i;
						if (ch == '\'')
							break;
						else if (iswspace(ch))
							nd->nodeName << ' ';
						else
							nd->nodeName << ch;
						}
					if (ch != '\'')
						{
						workspace.clear();
						workspace << "expecting single quote to mark the end of node name (";
						//workspace << abbreviate(nd->GetNodeName(), 20);
						std::string s = nd->GetNodeName();
                        if (s.length() < 20)
                            workspace << s;
                        else
                            {
                            std::string ss(s.substr(0, 20));
                            ss << "...";
                            workspace << ss;
                            }
						workspace << ")";
						throw XPhylogeny(workspace);
						}

					if (nd->IsTip())
						{
						if (first_tip == NULL)
							first_tip = nd;
						if (tip_numbers_equal_names)
							{
							bool ok = SetNumberFromName(nd, numbers_used);	// may throw XPhylogeny
							if (!ok)
								{
								nd->nodeNum = nTips;
								tip_numbers_equal_names = false;
								}
							}
						else
							nd->nodeNum = nTips;
						++nTips;
						nd->SetObservable(true);
						}

					previous = Prev_Tok_Name;
					break;

				default:
					// Expecting either an edge length or an unquoted node name
					if (previous == Prev_Tok_Colon)
						{
						// Edge length expected (e.g. "235", "0.12345", "1.7e-3")
						workspace.clear();
						for (; i != newick.end(); ++i)
							{
							ch = *i;
							if (ch == ',' || ch == ')' || iswspace(ch))
								{
								--i;
								break;
								}

                            // handle embedded comments
                            if (ch == '[')
                                {
                                for (++i; i != newick.end(); ++i)
                                    {
                                    ch = *i;
                                    if (ch == '[')
                                        throw XPhylogeny("nested comments are not allowed in tree definitions");
                                    if (ch == ']')
                                        {
                                        ++i;
                                        ch = *i;
                                        break;
                                        }
                                    }
                                }

							bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
							if (!valid)
								throw XPhylogeny(boost::str(boost::format("invalid branch length character (%c)") % ch));
							workspace << ch;
							}
						double brlen = atof(workspace.c_str());
						nd->SetEdgeLen(brlen);

						++nEdgeLengths;
						previous = Prev_Tok_EdgeLen;
						}
					else
						{
						// Get the node name
						nd->nodeName.clear();
						for (; i != newick.end(); ++i)
							{
							ch = *i;
							if (ch == '(')
								throw XPhylogeny("unexpected left parenthesis inside node name");
							if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')')
								{
								--i;
								break;
								}
							nd->nodeName << ch;
							}

						// Expect node name only after a left paren (child's name), a comma (sib's name)
						// or a right paren (parent's name)
						if (!(previous & Name_Valid))
							{
							workspace.clear();
							workspace << "unexpected node name (";
							workspace << nd->nodeName;
							workspace << ")";
							throw XPhylogeny(workspace);
							}

						if (nd->IsTip())
							{
							if (first_tip == NULL)
								first_tip = nd;
							if (tip_numbers_equal_names)
								{
								bool ok = SetNumberFromName(nd, numbers_used);	// may throw XPhylogeny
								if (!ok)
									{
									nd->nodeNum = nTips;
									tip_numbers_equal_names = false;
									}
								}
							else
								nd->nodeNum = nTips;
							++nTips;
							nd->SetObservable(true);
							}

						previous = Prev_Tok_Name;
						}
				}
				if (i == newick.end())
					{
					break;
					}

			}

		if (firstPreorder->NumberNotYetAssigned())
			{
			firstPreorder->nodeNum = INTERNAL_NODENUM_OFFSET + nInternals++;
			}
		firstPreorder->SetEdgeLen(0.0);
		if (nEdgeLengths > 0)
			{
			if (nEdgeLengths < nTips + nInternals - 1)
				{
				std::string msg = boost::str(boost::format("some but not all edge lengths were specified: %d found, %d expected (= %d tips + %d internals - 1)") % nEdgeLengths % (nTips + nInternals - 1) % nTips % nInternals);
				throw XPhylogeny(msg);
				}
			hasEdgeLens = true;
			}
		else
			hasEdgeLens = false;
		if (!nd->NoParent())
			{
			throw XPhylogeny("too many left parentheses");
			}

        if (isRooted)
            {
            // Add a fake tip node as the left child of the curent root node and let first_tip point to it
            first_tip = GetNewNode();
            first_tip->par = firstPreorder;
            first_tip->rSib = firstPreorder->lChild;
            first_tip->lChild = NULL;
            #if 0//@TEMP DEL
            first_tip->nodeNum = UINT_MAX;
            #else
            first_tip->nodeNum = 4;
            #endif
            first_tip->nodeName = "root";
            first_tip->SetUnobservable();
            first_tip->SetEdgeLen(0.0);
            firstPreorder->lChild = first_tip;
            }

		RerootAtThisTip(first_tip);

	    if (tip_numbers_equal_names)
		    {
		    PHYCAS_ASSERT(!numbers_used.empty());
		    if (!zero_based_tips)
			    {
			    // Assuming that tip node numbers start at 1
			    for (preorder_iterator nd = begin(); nd != end(); ++nd)
				    {
				    if (nd->IsObservable())
					    {
                        if (nd->nodeNum == 0)
                            {
                            throw XPhylogeny("assuming that tip node numbers start at 1, but found a tip node number equal to 0");
                            }
					    nd->nodeNum = nd->nodeNum - 1;
					    }
				    }
			    }
		    }
		}
	catch(XPhylogeny x)
		{
		Clear();
		throw x;
		}

	// Renumber internal nodes so that they do not overlap with numbers of tip nodes
	//@POL should use for_each
	//std::cerr << "ADDING " << nTips << " to each internal node number" << std::endl; //POL temp
 	//for (preorder_iterator nd = begin(); nd != end(); ++nd)
 	//	{
 	//	if (nd->IsInternal())
 	//		{
 	//		//std::cerr << "RENUMBERING internal node " << nd->nodeNum << " (new number is " << (nd->nodeNum + nTips) << ")" << std::endl; //POL temp
 	//		nd->nodeNum = nd->nodeNum + nTips;
 	//		}
 	//	}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a parenthetical description of the tree (Newick format). If the tree has edge lengths, these are included in
|	the Newick description. Any nodes having non-empty node names are named in the Newick description. If a tip node
|	having node number x lacks a name, the "name" used for it in the tree description is x + 1.
*/
std::string & Tree::AppendNewick(
  std::string & s,                  /**< is the string to which the tree description will be appended */
  unsigned float_precision,         /**< is the number of decimal places to use in representing edge lengths */
  bool useNumbers,                  /**< if true, use node number (plus 1) as node name (only for tip nodes); if false, use node name if available, otherwise use node number for tips */
  bool showSupportAsNodeName,       /**< if true, support values will be used as node names for (internal nodes only) */
  bool showSupportAsEdgeLength)     /**< if true, edge lengths will equal support values for both tips (always 1.0) and internal nodes (used only in generating reference distributions for stepping-stone analyses when tree topology is variable) */
	{
	if (preorderDirty)
		RefreshPreorder();

	// Create a tool for formatting floating point values
	const unsigned floatWidth = UINT_MAX;	// make output width wide enough to hold the formatted value
	DoubleFormatter df(floatWidth, float_precision);

	bool showEdgeLens = HasEdgeLens() || showSupportAsEdgeLength;
	bool rooted = IsRooted();

	s << '(';
	unsigned openParens = 1;

	// Start with root node, which may actually represent an extant tip (unrooted trees are
	// rooted at one of the tips)
	TreeNode *nd = GetFirstPreorder();
	PHYCAS_ASSERT(nd->IsTipRoot());

	if (!rooted)
		{
        // In this case (unrooted tree), root node is actually a tip
        if (useNumbers || nd->nodeName.empty())
            {
            // output node number (plus 1) as the node name
            s << boost::str(boost::format("%d") % (nd->GetNodeNumber() + 1));
            if (useNumbers && nd->IsSelected())
                s << "[&!color=#-1241899]";
            }
        else
            {
            // output node name (enclosing in single quotes if name includes embedded spaces)
            if (nd->nodeName.find(" ") == std::string::npos)
			    s << nd->nodeName;
            else
			    s << "'" << nd->nodeName << "'";
            }
		}

	if (showEdgeLens && !rooted)
		{
        // In an unrooted tree, where the root node is actually an upside-down extant tip,
        // the root node's edge length is actually held by its only child
        s << ':';
        TreeNode *onlyChild = nd->GetLeftChild();
        PHYCAS_ASSERT(onlyChild != NULL);
        PHYCAS_ASSERT(onlyChild->rSib == NULL);
        if (showSupportAsEdgeLength)
            {
            double e = onlyChild->GetSupport();
            df.FormatDouble(s, e);
            }
        else
            {
            double e = onlyChild->GetEdgeLen();
            df.FormatDouble(s, e);
            }
		}

	nd = nd->GetNextPreorder();
	if (nd == NULL)
		{
		// Presumably this situation (a one-node tree) will never happen, but then it is
		// not good to presume anything
		return s << ')';
		}

    if (!rooted)
	    s << ',';

	if (nd->GetNextPreorder() == NULL)
		{
		// Two taxon tree, presumably this won't happen, but it is easy to deal with
		if (useNumbers || nd->nodeName.empty())
		    {
            // use node number (plus 1) as the node name
			s << boost::str(boost::format("%d") % (nd->GetNodeNumber() + 1));
            if (useNumbers && nd->IsSelected())
                s << "[&!color=#-1241899]";
			}
		else
            {
            // output node name (enclosing in single quotes if name includes embedded spaces)
            if (nd->nodeName.find(" ") == std::string::npos)
			    s << nd->nodeName;
            else
			    s << "'" << nd->nodeName << "'";
            }
		if (showEdgeLens)
			{
			s << ":";
			// 0.0 branch length because we gave the root's child's edge length
			// to the root because we print as a basal polytomy (not rooted at a leaf)
			df.FormatDouble(s, 0.0);
			}
		return s << ")";
		}

	// Note: no left parenthesis here; this is how we get the apparent polytomy at the root
	nd = nd->GetNextPreorder();

	// If we trip this we have a two taxon tree with 3 nodes (degree 2 node in the middle)
	PHYCAS_ASSERT(nd->GetNextPreorder());

	// Now visit all other nodes
	for (; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			// Output taxon names for tips if they are available, otherwise output tip node number plus 1
			if (useNumbers || nd->nodeName.empty())
			    {
				// use node number (plus 1) as the node name
    			s << boost::str(boost::format("%d") % (nd->GetNodeNumber() + 1));
                if (useNumbers && nd->IsSelected())
                    s << "[&!color=#-1241899]";
    			}
			else
                {
				// output node name (enclosing in single quotes if name includes embedded spaces)
                if (nd->nodeName.find(" ") == std::string::npos)
			        s << nd->nodeName;
                else
			        s << "'" << nd->nodeName << "'";
                }

			// Output edge length if they were provided
			if (showEdgeLens)
				{
				s << ':';
                if (showSupportAsEdgeLength)
                    {
                    // This is a tip node, so the support is always 1.0
                    double e = 1.0;
                    df.FormatDouble(s, e);
                    }
                else
                    {
                    double e = nd->GetEdgeLen();
                    df.FormatDouble(s, e);
                    }
				}

			if (nd->GetRightSib() != NULL)
				s << ',';
			else if (nd->nextPreorder != NULL)
				{
				// Go down until we find the left sib of the next preorder node,
				// outputting edge lengths as we go
				const TreeNode *tempAncNode = nd;
				do
					{
					s << ')';
					--openParens;
					PHYCAS_ASSERT(openParens > 0);
					tempAncNode = tempAncNode->par;
					PHYCAS_ASSERT(tempAncNode != NULL);

					if (showSupportAsNodeName && tempAncNode->GetSupport() > 0.0)
                        {
                        df.FormatDouble(s, tempAncNode->GetSupport());
                        }
                    else
                        {
                        // Output names for internal nodes if they are available, otherwise do not output anything
                        if (!(useNumbers || tempAncNode->nodeName.empty()))
                            {
                            // output node name (enclosing in single quotes if name includes embedded spaces)
                            if (tempAncNode->nodeName.find(" ") == std::string::npos)
                                s << tempAncNode->nodeName;
                            else
                                s << "'" << tempAncNode->nodeName << "'";
                            }
                        if (useNumbers && tempAncNode->IsSelected())
                            s << "[&!color=#-1241899]";
                        }

					if (showEdgeLens)
						{
						s << ":";
                        if (showSupportAsEdgeLength)
                            {
                            // This is an internal node
                            double e = tempAncNode->GetSupport();
                            df.FormatDouble(s, e);
                            }
                        else
                            {
                            double e = tempAncNode->GetEdgeLen();
                            df.FormatDouble(s, e);
                            }
						}
					}
				while (tempAncNode->rSib != nd->nextPreorder);
				s << ',';
				}
			}
		else
			{
			// nd is internal
			s << '(';
			++openParens;
			}
		}

	if (openParens > 0)
		{
		// We are at the upper right tip node of the tree now, so work back down
		// to the root, outputting branch lengths along the way
		const TreeNode *tmpNode = GetLastPreorder();
		while (openParens > 0)
			{
			s << ')';
			--openParens;
			PHYCAS_ASSERT(tmpNode != NULL);
			tmpNode = tmpNode->par;
			TreeNode *tmpNodePar = tmpNode->par;
			bool tmpNodeParIsRoot = tmpNodePar->IsTipRoot();
            if (showSupportAsNodeName && tmpNode->GetSupport() > 0.0)
                {
                df.FormatDouble(s, tmpNode->GetSupport());
                }
            else
                {
                if (!(useNumbers || tmpNode->nodeName.empty()))
                    {
                    // output node name (enclosing in single quotes if name includes embedded spaces)
                    if (tmpNode->nodeName.find(" ") == std::string::npos)
                        s << tmpNode->nodeName;
                    else
                        s << "'" << tmpNode->nodeName << "'";
                    }
                if (useNumbers && tmpNode->IsSelected())
                    s << "[&!color=#-1241899]";
                }
			if (showEdgeLens && !tmpNodeParIsRoot)
				{
				s << ":";
                if (showSupportAsEdgeLength)
                    {
                    // This is an internal node
                    double e = tmpNode->GetSupport();
                    df.FormatDouble(s, e);
                    }
                else
                    {
                    double e = tmpNode->GetEdgeLen();
                    df.FormatDouble(s, e);
                    }
				}
			}
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a postorder traversal, recalculating every split to reflect the nodes above and below the node being
|   visited.
*/
void Tree::RecalcAllSplits(
  unsigned max_nbits)   /**< is the maximum number of bits to allow in each split */
    {
	tree_id.clear();
    std::vector<unsigned> tips_seen;
    for (TreeNode * nd = GetLastPreorder(); nd != NULL; nd = nd->GetNextPostorder())
        {
        bool is_tip = nd->IsTip();
        Split & s = nd->GetSplit();
        if (s.GetNTaxa() != max_nbits)
            s.SetNTaxa(max_nbits);
        unsigned nn = nd->GetNodeNumber();
        tips_seen.push_back(nn);
        TreeNode * par = nd->GetParent();
        TreeNode * rsib = nd->GetRightSib();
        if (par)
            {
            Split & par_s = par->GetSplit();
            if (par_s.GetNTaxa() != max_nbits)
                par_s.SetNTaxa(max_nbits);
            if (is_tip && !rsib)
                {
                s.Reset();
                s.SetBit(nn);
                par_s.Reset();
                par_s.SetBit(nn);
                }
            else if (is_tip)
                {
                s.Reset();
                s.SetBit(nn);
                par_s.SetBit(nn);
                }
            else if (!rsib)
                {
                par_s.Reset();
                par_s |= s;
                }
            else
                {
                par_s |= s;
                }
            }
        else
            {
            PHYCAS_ASSERT(nd->IsTipRoot());
            if (nd->IsObservable())
                s.SetBit(nn);
            }
		tree_id.insert(s);
        }
    }

// below here lies previous contents of basic_tree.inl

/*----------------------------------------------------------------------------------------------------------------------
|	Attempts to locate a tip having the name `tipname'. If one can be found, its node number is returned. If no tip by
|   that name can be found, returns UINT_MAX.
*/
unsigned Tree::FindTipByName(
  std::string tipname)  /**< is the name to look for amongst tip nodes in the tree */
    {
    //std::cerr << boost::str(boost::format("*** Tree::FindTipByName ***\n*** tipname = \"%s\"") % tipname) << std::endl;
	for (preorder_iterator nd = this->begin(); nd != this->end(); ++nd)
		{
		if (nd->IsTip())
			{
            //std::cerr << boost::str(boost::format("*** %s <-> %s = %d") % tipname % nd->GetNodeName() % tipname.compare(nd->GetNodeName())) << std::endl;
            if (tipname.compare(nd->GetNodeName()) == 0)
                {
                //std::cerr << boost::str(boost::format("Found! node number is %d") % nd->GetNodeNumber()) << std::endl;
                return nd->GetNodeNumber();
                }
			}
        }
    return UINT_MAX;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Attempts to convert the name of `nd' to an integer (call it x) and then use x as the node number for `nd'. The
|	supplied set `used' is consulted to ensure that no two nodes are given the same number. If conversion of the name to
|	the integer value x is successful, and if x is not in the set `used', then: (1) the node number of `nd' is set to x;
|	(2) x is added to the set `used'; and (3) the function returns true. If the conversion of name to number fails
|	because name is not able to be converted to an integer, then one of two things can happen: (1) if `used' is empty,
|	then the function returns false, which serves as a signal that this function should not be called for future nodes
|	because apparently node names do not represent integers in this tree description; (2) if `used' is not empty, then
|	an XPhylogeny exception is thrown because some node numbers have already been set according to node names and thus
|	were we to discontinue converting names to numbers, the tree would end up with a mixture of node numbers derived
|	from node names and node numbers assigned arbitrarily. If the conversion of name to number fails because the number
|	has already been used, then an XPhylogeny exception is thrown because this is clearly an invalid tree description.
|	This function sets the `numbers_from_names' data member to true if it succeeds, and leaves it false if an exception
|	is thrown or the function returns false. Note: this function has no way of knowing whether node names (=numbers) in
|	the tree description start at 0 or 1, so `used' should be consulted after the tree is built and if 0 is not in the
|	set, a sweep of the tree should be made to decrement each tip node number.
*/
bool Tree::SetNumberFromName(TreeNode * nd, std::set<unsigned> & used)
	{
	PHYCAS_ASSERT(nd != NULL);
	std::string name = nd->GetNodeName();
	boost::trim(name);
	unsigned x = 0;
	try
		{
		x = boost::lexical_cast<unsigned>(name);
		}
	catch(boost::bad_lexical_cast &)
		{
		// signal that the node name could not be converted to an integer value
		x = UINT_MAX;
		}

	// This data member is set to false initially for every node examined.
	// It will thus end up being true after a tree is built only if the
	// conversion from names to numbers succeeds for every node.
	numbers_from_names = false;

	if (x == UINT_MAX)
		{
		if (used.empty())
			{
			// First node examined, so return false to signal that this function should
			// not be called for future node number assignments
			return false;
			}
		else
			{
			// This is not the first node examined, and previously examined nodes did have
			// names that could be converted to numbers.
			std::string s = "node name (";
			s << name;
			s << ") could not be converted to an integer to be used as the node number";
			throw XPhylogeny(s);
			}
		}

	// It is an error if x has already been used for a different node
	std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
	if (insert_result.second)
		{
		// insertion was made, so x has NOT already been used
		nd->SetNodeNum(x);
		}
	else
		{
		// insertion was not made, so set already contained x
		std::string s = "node number (";
		s << x;
		s << ") used for more than one node";
		throw XPhylogeny(s);
		}

	numbers_from_names = true;
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used to find out where we are in Python doctest examples.
*/
void Tree::DebugHere(
  std::string s)	/**< is a string indicating the current location */
	{
	std::cerr << "~~~~~~~ " << s << " ~~~~~~~" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the data member `tree_id'.
*/
const TreeID & Tree::getTreeID() const
	{
	return tree_id;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls RecalcAllSplits to build tree id by adding the split associated with every internal edge to the data member
|	`tree_id'.
*/
void Tree::buildTreeID()
	{
	RecalcAllSplits(GetNTips());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates and returns the Robinson-Foulds distance between this tree and `other'.
*/
unsigned Tree::robinsonFoulds(TreeShPtr other)
	{
	buildTreeID();
	other->buildTreeID();
	const TreeID & other_tree_id = other->getTreeID();

	TreeID::const_iterator first1 = tree_id.begin();
	TreeID::const_iterator first2 = other_tree_id.begin();
	TreeID::const_iterator last1 = tree_id.end();
	TreeID::const_iterator last2 = other_tree_id.end();

	// The tree id is a std::set of Split objects. We can determine if a split is in one tree and not the other
	// by comparing the two split objects with the < operator. Suppose there are 5 splits in each of two trees
	// (tree A and tree B), but the trees are not identical and only splits 5 and 7 are in both trees:
	//
	// A  1 2 3 - 5 - 7 -
	// B  - - - 4 5 6 7 8
	//
	// Define dAB to be the number of splits in A that are absent in B. Thus, dAB = 3. Likewise, dBA equals
	// the number of splits in B that are absent in A (3). In this case, dAB = dBA because both trees are
	// fully resolved. If one of the two trees was less resolved than the other, then dAB would differ from
	// dBA. This function calculates dAB and dBA separately and adds them together to obtain the return value.
	// Thus, for the example above, the result returned would be 6 (splits 1, 2, 3, 4, 6 and 8 are in one
	// tree and not the other).

	// times
	// thru
	// loop   first1  first2
	// ---------------------
	//  1       1       4     <- advance first1, num_differences = 1
	//  2       2       4     <- advance first1, num_differences = 2
	//  3       3       4     <- advance first1, num_differences = 3
	//  4       5       4     <- advance first2, num_differences = 4
	//  5       5       5     <- advance both
	//  6       7       6     <- advance first2, num_differences = 5
	//  7       7       7     <- advance both
	//  8       -       8     <- stop because first1 == last1, num_differences = 6

	unsigned num_differences = 0;
	while (first1 != last1 && first2 != last2)
		{
		if (*first1 < *first2)
			{
			++num_differences;
			++first1;
			}
		else if (*first2 < *first1)
			{
			++num_differences;
			++first2;
			}
		else
			{
			++first1;
			++first2;
			}
		}

	// Take care of incrementing num_differences to account for the fact that one
	// split set extends further than the other (i.e. step 8 in the example above)
	if (first1 == last1 && first2 != last2)
		{
		// tree id 2 (i.e. other tree) still has some elements left
		num_differences += std::distance(first2, last2);
		}
	else if (first2 == last2 && first1 != last1)
		{
		// tree id 1 (i.e. this tree) still has some elements left
		num_differences += std::distance(first1, last1);
		}
	return num_differences;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Scales all edge lengths in the tree by multiplying them by the supplied value `scaling_factor'.
*/
void Tree::ScaleAllEdgeLens(
  double scaling_factor)    /**< is the value by which all edge lengths will be multiplied */
	{
	std::for_each(begin(), end(), boost::lambda::bind(&TreeNode::ScaleEdgeLen, boost::lambda::_1, scaling_factor));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns all edge lengths in the tree to the supplied value `v'.
*/
void Tree::SetAllEdgeLens(double v)
	{
	std::for_each(begin(), end(), boost::lambda::bind(&TreeNode::SetEdgeLen, boost::lambda::_1, v));
	}

#if 0//use setRootedness instead

/*----------------------------------------------------------------------------------------------------------------------
|	Converts the tree to a rooted tree (if it is not already a rooted tree) and sets `isRooted' data member to true.
*/
void Tree::setRooted()
	{
    // Right now, this function can only be called immediately after construction or immediately after a call to Clear()
    // Much work needs to be done before an existing tree can be safely converted from unrooted to rooted
    PHYCAS_ASSERT(firstPreorder == NULL);
    isRooted = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Converts the tree to an unrooted tree (if it is not already unrooted) and sets `isRooted' data member to false.
*/
void Tree::setUnrooted()
	{
    // Right now, this function can only be called immediately after construction or immediately after a call to Clear()
    // Much work needs to be done before an existing tree can be safely converted from rooted to unrooted
    PHYCAS_ASSERT(firstPreorder == NULL);
    isRooted = false;
    }

#else//new

void Tree::setRootedness(bool rooted)
	{
	PHYCAS_ASSERT(firstPreorder == NULL);		// should only come here for an uninitialized Tree
	isRooted = rooted;

	#if 0//
	if (rooted != isRooted)
		{
		isRooted = rooted;
        if (rooted)
            {
            // Add a fake tip node as the left child of the curent root node and let first_tip point to it
            first_tip = GetNewNode();
            first_tip->par = firstPreorder;
            first_tip->rSib = firstPreorder->lChild;
            first_tip->lChild = NULL;
            first_tip->nodeNum = nTips;
            first_tip->nodeName = "root";
            first_tip->SetUnobservable();
            first_tip->SetEdgeLen(0.0);
            firstPreorder->lChild = first_tip;
            }
        else
        	{
        	//need to deroot
        	}

		RerootAtThisTip(first_tip);
    #endif//

    }

#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Pulls a node out of storage, if `internalNodeStorage' is not empty; otherwise, allocates memory for a new TreeNode. If it
|	is known in advance how many nodes will be needed, Reserve() can be called to create all nodes needed, storing them
|	in `internalNodeStorage'.
*/
TreeNode * Tree::GetNewNode()
	{
	if (!internalNodeStorage.empty())
		return PopInternalNode();
	return AllocNewNode();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for a new TreeNode object, sets the new node's tree pointer to this Tree object, and returns a
|   TreeNode pointer to the new node.
*/
TreeNode * Tree::AllocNewNode()
	{
	TreeNode * nd = new TreeNode();
#if 0
#	error Should not be compiled
	nd->SetTreeShPtr(TreeShPtr(this));
#endif
	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for `n' TreeNode objects and stores pointers to these in the `internalNodeStorage' vector. If some nodes
|	are already in `internalNodeStorage', only allocates enough new nodes to make the length of `internalNodeStorage' equal to `n'.
|	If the length of `internalNodeStorage' already equals or exceeds `n', returns immediately without allocating any new nodes.
*/
void Tree::Reserve(
  unsigned n)	/**< is the desired length of the `internalNodeStorage' vector */
	{
	PHYCAS_ASSERT(n > 0);
	unsigned num_existing_nodes = (unsigned)internalNodeStorage.size();
	if (num_existing_nodes >= n)
		return;

	unsigned num_nodes_needed = n - num_existing_nodes;
	for (unsigned i = 0; i < num_nodes_needed; ++i)
		{
        TreeNode * nd = new TreeNode();
#if 0
        nd->SetTreeShPtr(TreeShPtr(this));
#endif
		internalNodeStorage.push(nd);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `firstPreorder' points to a node having no parent, no right sibling and only one child.
*/
bool Tree::RootValid() const
	{
	return (firstPreorder != NULL && firstPreorder->IsTipRoot());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `isRooted' data member.
*/
bool Tree::IsRooted() const
	{
	return isRooted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `hasEdgelens' data member.
*/
bool Tree::HasEdgeLens() const
	{
	return hasEdgeLens;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `preorderDirty' data member.
*/
bool Tree::PreorderDirty() const
	{
	return preorderDirty;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `numbers_from_names' data member. Set member function Tree::SetNumberFromName for
|	details.
*/
bool Tree::TipNumbersSetUsingNames() const
	{
	return numbers_from_names;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of nodes currently composing the tree (i.e. `nTips' plus `nInternals'). For unrooted
|	trees of n taxa, one of the observable nodes serves as the root node, so this number can range from n+1 (star tree)
|	to 2n-2 (fully-resolved). For rooted trees, the range is n+2 (star tree) to 2n-1 (fully-resolved) because the root
|	node exists but is not identical to one of the observahble nodes. If `firstPreorder' is NULL, returns 0 because this
|	means that the tree does not yet have any nodes. For trees that have at least one node, calls RefreshNodeCounts() if
|	`nodeCountsValid' is false.
*/
unsigned Tree::GetNNodes()
	{
	// First check to make sure that a tree really exists
	if (firstPreorder == NULL)
		return 0;

	if (!nodeCountsValid)
		RefreshNodeCounts();
	return (nTips + nInternals);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	return
*/
unsigned Tree::GetNInternalsAllocated()
	{
	RefreshNodeCounts(); //temp
	return (unsigned)(internalNodeStorage.size() + GetNInternals());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of node objects stored in the data member `internalNodeStorage'
*/
unsigned Tree::NumInternalNodesStored()
	{
	return (unsigned)internalNodeStorage.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of node objects stored in the data member `tipStorage'
*/
unsigned Tree::NumTipNodesStored()
	{
	return (unsigned)tipStorage.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of tip nodes currently composing the tree. The data member `nTips' stores the number of degree
|	one nodes in the tree (tip means degree one). Note that one of these degree one nodes is the root node. The number
|	of tips equals the number of taxa for unrooted trees. For rooted trees of n taxa, there are n+1 tips because the
|	root node is a tip node but not an observable node. If `firstPreorder' is NULL, returns 0 because this means that
|	the tree does not yet have any nodes. For trees that have at least one node, calls RefreshNodeCounts() if
|	`nodeCountsValid' is false.
*/
unsigned Tree::GetNTips()
	{
	// First check to make sure that a tree really exists
	if (firstPreorder == NULL)
		return 0;

	if (!nodeCountsValid)
		RefreshNodeCounts();
	return nTips;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of observable tip nodes currently composing the tree. The data member `nTips' stores the number
|	of degree one nodes in the tree (tip means degree one). Note that one of these degree one nodes is the root node.
|	The number of observables equals the number of tips for unrooted trees. For rooted trees of n taxa, however, the
|	number of observables equals one fewer than the number of tips because the root node is a tip but not an observable.
|	If `firstPreorder' is NULL, returns 0 because this means that the tree does not yet have any nodes. For trees that
|	have at least one node, calls RefreshNodeCounts() if `nodeCountsValid' is false.
*/
unsigned Tree::GetNObservables()
	{
	return GetNTips() - (IsRooted() ? 1 : 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of internal nodes currently composing the tree (i.e. `nInternals'). The data member `nInternals'
|	stores the number of nodes in the tree of degree 2 or higher. Its value is 1 for a (rooted or unrooted) star tree,
|	n-2 for a fully-resolved unrooted tree, and n-1 for a fully-resolved rooted tree. If `firstPreorder' is NULL,
|	returns 0 because this means that the tree does not yet have any nodes. For trees that have at least one node, calls
|	RefreshNodeCounts() if `nodeCountsValid' is false.
*/
unsigned Tree::GetNInternals()
	{
	// First check to make sure that a tree really exists
	if (firstPreorder == NULL)
		return 0;

	if (!nodeCountsValid)
		RefreshNodeCounts();
	return nInternals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns first node in the preorder sequence, which is equivalent to the tip root in an unrooted tree.
*/
TreeNode * Tree::GetRoot()
	{
	return GetFirstPreorder();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the subroot node, which is the only child of the root tip node in an unrooted tree.
*/
TreeNode * Tree::GetSubroot()
	{
	TreeNode * root_tip = GetFirstPreorder();
	PHYCAS_ASSERT(root_tip);
	TreeNode * subroot  = root_tip->GetLeftChild();
	PHYCAS_ASSERT(subroot);
	PHYCAS_ASSERT(subroot->GetNextPreorder() == subroot->GetLeftChild());
	return subroot;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return value of data member `firstPreorder'. Note: makes a (rather expensive) call to RefreshPreorder() if the
|	preorder pointers need to be updated (i.e., `preorderDirty' is true).
*/
TreeNode * Tree::GetFirstPreorder()
	{
	if (preorderDirty)
		RefreshPreorder();
	return firstPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return value of data member `lastPreorder'. Note: makes a (rather expensive) call to RefreshPreorder() if the
|	preorder pointers need to be updated (i.e., `preorderDirty' is true).
*/
TreeNode * Tree::GetLastPreorder()
	{
	if (preorderDirty)
		RefreshPreorder();

	return lastPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return value of data member `firstPreorder'. Note: makes a (rather expensive) call to RefreshPreorder() if the
|	preorder pointers need to be updated (i.e., `preorderDirty' is true).
*/
const TreeNode * Tree::GetFirstPreorderConst() const
	{
	if (preorderDirty)
		RefreshPreorder();
	return firstPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return value of data member `lastPreorder'. Note: makes a (rather expensive) call to RefreshPreorder() if the
|	preorder pointers need to be updated (i.e., `preorderDirty' is true).
*/
const TreeNode * Tree::GetLastPreorderConst() const
	{
	if (preorderDirty)
		RefreshPreorder();

	return lastPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a vector of pointers to all TreeNode objects that are responsible for edge lengths. The nodes
|	that qualify are all nodes in the tree except the root node, which is a tip node whose edge is managed by the
|	internal node that represents its only child.
*/
Tree::TreeNodeVec Tree::GetNodesWithEdges()
	{
	TreeNodeVec v;

	TreeNode * nd = GetFirstPreorder();	// this is the root node, which we skip
	for (nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		v.push_back(nd);
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return preorder_iterator representing beginning of the preorder sequence of nodes. Note: this function calls
|	GetFirstPreorder(), which will make the (expensive) call to RefreshPreorder() if preorder pointers need to be
|	updated.
*/
preorder_iterator Tree::begin()
	{
	return preorder_iterator(GetFirstPreorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return preorder_iterator representing end of the preorder sequence of nodes.
*/
preorder_iterator Tree::end()
	{
	return preorder_iterator();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a parenthetical description of the tree (Newick format). If the tree has edge lengths, these are included in
|	the Newick description. Any nodes having non-empty node names are named in the Newick description. Tip nodes that
|	do not have stored names are given a name that equals 1 plus their node number. Calls AppendNewick() to do the work,
|	but unlike AppendNewick() does not have a std::string reference parameter and returns a std::string by value rather
|	than by reference.
*/
std::string Tree::MakeNewick(
  unsigned float_precision)
	{
	std::string s;
	AppendNewick(s, float_precision, false, true);   //POLNEW: added false, true
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a parenthetical description of the tree (Newick format). If the tree has edge lengths, these are included in
|	the Newick description. Any nodes having non-empty node names are named in the Newick description. Tip nodes that
|	do not have stored names are given a name that equals 1 plus their node number. Calls AppendNewick() to do the work,
|	but unlike AppendNewick() does not have a std::string reference parameter and returns a std::string by value rather
|	than by reference.
*/
std::string Tree::MakeNumberedNewick(
  unsigned float_precision)
	{
	std::string s;
	AppendNewick(s, float_precision, true, true);   //POLNEW: added last true
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a parenthetical description of the tree (Newick format) for use in a stepping-stone variable-topology
|   reference distribution. This is an unsual tree description in that support values are substituted for edge lengths.
*/
std::string Tree::MakeNewickForRefDist(
  unsigned float_precision)
	{
	std::string s;
	AppendNewick(s, float_precision, true, false, true);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `treeid_valid' to false, indicating that this tree's TreeID is no longer valid because of a recent
|	change in the tree topology.
*/
void Tree::InvalidateTreeID()
	{
	treeid_valid = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Flags the `tree_id' data member as being in need of recalculation. Lazy recalculation is used because sometimes
|	several operations in quick succession all make changes to the topology of the tree without needing to reference
|	the `tree_id' data member. Recalculating `tree_id' for each of these would slow down the program unecessarily.
*/
void Tree::InvalidateID()
	{
	treeid_valid = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `nodeCountsValid' to false, indicating that the node counts reported by member functions such as
|	GetNNodes(), GetNTips(), GetNInternals(), etc., may no longer be valid due to a recent change in the tree topology.
*/
void Tree::InvalidateNodeCounts()
	{
	nodeCountsValid = false;
	}

