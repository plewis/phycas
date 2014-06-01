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

#include <cfloat>
#include "basic_tree_node.hpp"
#include "phycas_string.hpp"
#include "basic_tree.hpp"

using namespace phycas;

const double TreeNode::edgeLenEpsilon = 1.e-12;
const double TreeNode::edgeLenDefault = 0.1;
const double TreeNode::edgeLenInitValue = DBL_MAX;
const unsigned TreeNode::nodeNumInitValue = UINT_MAX;

//static unsigned debug_treenode_number = 0;

/*----------------------------------------------------------------------------------------------------------------------
|	Begins with left child of parent of `this' and calls GetRightSib() until the left sibling of `this' is located.
*/
TreeNode * TreeNode::FindLeftSib()
	{
	TreeNode * nd = this->GetParent();

	// If this has no parent, then there's no way it can have a left sibling
	if (nd == NULL)
		return NULL;

	nd = nd->GetLeftChild();

	// Parent of this should have at least one child
	PHYCAS_ASSERT(nd != NULL);

	// If left child of this's parent equals this, then this is an only child and has no siblings
	if (nd == this)
		return NULL;

	TreeNode * leftsib = NULL;
	while (nd != this)
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
TreeNode * TreeNode::FindRightmostChild()
	{
	TreeNode * curr = this->GetLeftChild();
	TreeNode * rightmost = NULL;
	while (curr != NULL)
		{
		rightmost = curr;
		curr = curr->GetRightSib();
		}
	return rightmost;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Zig-zags up the right side of the clade starting with the node `start' until reaching a tip. This tip node is the
|	last node in the clade according to the preorder sequence. Zig-zagging is necessary because there is no pointer to
|	the right-most child of a node, so we must use FindRightmostChild to move up the clade. Assumes `start' is non-NULL.
*/
TreeNode * TreeNode::FindLastPreorderInClade()
	{
	TreeNode * curr = this;
	TreeNode * rChild = curr->FindRightmostChild();
	while (rChild != NULL)
		{
		curr = rChild;
		rChild = curr->FindRightmostChild();
		}
	return curr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a vector of pointers to the children of this node.
*/
std::vector<TreeNode *> TreeNode::GetChildren() const
    {
    std::vector<TreeNode *> v;
    TreeNode * c = lChild;
    while (c)
        {
        v.push_back(c);
        c = c->rSib;
        }
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Only sets par pointer and either lChild (of parent) or rSib (of left sib).
*/
void TreeNode::AddChild(TreeNode *nd)
    {
    TreeNode * leftSib = FindRightmostChild();
    if (leftSib == 0L)
        lChild = nd;
    else
        leftSib->rSib = nd;
    nd->par = this;

    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all pointers to 0, nodeNum to TreeNode::nodeNumInitValue, and edgeLen to TreeNode::edgeLenInitValue.
*/
TreeNode::TreeNode()
  :nodeNum(TreeNode::nodeNumInitValue),
  edgeLen(TreeNode::edgeLenInitValue),
  lChild(0),
  par(0),
  rSib(0),
  nextPreorder(0),
  prevPreorder(0),
  adhocTreeNodePtr(0),
  //observable(false),
  support(0.0),
  tmp(0.0),
  x(0.0),
  y(0.0),
  selected(false),
  tipData(0),
  tipDataDeleter(0),
  internalData(0),
  internalDataDeleter(0),
  correspondingNd(0),
  ptr(NULL)
  	{
  	//std::cerr << "=====> creating TreeNode " << (++debug_treenode_number) << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor for TreeNode.
*/
TreeNode::~TreeNode()
	{
	//std::cerr << "In node destructor" << std::endl;
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all pointers to 0, nodeNum to TreeNode::nodeNumInitValue and edgeLen to TreeNode::edgeLenInitValue.
|	Also deletes any structures assigned to `tipData' or `internalData' using the callbacks provided when these
|	structures were allocated. Basically, returns node to its just-constructed state.
*/
void TreeNode::Clear()
	{
	//std::cerr << "tree use count = " << tree.use_count() << std::endl;
	//tree.reset();
	ResetInternalData();
	ResetTipData();

	lChild			= 0;
	par				= 0;
	rSib			= 0;
	nextPreorder	= 0;
	prevPreorder	= 0;
    adhocTreeNodePtr = 0;
	nodeNum			= TreeNode::nodeNumInitValue;
	edgeLen			= TreeNode::edgeLenInitValue;
	nodeName		= "";
	//observable		= false;
	support			= 0.0;
	tmp				= 0.0;
	x				= 0.0;
	y				= 0.0;
	selected		= false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the edge length (value of `edgeLen' data member multiplied by the tree's scaler value).
*/
double TreeNode::GetEdgeLen() const
	{
	return edgeLen;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the split associated with this node.
*/
const Split & TreeNode::GetSplitConst() const
	{
	return split;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the split associated with this node.
*/
Split & TreeNode::GetSplit()
	{
	return split;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `edgeLen' to the product of the current value of `edgeLen' and the supplied `scaling_factor'.
*/
void TreeNode::ScaleEdgeLen(
  double scaling_factor)			/**< is the value by which the edge length will be multiplied */
	{
    PHYCAS_ASSERT(scaling_factor > 0.0);
    double new_edgeLen = edgeLen*scaling_factor;
    SetEdgeLen(new_edgeLen);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Allows write access to protected data member `edgeLen'.
*/
void TreeNode::SetEdgeLen(
  double x)							/**< is the new edge length value */
	{
	edgeLen = (x < TreeNode::edgeLenEpsilon ? TreeNode::edgeLenEpsilon : x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of immediate descendants (should be 0 for a tip node, 1 for the root (tip) node, and 2 or more for
|	internal nodes.
*/
unsigned TreeNode::CountChildren() const
	{
	unsigned nDescendants = 0;
	for (const TreeNode * child = GetLeftChildConst(); child != NULL; child = child->GetRightSibConst())
		 ++nDescendants;
	return nDescendants;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Not yet documented.
*/
void TreeNode::AppendNodeInfo(std::string &s, bool num_and_name_only) const
	{
	std::string tmpstr;
	const std::string &nm = GetNodeName();

	if (nodeNum == TreeNode::nodeNumInitValue)
		{
		tmpstr << "<no number>";
		}
	else
		{
		tmpstr = str(boost::format("%d") % nodeNum);
		}

	if (num_and_name_only)
		{
		if (IsTipRoot())
			{
			s << "{";
			s << tmpstr;
			s << "}";
			}
		else if (IsTip())
			{
			s << "(";
			s << tmpstr;
			s << ")";
			}
		else if (IsInternal())
			{
			s << "[";
			s << tmpstr;
			s << "]";
			}
		else
			{
			s << "<***error***";
			s << tmpstr;
			s << "***error***>";
			}

		if (nm.length() > 0)
			{
			s << " \'";
			s << nm;
			s << "\'";
			}
		}
	else
		{
		// Node number
		s << "\nNode number: ";
		s << tmpstr;

		// Node status
		s << "\nNode status: ";
		if (IsTip())
			s << "tip";
		else if (IsTipRoot())
			s << "root";
		else if (IsInternal())
			s << "internal";
		else
			s << "unknown (this is a bug)";

		// Node name
		s << "\nNode name: ";
		if (nm.length() > 0)
			s << nm;
		else
			s << "<unnamed>";

		// lChild
		s << "\nLeft child: ";
		if (GetLeftChildConst() != NULL)
			lChild->AppendNodeInfo(s, true);
		else
			s << "<no left child>";

		// rSib
		s << "\nRight sibling: ";
		if (GetRightSibConst() != NULL)
			rSib->AppendNodeInfo(s, true);
		else
			s << "<no right sib>";

		// par
		s << "\nParent: ";
		if (GetParentConst() != NULL)
			par->AppendNodeInfo(s, true);
		else
			s << "<no parent>";

		// nextPreorder
		s << "\nNext preordert: ";
		if (GetNextPreorderConst() != NULL)
			nextPreorder->AppendNodeInfo(s, true);
		else
			s << "<no next preorder>";

		// prevPreorder
		s << "\nNext postorder: ";
		if (GetNextPostorderConst() != NULL)
			prevPreorder->AppendNodeInfo(s, true);
		else
			s << "<no next postorder>";

		// Edge length
		s << "\nEdge length: ";
		if (edgeLen != edgeLenInitValue)
			{
			tmpstr = str(boost::format("%f") % GetEdgeLen());
			s << tmpstr;
			}
		else
			s << "<no edge length>";
		}
	}

std::string	TreeNode::briefDebugReport(
  unsigned verbosity) const
	{
	std::string nm = GetNodeName();
	std::string tmpstr;
    unsigned namelen = (unsigned)nm.length();

    if (verbosity == 0)
        {
        if (namelen == 0)
            {
		    if (IsTip())
                tmpstr = str(boost::format("(%d)") % GetNodeNumber());
            else
                tmpstr = str(boost::format("[%d]") % GetNodeNumber());
            }
        }
    else    // verbosity > 0
        {
        //     a (0) -> ? [5] -> b (1) -> ? [7] -> c (2) -> ? [6] -> d (3) -> e (4)
        if (namelen == 0)
            nm = "?";
		if (IsTip())
            tmpstr = str(boost::format(" (%d)") % GetNodeNumber());
        else
            tmpstr = str(boost::format(" [%d]") % GetNodeNumber());
        }

	return str(boost::format("%s%s") % nm % tmpstr);
	}

std::string	TreeNode::oneLineDebugReport() const
	{
	const std::string self_str = this->briefDebugReport();
	std::string lch_str;
	if (lChild)
		lch_str = lChild->briefDebugReport();
	else
		lch_str = "NULL";
	std::string par_str;
	if (par)
		par_str = par->briefDebugReport();
	else
		par_str = "NULL";
	std::string sib_str;
	if (rSib)
		sib_str = rSib->briefDebugReport();
	else
		sib_str = "NULL";
	std::string n_str;
	if (nextPreorder)
		n_str = nextPreorder->briefDebugReport();
	else
		n_str = "NULL";
	std::string prev_str;
	if (prevPreorder)
		prev_str = prevPreorder->briefDebugReport();
	else
		prev_str = "NULL";
	return str(boost::format("%s | lchild=%s | par=%s | rsib=%s | next=%s | prev=%s") % self_str % lch_str % par_str % sib_str % n_str % prev_str);
	}

//////////////////////////// below here previously in basic_tree_node.inl /////////////////////////////////
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the rSib or the parent's left child.
|	Returns NULL if the parent has only one child.
*/
TreeNode * TreeNode::FindNextSib()
	{
	if (rSib != NULL)
		return rSib;
	PHYCAS_ASSERT(par);
	if (par->lChild == this)
		return NULL;
	return par->lChild;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' share pointer data member to the supplied value `t'. This allows subsequent access to the tree of
|   which this node is a part in order to, for example, gain access to the whole-tree scaling factor.
*/
void TreeNode::SetTreeShPtr(
  TreeShPtr t) /**> is a pointer to the tree (should be called just after a node is created) */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `support'.
*/
float TreeNode::GetSupport() const
	{
	return support;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `x'.
*/
float TreeNode::GetX()
	{
	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `observable' to true.
*/
void TreeNode::SetObservable()
	{
	// observable has been deprecated
//	observable = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `observable' to false.
*/
void TreeNode::SetUnobservable()
	{
	// observable has been deprecated
//	observable = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `support' to the value `x'.
*/
void TreeNode::SetSupport(
  float x)						/**< is the new support value */
	{
	support = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `x' to the value `xx'.
*/
void TreeNode::SetX(
  float xx)						/**< is the new x-coordinate */
	{
	x = xx;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the data member `y'.
*/
float TreeNode::GetY()
	{
	return y;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `y' to the value `yy'.
*/
void TreeNode::SetY(
  float yy)						/**< is the new y-coordinate */
	{
	y = yy;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets `nodeNum' data member to the supplied value.
*/
void TreeNode::SetNodeNum(
  unsigned num)						/**< is the new node number */
	{
	nodeNum = num;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets `nodeName' data member to the supplied value.
*/
void TreeNode::SetNodeName(
  std::string name)					/**< is the new node name */
	{
	nodeName = name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `lChild' data member is not zero.
*/
bool TreeNode::HasChildren() const
	{
	return (lChild != 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `lChild' data member equals zero.
*/
bool TreeNode::NoChildren() const
	{
	return (lChild == 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `par' data member is not zero.
*/
bool TreeNode::HasParent() const
	{
	return (par != 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `par' data member equals zero.
*/
bool TreeNode::NoParent() const
	{
	return (par == 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this node is an internal node that is the only child of the root tip node.
*/
bool TreeNode::IsSubroot() const
	{
	return IsInternal() && par->IsTipRoot();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if this node is an internal node that is the only child of the root tip node.
*/
bool TreeNode::IsExternalEdge() const
	{
	return IsTip() || par->IsTipRoot();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if either NoChildren() or IsTipRoot() returns true. A node is a tip node if it has degree
|	one, which is true if only one edge connects it to the rest of the tree. This is true for nodes that have no
|	children (their only connection to the rest of the tree is through their parent), but is also true for the root node
|	(the only node in a tree that has only one child).
*/
bool TreeNode::IsTip() const
	{
	return NoChildren() || IsTipRoot();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if the node can potentially have associated observations. For unrooted trees, this will be
|	true for any IsTip() node, but for rooted trees, it is only true for NoChildren() nodes because the root is
|	not observable even though it is degree one. Note that a TreeNode object has no way of knowing whether the tree to
|	which it is connected is rooted or unrooted, so this function simply returns the value of the `observable' data
|	member, and depends on the tree to set that data member correctly when the node is incorporated into the tree.
*/
bool TreeNode::IsObservable() const
	{
	#if 0//
	return observable;
	#else//
	return (IsTip() && (nodeNum != UINT_MAX));
	#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets flag that serves to indicate whether the node can potentially have associated observations. For unrooted trees,
|	this should be true for any IsTip() node, but for rooted trees, it should only be true for NoChildren() nodes
|	because the root is	not observable even though it is degree one. Note that a TreeNode object has no way of knowing
|	whether the tree to which it is connected is rooted or unrooted, so this function should be used by functions that
|	build trees to set the `observable' data member correctly. One use of this data member is in deciding which nodes
|	should be equipped with structures that store a copy of a particular row in the data matrix when preparing the tree
|	for likelihood calculations.
*/
void TreeNode::SetObservable(bool is_observable)
	{
//	observable = is_observable;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if IsTip() returns false.
*/
bool TreeNode::IsInternal() const
	{
	return !IsTip();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `nodeNum' equals the value `TreeNode::nodeNumInitValue', which is the value assigned to
|	the `nodeNum' data member of all newly-created nodes.
*/
bool TreeNode::NumberNotYetAssigned() const
	{
	return nodeNum == TreeNode::nodeNumInitValue;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `edgeLen' equals the value `TreeNode::edgeLenInitValue', which is the value assigned to
|	the `edgeLen' data member of all newly-created nodes.
*/
bool TreeNode::EdgeLenNotYetAssigned() const
	{
	return edgeLen == TreeNode::edgeLenInitValue;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this is a node with no parent and no siblings. This function does not pay attention to
|	the number of children. If you know that the root node should be a tip, use TreeNode::IsTipRoot instead (that
|	function also requires root to have just one child). If you know that the root node should be an internal node, use
|	TreeNode::IsInternalRoot instead (that function also requires that the root have at least two children).
*/
bool TreeNode::IsAnyRoot() const
	{
	bool parentless = !par;
	bool only_child = !rSib;
	return (parentless && only_child);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this is a node with no parent, no siblings and only one child. Use this function if
|	you specifically want to test for the case where a tip node is serving as the root of the tree. If an internal node
|	is serving as the root of the tree, it is expected to have more than one child, and the method IsInternalRoot()
|	should be used instead.
*/
bool TreeNode::IsTipRoot() const	//POL: formerly IsRoot()
	{
	bool parentless = !par;
	bool only_child = !rSib;
	bool one_child = (lChild && !(lChild->rSib));
	return (parentless && only_child && one_child);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if this is a node with no parent, no siblings and at least 2 children.
*/
bool TreeNode::IsInternalRoot() const
	{
	bool parentless = !par;
	bool only_child = !rSib;
	bool at_least_two_children = (lChild && lChild->rSib);
	return (parentless && only_child && at_least_two_children);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the edge length (value of `nodeNum' data member).
*/
unsigned TreeNode::GetNodeNumber() const
	{
	return nodeNum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the name of this node (value of `nodeName' data member).
*/
const std::string & TreeNode::GetNodeName() const
	{
	return nodeName;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `selected' data member. This function (as well as its cousings SelectNode and UnselectNode) may be
|	temporary. They are intended to be used for situations in which some nodes need to be selected, either in a
|	graphical context when manipulating trees via a GUI, or in a non-graphical context, such as identification of nodes
|	for lazy recalculation. Both of these applications are probably better served by having more specific flags,
|	however, so these functions will probably eventually be moved into other structures.
*/
bool TreeNode::IsSelected() const
	{
	return selected;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to true. (See notes under IsSelected.)
*/
void TreeNode::SelectNode()
	{
	selected = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to false. (See notes under IsSelected.)
*/
void TreeNode::UnselectNode()
	{
	selected = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `tipData' data member.
*/
TipData * TreeNode::GetTipData()
	{
	return tipData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `internalData' data member.
*/
InternalData * TreeNode::GetInternalData()
	{
	return internalData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `tipData' data member.
*/
const TipData * TreeNode::GetTipData() const
	{
	return tipData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `internalData' data member.
*/
const InternalData * TreeNode::GetInternalData() const
	{
	return internalData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `tipData' to its just-constructed state, namely 0. If `tipData' is non-zero, calls the associated
|	`tipDataDeleter' function to delete the object. Assumes that if `tipData' is non-zero, then `tipDataDeleter' is
|	valid.
*/
void TreeNode::ResetTipData()
	{
	if (tipData)
		{
		PHYCAS_ASSERT(tipDataDeleter);
		tipDataDeleter(tipData);
		}
	tipData = 0;
	tipDataDeleter = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `internalData' to its just-constructed state, namely 0. If `internalData' is non-zero, calls the
|	associated `internalDataDeleter' function to delete the object. Assumes that if `internalData' is non-zero, then
|	`internalDataDeleter' is valid.
*/
void TreeNode::ResetInternalData()
	{
	if (internalData)
		{
		PHYCAS_ASSERT(internalDataDeleter);
		internalDataDeleter(internalData);
		}
	internalData = 0;
	internalDataDeleter = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `tipData' to `d' and `tipDataDeleter' to `f'
*/
 void TreeNode::SetTipData(
  TipData * d,						/**< is a pointer to the TipData data structure to be assigned to `tipData' */
  TreeNode::TipDataDeleter f)		/**< is the function object that knows how to delete a TipData object */
	{
	ResetTipData();
	tipData = d;
	tipDataDeleter = f;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `internalData' to `c' and `internalDataDeleter' to `f'
*/
void TreeNode::SetInternalData(
  InternalData * c,						/**< is a pointer to the InternalData data structure to be assigned to `internalData' */
  TreeNode::InternalDataDeleter f)		/**< is the function object that knows how to delete a InternalData object */
	{
	ResetInternalData();
	internalData = c;
	internalDataDeleter = f;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `lChild'.
*/
TreeNode	* TreeNode::GetLeftChild()
	{
	return lChild;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `lChild'. Useful for iterating through children of a node inside a const
|	function.
*/
const TreeNode * TreeNode::GetLeftChildConst() const
	{
	return lChild;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `rSib'.
*/
TreeNode	* TreeNode::GetRightSib()
	{
	return rSib;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `rSib'. Useful for iterating through children of a node inside a const
|	function.
*/
const TreeNode * TreeNode::GetRightSibConst() const
	{
	return rSib;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `par'.
*/
TreeNode	* TreeNode::GetParent()
	{
	return par;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `par'.
*/
const TreeNode	* TreeNode::GetParentConst() const
	{
	return par;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `nextPreorder'.
*/
TreeNode	* TreeNode::GetNextPreorder()
	{
	return nextPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `prevPreorder'.
*/
TreeNode	* TreeNode::GetNextPostorder()
	{
	return prevPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `nextPreorder'.
*/
const TreeNode	* TreeNode::GetNextPreorderConst() const
	{
	return nextPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `prevPreorder'.
*/
const TreeNode	* TreeNode::GetNextPostorderConst() const
	{
	return prevPreorder;
	}
