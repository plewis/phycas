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

#ifndef PYPHY_BASIC_TREE_HPP
#define PYPHY_BASIC_TREE_HPP

#include <string>
#include <stack>
#include <vector>
#include <set>
#include <iostream>
#include <boost/algorithm/string.hpp>	// used by SetNumberFromName member function
#include <boost/lexical_cast.hpp>		// used by SetNumberFromName member function
#include "basic_tree_node.hpp"
#include "tree_iterators.hpp"
#include "xphylogeny.hpp"
#include "phycas_string.hpp"
#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this

//	To do next:
//		o ASCII tree drawing
//		o splits
//		o internal nodes should be numbered starting with nTips
//		o should ensure that tip node numbering sequence starts and 0 and ends at ntax-1 with no gaps
//			(see bottom of BuildFromString). (Plant to deal with valid exceptions to this rule later.)

namespace phycas{

#include "split.hpp"
typedef std::set<Split> TreeID;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the notion of a phylogenetic tree. This class has only methods necessary for representing the tree,
|	copying trees from other trees, manipulating the tree, and representing the tree graphically. It does not perform
|	any specialized activities, such as computing its own likelihood: these sorts of activities are left to member
|	functions of other classes.
*/
class Tree : public boost::enable_shared_from_this<Tree>
	{
	public:
		//static bool gDebugOutput;
		bool debugOutput;

		friend class TreeManip;
		friend class BushMove;
		friend class TreeLikelihood;

		typedef std::vector<TreeNode *> TreeNodeVec;
		typedef std::stack<TreeNode *> TreeNodeStack;

		// Constructors/destructor
		//
								Tree();
		virtual					~Tree();

		// Accessors
		//
		unsigned				GetNNodes();
		unsigned				GetNTips();
		unsigned				GetNInternals();
		unsigned 				GetNInternalsAllocated();
		unsigned				GetNObservables();
		TreeNode *				GetFirstPreorder();
		TreeNode *				GetRoot();
		TreeNode *				GetSubroot();
		TreeNode *				GetLastPreorder();
		const TreeNode *		GetFirstPreorderConst() const;
		const TreeNode *		GetLastPreorderConst() const;
		std::vector<double>		EdgeLens();
		std::vector<double>		GetInternalEdgeLens();
		std::vector<double>		GetExternalEdgeLens();
		preorder_iterator		begin();
		preorder_iterator		end();

        TreeNode * FindNodeBySplit(const Split & s);
		std::vector<TreeNode *> GetTips();

        unsigned                NumInternalNodesStored();
        unsigned                NumTipNodesStored();

		const TreeID & 			getTreeID() const;


		// Predicates
		//
		bool					RootValid() const;
		bool					IsRooted() const;
		bool					IsPolytomous() const;
		bool					HasEdgeLens() const;
		bool					PreorderDirty() const;
		bool					TipNumbersSetUsingNames() const;

        // Modifiers
#       if 0//
        void                    setRooted();
        void                    setUnrooted();
#       else
        void					setRootedness(bool rooted);
#       endif


		// Utilities
		//
		std::string				KeyToEdges();
		void					replaceEdgeLens(std::vector<double> new_edgelens);
		unsigned 				robinsonFoulds(TreeShPtr other);
		void 					buildTreeID();
        void                    Forget();
		void					Clear(bool forget = false);
		void					DeleteAllStoredNodes();
        void                    stripNodeNames();
        void                    renumberInternalNodes(unsigned start_at);
		void					BuildFromString(const std::string & newick, bool zero_based_tips = false); // throws XPhylogeny
		void					RectifyNumbers(std::vector<std::string> name_vector); // throws XPhylogeny
		void					RectifyNames(std::vector<std::string> name_vector); // throws XPhylogeny
		double					internalEdgeLenSum();
		double					externalEdgeLenSum();
		double					EdgeLenSum();
		double					calcTotalHeight();
		void					SetAllEdgeLens(double v);
		void					ScaleAllEdgeLens(double scaling_factor);
		void					RecalcAllSplits(unsigned max_nbits);
		void					RerootAtThisTip(TreeNode * nd);
		void					RerootAtThisInternal(TreeNode * nd);
		void					RerootAtTip(unsigned num);
		void					RefreshPreorder(TreeNode * nd = NULL) const;
		std::string &			AppendNewick(std::string &, unsigned float_precision = 5, bool useNumbers=false, bool showSupportAsNodeName = false, bool showSupportAsEdgeLength = false);
		std::string				MakeNewick(unsigned float_precision = 5);
		std::string				MakeNumberedNewick(unsigned float_precision = 5);
		std::string				MakeNewickForRefDist(unsigned float_precision = 5);
		bool					SetNumberFromName(TreeNode * nd, std::set<unsigned> & used);
        unsigned                FindTipByName(std::string tipname);
		TreeNodeVec				GetNodesWithEdges();
		void					SelectAllNodes();
		void					UnselectAllNodes();
		void					clearAdHocTreeNodePointers();
		//TreeNode *				FindMRCA(unsigned tip1, unsigned tip2);
		void					Ladderize(bool right);
		unsigned 				deroot();

		// Debugging
		//
		std::string				DebugWalkTree(bool preorder = true, unsigned verbosity = 0);
		bool					DebugCheckTree(bool allowDegTwo, bool checkDataPointers, int verbosity) const;
		void					DebugHere(std::string s);
		void		 			debugListTree();

        void                    debugMode(bool turn_on);

	protected:
		TreeNode * 				AllocNewNode();
		TreeNode *				FindTipNode(unsigned num);
		void					RerootHelper(TreeNode *m, TreeNode *t);
		void					GetNextNewickToken(const std::string &newick, unsigned start_pos);
		TreeNode *				GetNewNode();
		TreeNode *				PopLeafNode();
		TreeNode *				PopInternalNode();
		void					StoreInternalNode(TreeNode * u);
        void                    StoreLeafNode(TreeNode * u);
		void					Reserve(unsigned n);
		void					InvalidateID();
		void					InvalidateNodeCounts();
		void					RefreshNodeCounts();
		void					InvalidateTreeID();
		//void					RefreshTreeID();

		TreeNode *				FindLeftSib(TreeNode * start);
		TreeNode *				FindRightmostChild(TreeNode * start);
		TreeNode *				FindLastPreorderInClade(TreeNode * start);
		void					DetachSubtree(TreeNode * s);
		void					InsertSubtree(TreeNode * s, TreeNode * u, bool on_right,  TreeNode * targetSib = NULL);

		void 					recDebugListTree(TreeNode *p, int nindent);

        void SetFirstPreorder(TreeNode *nd)
            {
            assert(firstPreorder == 0L);
            firstPreorder = nd;
            }

        void MirrorTopology(Tree &source);
        void RebuildTopologyFromMirror(const Tree & source);
	protected:

		TreeID					tree_id;			/**< A vector of splits that uniquely identify the tree topology */
		bool					treeid_valid;			/**< True if the tree_id data member is valid; if false, call RefreshTreeID to make it valid again */
		TreeNodeStack			tipStorage;			    /**< A stack of pointers to (tip) TreeNode objects */
		TreeNodeStack			internalNodeStorage;	/**< A stack of pointers to (internal) TreeNode objects */
		mutable TreeNode *		firstPreorder;			/**< Pointer to the first preorder node (equals last postorder node) (mutable because it is not kept up-to-date, and may have to be recalculated on the fly)*/
		mutable TreeNode *		lastPreorder;			/**< Pointer to the last preorder node (equals first postorder node) (mutable because it is not kept up-to-date, and may have to be recalculated on the fly)*/
		unsigned				nTips;					/**< Total number of tip (degree = 1) nodes in the tree */
		unsigned				nInternals;				/**< Total number of internal (degree > 1) nodes in the tree */
		bool					hasEdgeLens;			/**< True if edge lengths have been specified */
		bool					isRooted;				/**< True if the tree is rooted */
		mutable bool			preorderDirty;			/**< Set to false when preorder traversal pointers are set, but a function should set to true if it modifies the tree and does not leave the preorder traversal pointers valid */
		bool					nodeCountsValid;		/**< If false, causes functions that depend on accurate node counts, such as GetNTips(), GetNNodes(), GetNInternals() and GetNObservables(), to recompute nTips and nInternals */
		bool					numbers_from_names;		/**< True if tip node numbers were set using tip node names in tree description. */

	private:

		std::string				workspace;				/**< Used by GetNextNewickToken for storing tokens read from newick tree descriptions */
		friend class FocalTreeTopoProbCalculator;
	};

typedef boost::shared_ptr<Tree> TreeShPtr;

}	// namespace phycas

//#include "phycas/src/basic_tree.inl"

#endif
