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

#ifndef TREE_MANIP_HPP
#define TREE_MANIP_HPP

//POL This file is a modified version of phycasdev/phycas/tree/tree_manip.hpp
// I am gradually incorporating the functions needed by the new Phycas

#include "boost/shared_ptr.hpp"
#include "basic_tree.hpp"

namespace phycas
{

class TreeLengthDistribution;
typedef boost::shared_ptr<TreeLengthDistribution> TreeLengthDistributionShPtr;

class ProbabilityDistribution;
typedef boost::shared_ptr<ProbabilityDistribution> ProbDistShPtr;

class TreeNode;

typedef boost::shared_ptr<Tree> TreeShPtr;

typedef boost::shared_ptr<TreeID> TreeIDShPtr;

class Lot;
typedef boost::shared_ptr<Lot> LotShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	The TreeManip class manipulates Tree objects (e.g. topological rearrangements needed for tree searching or MCMC).
*/
class TreeManip
	{
	public:
					TreeManip();
					TreeManip(TreeShPtr t);
					~TreeManip();

		void		setTree(TreeShPtr t);

		void		starTree(unsigned ntips, ProbDistShPtr edge_len_dist);
		void		randomTree(unsigned ntips, LotShPtr rng, ProbDistShPtr edge_len_dist, bool yule);

        void        buildEquiprobTree(unsigned ntips, LotShPtr rng);
		void		equiprobTree(unsigned ntips, LotShPtr rng, ProbDistShPtr internal_edge_dist, ProbDistShPtr external_edge_dist);
        void        setRandomEdgeLensFromTreeLengthDist(TreeLengthDistributionShPtr tree_length_dist);

        void        deleteRandomInternalEdge(LotShPtr rng);

		void		setRandomEdgeLens(ProbDistShPtr d);
		void		setRandomInternalExternalEdgeLens(ProbDistShPtr internal_dist, ProbDistShPtr external_dist);

        void        buildTreeFromSplitVector(const std::vector<std::string> & split_vect, ProbDistShPtr edge_len_dist);

		void		NNISwap(TreeNode * swap1, TreeNode * swap2);
		void		NNISwapSpecial(TreeNode * swap1);

		enum InsertMode		/**< Determines whether nodes are inserted on right or left */
			{
			kOnRight,	/**< Insert new node as rightmost child (or as right sibling of targetSib, if targetSib specified) */
			kOnLeft		/**< Insert new node as lstmost child (or as left sibling of targetSib, if targetSib specified) */
			};

		TreeNode *	FindLeftSib(TreeNode * start)
		    {
		    PHYCAS_ASSERT(start != NULL);
            return start->FindLeftSib();
		    }
		TreeNode *	FindRightmostChild(TreeNode * start)
		    {
		    PHYCAS_ASSERT(start != NULL);
            return start->FindRightmostChild();
		    }
		TreeNode *	FindLastPreorderInClade(TreeNode * start);

        void addRightChild(TreeNode * parent_node, TreeNode * new_rightmost_child)
            {
            InsertSubtree(new_rightmost_child, parent_node, kOnRight);
            }

		void		DetachSubtree(TreeNode * s);
		void		InsertSubtree(TreeNode * s, TreeNode * u, InsertMode m, TreeNode * targetSib = NULL);

		void		DeleteLeaf(TreeNode * u, bool store_as_internal = false);
		void		InsertSubtreeIntoEdge(TreeNode * subtree, TreeNode * edge_nd);

		void		LChildToLSib(TreeNode * u, TreeNode * w);
		void		RChildToRSib(TreeNode * u, TreeNode * w);

		void		SibToChild(TreeNode * u, TreeNode * s, InsertMode m, TreeNode * targetSib = NULL);

    protected:

		TreeShPtr	tree;	/**< The tree to be manipulated */

	};

TreeNode * randomInternalAboveSubroot(const Tree &, Lot &);

} // namespace phycas

#endif
