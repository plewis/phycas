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
#include <limits>
#include "topo_prior_calculator.hpp"
#include "basic_tree.hpp"
#include "basic_lot.hpp"
#include "tree_manip.hpp"
#include "probability_distribution.hpp"
#include "ncl/nxsallocatematrix.h"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Returns -log(number of unrooted binary tree topologies), where number of unrooted binary tree topologies equals
|   (2n-5)!! (n is number of tips).
*/
double TopoProbCalculator::GetLnTopoProb(TreeShPtr t)
    {
    // Note: this function assumes t is unrooted!
	double ntips = (double)(t->GetNObservables());
	PHYCAS_ASSERT(ntips >= 4.0);
	double ntips_minus_three = ntips - 3.0;
	double numer = lgamma(2.0*ntips - 5.0 + 1.0);
	double denom1 = ntips_minus_three*log(2.0);
	double denom2 = lgamma(ntips_minus_three + 1.0);
	double log_num_binary_trees = numer - denom1 - denom2;
    return -log_num_binary_trees;
    }

const double kLog2 = 0.6931471805599453;    // log 2  <-- 2 trees at maximum RF distance for 4-taxon tree (3 trees total)
const double kLog10 = 2.302585092994046;    // log 10 <-- 10 trees at maximum RF distance for 5-taxon tree (15 trees total)
const double kLog68 = 4.219507705176107;    // log 68 <-- 68 trees at maximum RF distance for symmetric 6-taxon tree (105 trees total)
const double kLog74 = 4.30406509320417;     // log 74 <-- 74 trees at maximum RF distance for pectinate 6-taxon tree (105 trees total)

std::vector<double> bs_b;
std::vector<double> bs_rs;

bool kAllowDegTwoInDebugCheck = false;

bool verboseMode = false;

double FocalTreeTopoProbCalculator::GetLnTopoProb(TreeShPtr t)
    {
    Tree & tree_ref = *(t.get());
    return CalcTopologyLnProb(tree_ref, false).first;
    }

// like CalcTopologyLnProb, but easier to export to Python
std::vector<double> FocalTreeTopoProbCalculator::CalcTopologyLnProbVector(Tree &tree, bool calcEdgeLenLnProb)
    {
    std::pair<double, double> r = CalcTopologyLnProb(tree, calcEdgeLenLnProb);
    std::vector<double> v(1, r.first);
    if (calcEdgeLenLnProb)
        v.push_back(r.second);
    return v;
    }

void FocalTreeTopoProbCalculator::SetDefaultEdgeLenDist(ProbDistShPtr d)
    {
    this->defEdgeLenDist = d;
    }

void FocalTreeTopoProbCalculator::DoTreeChecks(Tree & destTree, bool doRefreshesFirst, const char * tag) const
    {
    if (verboseMode && tag)
        std:: cerr << "DoTreeChecks " << tag << ":\n";
    if (doRefreshesFirst)
        {
        destTree.preorderDirty = true;
        destTree.RefreshPreorder();
        if (verboseMode)
            std::cerr << "Done Refreshing the preorder\n";
        destTree.RefreshNodeCounts();
        if (verboseMode)
            std::cerr << "Done Refreshing node counts\n";
        }
    destTree.debugMode(true);
    destTree.DebugCheckTree(kAllowDegTwoInDebugCheck, false, (verboseMode ? 2 : 0));
    if (verboseMode)
        std::cerr << "Done DebugCheckTree\n";
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Samples a tree using the reference distribution described in Holder, MT, Lewis, PO, and Swofford, DL. 2013.
|   Variable tree topology stepping-stone marginal likelihood estimation. In Chen, M.-H., Kuo, L, and Lewis PO (eds.),
|   Bayesian Phylogenetics: Methods, Algorithms, and Applications. Chapman/Hall.
*/
void FocalTreeTopoProbCalculator::SampleTree(TreeShPtr destTree, LotShPtr rng) const
    {
    if (verboseMode)
        std::cerr << "Here we are in SampleTree\n";
    DoTreeChecks(*destTree, true, "Incoming destTree");
    DoTreeChecks(*focalTree, true, "focalTree");

    // Get vector of tips sorted by the node number
    std::vector<TreeNode *> destTips = destTree->GetTips();
    if (verboseMode)
        {
        unsigned j = 0;
        for (std::vector<TreeNode *>::const_iterator dtIt = destTips.begin(); dtIt != destTips.end(); ++dtIt, ++j)
            {
            std::cerr << " destTip[" << j << "] = ";
            if (*dtIt)
                std::cerr << (*dtIt)->GetNodeNumber() << "\n";
            else
                std::cerr << "NULL" << "\n";
            }
        }

    // Create focalToDest, which maps nodes in focalTree to nodes in destTree
    std::map<TreeNode *, TreeNode *> focalToDest;

    // Add the root tip node to focalToDest map
    TreeNode * destNd = destTree->GetFirstPreorder();
    TreeNode * fnd = focalTree->GetFirstPreorder();
    if (destNd->GetNodeNumber() != fnd->GetNodeNumber())
        {
        // reroot destTree at the same node as focalTree
        PHYCAS_ASSERT(fnd->IsTipRoot());
        destTree->RerootAtTip(fnd->GetNodeNumber());
        DoTreeChecks(*destTree, true, "Incoming destTree after rerooting");
        destNd = destTree->GetFirstPreorder();
        }

    focalToDest[fnd] = destNd;

    // Add the immediate descendant of the root tip node to focalToDest map
    destNd = destNd->GetNextPreorder();
    fnd = fnd->GetNextPreorder();
    focalToDest[fnd] = destNd;

    // Collapse all internal edges in destTree, creating a star tree
    TreeNode * toCollapse = destNd->GetNextPreorder();
    while (toCollapse)
        {
        TreeNode * nextDestNd = toCollapse->GetNextPreorder();
        if (!toCollapse->IsExternalEdge())
            {
            toCollapse->CollapseEdge();
            destTree->StoreInternalNode(toCollapse);
            }
        toCollapse = nextDestNd;
        }
    destTree->nInternals = 1;
    DoTreeChecks(*destTree, true, "after collapse");

    // Detach all children above the hub node in destTree by setting hub node's lChild pointer to NULL
    destNd->lChild = 0L;

    // Move current node in focalTree to left child of basalmost internal node
    fnd = fnd->GetNextPreorder();

    // Create polytomyMap, which maps each polytomy in destTree to the node in focalTree that failed to be included
    std::map<TreeNode *, TreeNode *> polytomyMap;

    if (verboseMode)
        std::cerr << "About to enter while loop...\n";
    kAllowDegTwoInDebugCheck = true;

    // Proceed to visit all other nodes in focalTree in preorder fashion
    // For tips in focalTree, add them as left child to the corresponding parental node in destTree
    // For internals in focalTree, randomly choose whether to include a corresponding node in destTree
    // based on the clade probability stored in the guise of a branch length in focalTree.
    // If focalTree internal node (focalNode) is chosen,
    //   - add a new corresponding internal node to the destTree
    // If focalTree internal node fails to be chosen,
    //   - select focalNode
    //   - polytomyMap[destNode] = focalNode
    //   - focalToDest[focalNode] = destNode
    while (fnd)
        {
        TreeNode * newPar = focalToDest[fnd->GetParent()];
        PHYCAS_ASSERT(newPar);
        if (verboseMode)
            std::cerr << "NewPar # = " << newPar->GetNodeNumber();
        if (fnd->IsExternalEdge())
            {
            if (verboseMode)
                std::cerr << "looking for leaf # " << fnd->GetNodeNumber() << "\n";
            TreeNode * leafNd = destTips.at(fnd->GetNodeNumber());
            PHYCAS_ASSERT(leafNd);
            if (verboseMode)
                std::cerr << "add new leafNd " << leafNd->GetNodeNumber() << "\n";
            leafNd->rSib = 0L;
            newPar->AddChild(leafNd);
            fnd->SetIsSelected(false);
            }
        else
            {
            double probInclusion = fnd->GetEdgeLen(); // hackety-hack: clade prob stored as branch length
            if (rng->Uniform() < probInclusion)
                {
                TreeNode * newInternal = destTree->GetNewNode();
                if (verboseMode)
                    std::cerr << "add new internal " << newInternal->GetNodeNumber() << "\n";
                newPar->AddChild(newInternal);
                destNd = newInternal;
                newInternal->split = fnd->split;
                fnd->SetIsSelected(false);
                }
            else
                {
                if (verboseMode)
                    std::cerr << "skipping edge\n";
                if (polytomyMap.find(newPar) == polytomyMap.end())
                    polytomyMap[newPar] = fnd->GetParent();
                destNd = newPar;
                fnd->SetIsSelected(true);
                }
            focalToDest[fnd] = destNd;
            }
        DoTreeChecks(*destTree, true, "after edge inclusion decision");
        fnd = fnd->GetNextPreorder();
        if (verboseMode)
            std::cerr << "Moving to the next node...\n";
        }
    if (verboseMode)
        std::cerr << "Done walking over the focal tree\n";
    kAllowDegTwoInDebugCheck = false;

    DoTreeChecks(*destTree, true, "after adding the split from focal tree");

    // Randomly resolve each polytomy created in destTree
    for (std::map<TreeNode *, TreeNode *>::const_iterator pmIt = polytomyMap.begin(); pmIt != polytomyMap.end(); ++pmIt)
        {
        TreeNode * destPolytomy = pmIt->first;
        TreeNode * focalPolytomy = pmIt->second;
        if (verboseMode)
            std::cerr << "resolving a polytomy...\n";
        ResolveToAvoidSharedSplits(destTree, destPolytomy, focalPolytomy, rng);
        }
    if (verboseMode)
        std::cerr << "Done resolving the variate\n";
    destTree->preorderDirty = true;
    destTree->RefreshPreorder();
    if (verboseMode)
        std::cerr << "Done Refreshing the preorder\n";
    destTree->RefreshNodeCounts();
    DoTreeChecks(*destTree, true, "after topology is finished");

    //std::cerr << "\nPOLPOL: -----------------------------\n";
    preorder_iterator ndIt = destTree->begin();
    ++ndIt;
    for (; ndIt != destTree->end(); ++ndIt)
        {
        TreeNode & nd = *ndIt;
        ProbDistShPtr d = GetEdgeLenProbDistForSplit(nd.split);
        PHYCAS_ASSERT(d);
        const double e = d->Sample();
        nd.SetEdgeLen(e);
        //ProbabilityDistribution * dptr = &*d;	//POLPOL
        //std::cerr << boost::str(boost::format("POLPOL: Drew %g for node %d using %s @") % e % nd.GetNodeNumber() % d->GetDistributionDescription()) << (dptr) << "\n";
        if (verboseMode)
            std::cerr << "Drew " << e << " for node # = " << nd.GetNodeNumber() << "\n";
        }
    if (verboseMode)
        std::cerr << "Returning from SampleTree\n";
    }


// inserts all "flagged" splits that are connected to correspondingFocalNode into
//  splitSet
void FocalTreeTopoProbCalculator::FindTabuSplitsFromSelectedNodes(TreeNode * correspondingFocalNode, std::set<Split> & splitSet) const
{
    if (verboseMode)
        std::cerr << "in FindTabuSplitsFromSelectedNodes\n";

    std::stack<TreeNode *> ndStack;
    PHYCAS_ASSERT(correspondingFocalNode);
    TreeNode * curr = correspondingFocalNode->GetLeftChild();
    PHYCAS_ASSERT(curr);
    if (!curr)
        return;
    PHYCAS_ASSERT(curr->rSib);
    if (curr->rSib)
        ndStack.push(curr->rSib);
    for (;;)
        {
        if (curr->IsSelected())
            {
            splitSet.insert(curr->split);
            curr = curr->GetLeftChild();
            PHYCAS_ASSERT(curr);
            PHYCAS_ASSERT(curr->rSib);
            if (curr->rSib)
                ndStack.push(curr->rSib);
            }
        else
            {
            if (ndStack.empty())
                return;
            curr = ndStack.top();
            ndStack.pop();
            }
        }

}

void FocalTreeTopoProbCalculator::ResolveToAvoidSharedSplits(
        TreeShPtr destTree,
        TreeNode * destPolytomy,
        TreeNode * correspondingFocalNode,
        LotShPtr rng) const
    {
    if (verboseMode)
        std::cerr << "Here we are in ResolveToAvoidSharedSplits\n";
    TreeNode * newInternal = destTree->GetNewNode();
    TreeNode * origLChild = destPolytomy->GetLeftChild();
    const unsigned polytomyDeg = 1 + destPolytomy->CountChildren();

    if (polytomyDeg == 4)
        {
        // correspondingFocalNode is parent of the node (now selected) that failed to be included in destTree
        // Start by finding the selected descendant of correspondingFocalNode
        destPolytomy->lChild = 0L;
        TreeNode * selectedNode = correspondingFocalNode->GetLeftChild();
        if (!selectedNode->IsSelected())
            {
            selectedNode = selectedNode->GetRightSib();
            PHYCAS_ASSERT(selectedNode->IsSelected());
            }

        // Choose one of the two children of selectedNode to leave out of the newly created clade
        // The one chosen will be selToMoveDown
        TreeNode * selLeft = selectedNode->GetLeftChild();
        TreeNode * selRight = selLeft->GetRightSib();
        TreeNode * selToMoveDown = 0L;
        selToMoveDown = (rng->Uniform() < 0.5 ? selLeft : selRight);
        if (verboseMode)
            {
            std::cerr << "selToMoveDown.split = " << selToMoveDown->split.CreatePatternRepresentation() << '\n';
            std::cerr << "origLChild.split = " << origLChild->split.CreatePatternRepresentation() << '\n';
            std::cerr << "origLChild->rSib.split = " << origLChild->rSib->split.CreatePatternRepresentation() << '\n';
            std::cerr << "origLChild->rSib->rSib.split = " << origLChild->rSib->rSib->split.CreatePatternRepresentation() << '\n';
            }

        // Now find selToMoveDown's counterpart in destTree and call it destLower
        // The two nodes that will be part of the new clade will be destU1 and destU2
        TreeNode *destLower, *destU1, *destU2;
        if (selToMoveDown->split == origLChild->split)
            {
            destLower = origLChild;
            destU1 = destLower->rSib;
            destU2 = destU1->rSib;
            }
        else
            {
            destU1 = origLChild;
            if (destU1->rSib->split == selToMoveDown->split)
                {
                destLower = destU1->rSib;
                destU2 = destLower->rSib;
                }
            else
                {
                destU2 = destU1->rSib;
                destLower = destU2->rSib;
                PHYCAS_ASSERT(destLower->split == selToMoveDown->split);
                }
            }
        destLower->rSib = 0L;
        destU1->rSib = 0L;
        destU2->rSib = 0L;

        destPolytomy->AddChild(destLower);
        destPolytomy->AddChild(newInternal);
        newInternal->AddChild(destU1);
        newInternal->AddChild(destU2);
        if (verboseMode)
            std::cerr << "leaving ResolveToAvoidSharedSplits in the deg = 4 branch\n";
        }
    else    // polytomyDeg > 4
        {
        if (verboseMode)
            std::cerr << "not in the deg = 4 branch\n";
        std::set<Split> tabuSplits;
        FindTabuSplitsFromSelectedNodes(correspondingFocalNode, tabuSplits);
        if (verboseMode)
            std::cerr << "back from FindTabuSplitsFromSelectedNodes\n";

        std::vector<TreeNode *> polytomyChildren = destPolytomy->GetChildren();
        for (;;)
            {

            std::vector<TreeNode *> newNodes = RandomlyResolve(destTree, destPolytomy, polytomyChildren, rng);
            if (!HasTabuSplit(destTree, destPolytomy, newNodes, tabuSplits))
                {
                if (verboseMode)
                    std::cerr << "Returning a tree because RandomlyResolve did not create a tree with a tabu split.\n";
                return;
                }
            std::vector<TreeNode *>::const_iterator nnIt = newNodes.begin();
            for (; nnIt != newNodes.end(); ++nnIt)
                {
                TreeNode * nd = *nnIt;
                nd->CollapseEdge();
                destTree->StoreInternalNode(nd);
                }
            if (verboseMode)
                std::cerr << "Tabu split generated trying again.\n";
            }
        }
    }

// returns newly added internals
std::vector<TreeNode *> FocalTreeTopoProbCalculator::RandomlyResolve(
        TreeShPtr destTree,
        TreeNode * destPolytomy,
        const std::vector<TreeNode *> & polytomyChildren,
        LotShPtr rng) const
    {
    if (verboseMode)
        std::cerr << "Entering RandomlyResolve with " << polytomyChildren.size() << " polytomy children. destPolytomy= << "<< destPolytomy->GetNodeNumber() << "\n";
    DoTreeChecks(*destTree, true, "rr");

    std::vector<TreeNode *> newNodes;
    if (polytomyChildren.size() < 3)
        return newNodes;
    std::vector<TreeNode *> currentEdges;
    currentEdges.reserve(2*polytomyChildren.size() + 1);
    currentEdges.push_back(destPolytomy);
    destPolytomy->lChild = polytomyChildren[0];
    destPolytomy->lChild->rSib = polytomyChildren[1];
    TreeNode * firstAdded = destPolytomy->lChild;
    TreeNode * secondAdded = firstAdded->rSib;
    currentEdges.push_back(firstAdded);
    currentEdges.push_back(secondAdded);
    secondAdded->rSib = 0L;
    DoTreeChecks(*destTree, true, "after trimming to 2 children");

    for (unsigned currentChildInd = 2; currentChildInd < polytomyChildren.size(); ++currentChildInd)
        {
        TreeNode * nextChild = polytomyChildren[currentChildInd];
        PHYCAS_ASSERT(nextChild);
        PHYCAS_ASSERT(nextChild != secondAdded && nextChild != firstAdded);

        unsigned edgeInd = rng->SampleUInt(currentEdges.size());
        TreeNode * newNd = destTree->GetNewNode();
        PHYCAS_ASSERT(newNd);
        if (verboseMode)
            std::cerr << "Adding to edge " << edgeInd << "\n";
        if (edgeInd == 0)
            {
            if (verboseMode)
                std::cerr << "In adding to edge 0 branch\n";
            TreeNode * edgeToBisect = destPolytomy->lChild;
            newNd->lChild = edgeToBisect;
            newNd->rSib = nextChild;
            newNd->par = destPolytomy;

            edgeToBisect->par = newNd;
            edgeToBisect->rSib->par = newNd;

            nextChild->par = destPolytomy;
            nextChild->rSib = 0L;

            destPolytomy->lChild = newNd;
            }
        else
            {
            TreeNode * edgeToBisect = currentEdges[edgeInd];
            PHYCAS_ASSERT(edgeToBisect);
            TreeNode * parent = edgeToBisect->GetParent();
            PHYCAS_ASSERT(parent);
            PHYCAS_ASSERT(parent->lChild);
            if (parent->lChild == edgeToBisect)
                {
                if (verboseMode)
                    std::cerr << "  Was lchild\n";
                parent->lChild = newNd;
                }
            else
                {
                if (verboseMode)
                    std::cerr << "  Was rchild\n";
                PHYCAS_ASSERT(parent->lChild->rSib);
                PHYCAS_ASSERT(parent->lChild->rSib->rSib == 0L);
                PHYCAS_ASSERT(parent->lChild->rSib == edgeToBisect);
                parent->lChild->rSib = newNd;
                }
            newNd->lChild = edgeToBisect;
            newNd->rSib = edgeToBisect->rSib;
            newNd->par = parent;

            edgeToBisect->rSib = nextChild;
            edgeToBisect->par = newNd;

            nextChild->par = newNd;
            nextChild->rSib = 0L;
            }

        if (verboseMode)
            std::cerr << "  resetting split\n";
        newNd->split.Reset();
        if (verboseMode)
            std::cerr << "  storing nodes\n";
        newNodes.push_back(newNd);
        currentEdges.push_back(newNd);
        currentEdges.push_back(nextChild);
        DoTreeChecks(*destTree, true, "another node added in RandomlyResolve");
        }

    std::vector<TreeNode *>::const_iterator pcIt = polytomyChildren.begin();
    for (; pcIt != polytomyChildren.end(); ++pcIt)
        {
        TreeNode * pc = *pcIt;
        TreeNode * anc = pc->par;
        while (anc != destPolytomy)
            {
            anc->split.CombineWith(pc->split);
            anc = anc->par;
            }
        }
    return newNodes;
    }

bool FocalTreeTopoProbCalculator::HasTabuSplit(
        TreeShPtr destTree,
        TreeNode * destPolytomy,
        const std::vector<TreeNode *> & newNodes,
        std::set<Split> & tabuSplits) const
    {
    if (verboseMode)
        std::cerr << "Entering HasTabuSplit.\n";

    std::vector<TreeNode *>::const_iterator nnIt = newNodes.begin();
    for (; nnIt != newNodes.end(); ++nnIt)
        {
        if (tabuSplits.find((*nnIt)->split) != tabuSplits.end())
            return true;
        }
    return false;
    }

double FocalTreeTopoProbCalculator::LnNumBinaryTrees()
	{
    PHYCAS_ASSERT(ntips > 3);
    double ntips_minus_three = (double)(ntips-3);
    double numer = lgamma(2.0*ntips - 5.0);
    double denom1 = ntips_minus_three*log(2.0);
    double denom2 = lgamma(ntips_minus_three);
    return numer - denom1 - denom2;
	}

void FocalTreeTopoProbCalculator::SetEdgeLenDist(const Split &s, ProbDistShPtr edgeLenDist)
    {
    splitToEdgeLenDistMap[s] = edgeLenDist;
    }

/*----------------------------------------------------------------------------------------------------------------------
|
|	Computes the number of possible rooted trees for 1, 2, ..., taxa.n taxa, storing the results as logs.
|
|	Make sure that 'b' is allocated for n+1 elements.
*/
static void enumTrees(double *b, unsigned n, bool rooted)
	{
	unsigned k, m;

	b[1] = b[2] = b[3] = 0.0;
	for (k = 3 + !rooted, m = 3; k <= n; k++, m += 2)
		b[k] = log((double)m) + b[k-1];
	}

struct bs_data
	{
	int							nv;		/**< number of internal edges above node v */
	ScopedTwoDMatrix<double>	r;		/**< R(v,.,.) for recursion */
	};

#define NV(p)	((bs_data *)((p)->ptr))->nv
#define R(p)	((bs_data *)((p)->ptr))->r


void zeroRMatrix(TreeNode *p, unsigned prevNV)
    {
    unsigned nv1 = NV(p) + 1;
    //nv1 = nmax + 1; //@TEMP
    if (prevNV < nv1)
        R(p).Initialize(nv1, nv1);
    double ** alias = R(p).GetAlias();
    for (unsigned i = 0; i < nv1 ; ++i)
        {
        for (unsigned j = 0; j < nv1 ; ++j)
            {
            alias[i][j] = 0.0;
            }
        }
    //std::cerr << "Zeroed matrix for node " << p->GetNodeNumber() << " ptr R(p)= " << R(p).GetAlias() << "\n" ;
    }

FocalTreeTopoProbCalculator::FocalTreeTopoProbCalculator(
        TreeShPtr t)
    :focalTree(t)
    {
    assert(bool(focalTree));
    ntips = focalTree->GetNObservables();
    focalTree->RefreshPreorder();
    focalTree->RefreshNodeCounts();
    focalTree->RecalcAllSplits(ntips);
    DoTreeChecks(*focalTree, false, "In FocalTreeTopoProbCalculator ctor");
    TreeNode * nd = focalTree->GetFirstPreorder();
    assert(nd->IsTipRoot());
    nd = nd->GetNextPreorder();
    while (nd)
        {
        if (!nd->IsExternalEdge())
            {
            double split_prob = nd->GetEdgeLen();
            assert(split_prob > 0.0 && split_prob < 1.0);
            Split & s = nd->GetSplit();
            if (s.IsBitSet(0))
                s.InvertSplit();
            splitToProbMap[s] = split_prob;
            }
        nd = nd->GetNextPreorder();
        }
    focalTree->RefreshPreorder();
    focalTree->RefreshNodeCounts();
    //buildScratchTree();

	unsigned n = focalTree->GetNObservables();
	int nmax = n - 3;

	bs_b.resize(n + 2);
	bs_rs.resize(nmax + 1);

	enumTrees(&bs_b[0], n, false);
	for (unsigned i = 0; i <= n; i++)
		bs_b[i] = exp(bs_b[i]);

    //bool prevVerboseMode = verboseMode;
    //verboseMode = true;
    //DoTreeChecks(*focalTree, true, "ctor");
    //verboseMode = prevVerboseMode;
	for (TreeNode *p = focalTree->GetLastPreorder(); p != NULL; p = p->GetNextPostorder())
		{
		PHYCAS_ASSERT(p->ptr == NULL);
		p->ptr = new bs_data();
		if (p->IsTip())
			NV(p) = 0;
		else
			{
			TreeNode * q = p->GetLeftChild();
			TreeNode * r = q->GetRightSib();
			NV(p) = NV(q) + NV(r) + !q->IsTip() + !r->IsTip();
			zeroRMatrix(p, 0);
            }
		}
    }

FocalTreeTopoProbCalculator::~FocalTreeTopoProbCalculator()
	{
	for (preorder_iterator p = focalTree->begin(); p != focalTree->end(); p++)
		{
		delete (bs_data *)p->ptr;
		p->ptr = NULL;
		}
	bs_b.clear();
	bs_rs.clear();
	}

/*
void FocalTreeTopoProbCalculator::buildScratchTree()
    {
    scratchTree.MirrorTopology(*focalTree);
    }

void Tree::MirrorTopology(Tree &source)
    {
    Clear();
    TreeNode * fnd = source.GetFirstPreorder();
    assert(fnd);
    assert(fnd->IsTipRoot());
    TreeNode * tmpSNd = GetNewNode();
    SetFirstPreorder(tmpSNd);


    tmpSNd->SetCorrespondingNode(fnd);
    fnd->SetCorrespondingNode(tmpSNd);
    tmpSNd->CopyNonPointerFields(*fnd);
    tmpSNd->par = 0L;
    tmpSNd->lChild = 0L;
    tmpSNd->rSib = 0L;
    tmpSNd->nextPreorder = 0L;
    tmpSNd->prevPreorder = 0L;

    TreeNode * snd = GetRoot();
    assert(snd);
    assert(snd->IsTipRoot());

    TreeNode * prevS = snd;
    fnd = fnd->GetNextPreorder();
    while (fnd)
        {
        tmpSNd = GetNewNode();
        tmpSNd->SetCorrespondingNode(fnd);
        fnd->SetCorrespondingNode(tmpSNd);
        tmpSNd->CopyNonPointerFields(*fnd);
        tmpSNd->lChild = 0L;
        tmpSNd->rSib = 0L;
        tmpSNd->nextPreorder = 0L;
        prevS->nextPreorder = tmpSNd;
        tmpSNd->prevPreorder = prevS;
        TreeNode * sPar = fnd->par->GetCorrespondingNode();
        assert(sPar);
        sPar->AddChild(tmpSNd);
        prevS = tmpSNd;
        fnd = fnd->GetNextPreorder();
        }

    }
// Assumes that all nodes in focalTree have a "CorrespondingNode" that is allocated in
//  "this."
// does not alter anything except the "key" navigational pointers (par, lChild, rSib), and the split field
void Tree::RebuildTopologyFromMirror(const Tree & source)
    {
    const TreeNode * fnd = source.GetFirstPreorderConst();
    assert(fnd);
    assert(fnd->IsTipRoot());
    TreeNode * tmpSNd = fnd->GetCorrespondingNode();
    SetFirstPreorder(tmpSNd);
    tmpSNd->par = 0L;
    tmpSNd->lChild = (fnd->lChild ? fnd->lChild->GetCorrespondingNode() : 0L);
    tmpSNd->rSib = 0L;
    tmpSNd->nextPreorder = tmpSNd->lChild;
    tmpSNd->prevPreorder = 0L;
    TreeNode * prevS = tmpSNd;
    fnd = fnd->GetNextPreorderConst();
    while (fnd)
        {
        tmpSNd = fnd->GetCorrespondingNode();
        tmpSNd->lChild = (fnd->lChild ? fnd->lChild->GetCorrespondingNode() : 0L);
        tmpSNd->rSib = (fnd->rSib ? fnd->rSib->GetCorrespondingNode() : 0L);
        tmpSNd->par = (fnd->par ? fnd->par->GetCorrespondingNode() : 0L);
        tmpSNd->nextPreorder = (fnd->nextPreorder ? fnd->nextPreorder->GetCorrespondingNode() : 0L);
        tmpSNd->prevPreorder = prevS;
        tmpSNd->split = fnd->split;
        TreeNode * sPar = fnd->par->GetCorrespondingNode();
        assert(sPar);
        prevS = tmpSNd;
        fnd = fnd->GetNextPreorderConst();
        }

    }
*/

// Replaces "this" with its children in the tree
// only fixes the "key" navigational pointers (par, lChild, rSib)
// does not alter "this" node's pointers
void TreeNode::CollapseEdge()
    {
    TreeNode * p = this->par;
    if (!p)
        return;
    TreeNode * ls = this->FindLeftSib();
    TreeNode * c = this->lChild;
    if (ls == 0L)
        p->lChild = c;
    else
        ls->rSib = c;
    if (!c)
        return;
    for (;;c = c->rSib)
        {
        c->par = p;
        if (c->rSib == 0L)
            {
            c->rSib = this->rSib;
            break;
            }
        }
    }

#if 1 //begin new DLS code

//#include "paup_util.h"


/*----------------------------------------------------------------------------------------------------------------------
|
|	Returns beta(m) as defined by Bryant and Steel.
*/
double bs_beta(unsigned m)
	{
	return bs_b[m + 3];
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Returns -1 raised to the power of k.
*/
double powerOfMinusOne(int k)
	{
	return (k & 1) ? -1.0 : 1.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Returns true iff p's edge either terminates at a leaf or at a node that is not selected.
*/
bool isLogicalExternalEdge(TreeNode *p)
	{
	return p->IsTip() || !p->IsSelected();
	}


/*----------------------------------------------------------------------------------------------------------------------
|
*/
bool isLogicalInternalEdge(TreeNode *p)
	{
	return p->IsInternal() && p->IsSelected();
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Counts the number of cherries (pairs of leaves on adjecent external edges) on a tree.
*/
// commented out to avoid "defined but not used" compiler warnings
//static int countCherries(TreeNode *pStart)
//	{
//	int ncherries = 0;
//	for (TreeNode *p = pStart; p != NULL; p = p->GetNextPreorder())
//		{
//		if (!isLogicalExternalEdge(p) && isLogicalExternalEdge(p->GetLeftChild())
//		  && isLogicalExternalEdge(p->GetLeftChild()->GetRightSib()))
//			ncherries++;
//		}
//	if (isLogicalExternalEdge(pStart->GetLeftChild()) || isLogicalExternalEdge(pStart->GetLeftChild()->GetRightSib()))
//		ncherries++;
//
//	return ncherries;
//	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Computes Poisson approximation to the proportion of trees at a given RF distance from a particular tree (shape).
|	The tree information is contained in the lambdaT value, which is a function of the number of cherries on the
|	tree.
*/
// commented out to avoid "defined but not used" compiler warnings
//static double poissonApprox(unsigned n, double m, double lambdaT)		/* m = RF distance, multiple of 2 */
//	{
//	double s = n - 3 - m/2;
//	double s_factorial = exp(lgamma(s + 1.0));
//	return exp(-lambdaT)*pow(lambdaT, s)/s_factorial;
//	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Computes distribution of Robinson-Foulds tree distances using method of Bryant and Steel (2009).
*/
double FocalTreeTopoProbCalculator::countDistancesUsingBryantSteel(
  TreeNode *ff, /**< root node of polytomous subtree */
  unsigned n)   /**< number of leaves in embedded tree */
	{
	TreeNode * v0 = ff;
	//std::cerr << "v0->num = " << v0->GetNodeNumber() <<'\n';
	//const bool prevVerboseMode = verboseMode;
	//verboseMode = true;
    DoTreeChecks(*focalTree, true, "countDistancesUsingBryantSteel");
	//verboseMode = prevVerboseMode;

	PHYCAS_ASSERT(!v0->IsTip());
	int n_max = n - 3;

    // Push all nodes in the embedded "tree" to postorderNodeStack, including
    // ff and descendants (but not ff ancestors)
	PHYCAS_ASSERT(postorderNodeStack.empty());
	postorderNodeStack.push(ff);
	for (TreeNode *v = ff->GetNextPreorder(); v != NULL; v = v->GetNextPreorder())
		{
		postorderNodeStack.push(v);
		if (isLogicalExternalEdge(v) && v->IsInternal())
			v = v->FindLastPreorderInClade();
		}

	while (!postorderNodeStack.empty())
		{
		TreeNode *v = postorderNodeStack.top();
		postorderNodeStack.pop();
		if ((v == ff) || isLogicalInternalEdge(v))
			{
            unsigned prevNV = NV(v);
			TreeNode * v1 = v->GetLeftChild();
			TreeNode * v2 = v1->GetRightSib();
			NV(v) = NV(v1) + NV(v2) + isLogicalInternalEdge(v1) + isLogicalInternalEdge(v2);
            zeroRMatrix(v, prevNV);
			int       n_v = NV(v);
			double ** r_v = R(v).GetAlias();

			unsigned nNonleafChildren = isLogicalInternalEdge(v1) + isLogicalInternalEdge(v2);

			/* Lemma 1 */
			r_v[0][n_v] = bs_beta(n_v);
			if (nNonleafChildren == 1)
				{
				if (isLogicalExternalEdge(v1))
					v1 = v2;
				PHYCAS_ASSERT(isLogicalInternalEdge(v1));

				int       n_v1 = NV(v1);
				double ** r_v1 = R(v1).GetAlias();

				for (int s = 1; s <= n_v; s++)
					{
					/* Lemma 2 step 3 for k=0 */
					double sum = 0.0;
					if (s-1 <= n_v1)
						{
						for (int k1 = 0; k1 <= n_v1; k1++)
							sum += r_v1[s-1][k1];
						}
					r_v[s][0] = sum;

					/* Lemma 2 step 3 for k>0 */
					if (s <= n_v1)
						{
						for (int k = 1; k <= n_v; k++)
							r_v[s][k] = r_v1[s][k-1]*(2*k + 1);
						}
					}
				}
			else if (nNonleafChildren == 2)
				{
				int       n_v1 = NV(v1);
				int       n_v2 = NV(v2);
				double ** r_v1 = R(v1).GetAlias();
				double ** r_v2 = R(v2).GetAlias();

				for (int s = 1; s <= n_v; s++)
					{
					/* Lemma 2 step 4 (k=0) */
					double sum_v = 0.0;
					for (int s1 = 0; s1 <= s-2; s1++)
						{
						if ((s1 <= n_v1) && (s-2-s1 <= n_v2))
							{
							double sum_v1 = 0.0;
							for (int k1 = 0; k1 <= n_v1; k1++)
								sum_v1 += r_v1[s1][k1];

							double sum_v2 = 0.0;
							for (int k2 = 0; k2 <= n_v2; k2++)
								sum_v2 += r_v2[s-2-s1][k2];

							sum_v += sum_v1*sum_v2;
							}
						}
					r_v[s][0] = sum_v;

					/* Lemma 2 step 5 (k>0) */
					for (int k = 1; k <= n_v; k++)
						{
						double beta_ratio = bs_beta(k)/bs_beta(k-1);

						double sum1 = 0.0;
						if (k-1 <= n_v2)
							{
							for (int s1 = 0; s1 <= std::min(s - 1, n_v1); s1++)
								{
								if (s-1-s1 <= n_v2)
									{
									double sum_v1 = 0.0;
									for (int k1 = 0; k1 <= n_v1; k1++)
										sum_v1 += r_v1[s1][k1];
									sum1 += sum_v1 * r_v2[s-1-s1][k-1] * beta_ratio;
									}
								}
							}

						double sum2 = 0.0;
						if (k-1 <= n_v1)
							{
							for (int s2 = 0; s2 <= std::min(s-1, n_v2); s2++)	//min probably not needed
								{
								if (s-1-s2 <= n_v1)
									{
									double sum_v2 = 0.0;
									for (int k2 = 0; k2 <= n_v2; k2++)
										sum_v2 += r_v2[s2][k2];
									sum2 += sum_v2 * r_v1[s-1-s2][k-1] * beta_ratio;
									}
								}
							}

						double sum3 = 0.0;
						for (int s1 = 0; s1 <= std::min(s, n_v1); s1++)
							{
							if (s-s1 <= n_v2)
								{
								for (int k1 = 0; k1 <= k-2; k1++)
									{
									if ((k1 <= n_v1) && (k-2-k1 <= n_v2))
										sum3 += r_v1[s1][k1] * r_v2[s-s1][k-2-k1]
												* bs_beta(k)/(bs_beta(k1)*bs_beta(k-2-k1));
									}
								}
							}

						r_v[s][k] = sum1 + sum2 + sum3;
						}
					}
				}
			}
		else
		    {
		    NV(v) = 0;
		    }

		}

	double ** r_v0 = R(v0).GetAlias();
	for (int s = 0; s <= n_max; s++)
		{
		double rs = 0.0;
		int n_v0 = NV(v0);
		for (int k = 0; k <= n_v0-s; k++)	/* POL changed "k <= n_v0" to "k <= n_v0-s" */	/* Mark's fix; Bryant and Steel incorrectly used "k <= s" */
			rs += r_v0[s][k];   // POL crashes here sometimes - bad access
		bs_rs[s] = rs;
		}

	double bmT = 0.0;
	for (int s = 0; s <= n_max; s++)
		bmT += bs_rs[s]*powerOfMinusOne(s);

//	myprintf("\nNumber of trees at maximum R-F distance calculated using Bryant-Steel method:\n\n");
//	myprintf("%4d%12g\n", d_max, bmT);

	return log(bmT);
	}

#endif	//end new DLS code

double FocalTreeTopoProbCalculator::CalcLnNumTreesMaxDistFromTreeInSelectedRegion(
  TreeNode *firstFork,  /**< the root node of the polytomy  */
  unsigned numLeaves)   /**< the number of leaves in the embedded tree */
	{
    PHYCAS_ASSERT(firstFork);
    PHYCAS_ASSERT(!firstFork->IsTip());
    double lnNumTrees = 0.0;

    if (numLeaves < 7)
        {
        if (numLeaves == 4)
            lnNumTrees = kLog2;
        else if (numLeaves == 5)
            lnNumTrees = kLog10;
		else
			{
			const TreeNode * lc = firstFork->lChild;
			const TreeNode * rc = lc->rSib;
			if (lc->IsSelected())
				{
				if (rc->IsSelected())
					lnNumTrees = kLog74; // pectinate 6-leaf
				else
					{
					lc = lc->lChild;
					rc = lc->rSib;
					if (lc->IsSelected() && rc->IsSelected())
						lnNumTrees = kLog68; // symmetric 6-leaf
					else
						lnNumTrees = kLog74; // pectinate 6-leaf
					}
				}
			else
				{
				PHYCAS_ASSERT(rc->IsSelected());
				lc = rc->lChild;
				rc = lc->rSib;
				if (lc->IsSelected() && rc->IsSelected())
					lnNumTrees = kLog68; // symmetric 6-leaf
				else
					lnNumTrees = kLog74; // pectinate 6-leaf
				}
			}
        PHYCAS_ASSERT(firstFork);
        PHYCAS_ASSERT(!firstFork->IsTip());
        return lnNumTrees;
/*        std::cerr << "lnNumTrees = " << lnNumTrees << '\n';
        const double countedLnNumTrees = countDistancesUsingBryantSteel(firstFork, numLeaves);
        std::cerr << "countedLnNumTrees = " << countedLnNumTrees << '\n';
        std::cerr << "kLog2 = " << kLog2 << '\n';
        std::cerr << "kLog10 = " << kLog10 << '\n';
        std::cerr << "kLog68 = " << kLog68 << '\n';
        std::cerr << "kLog74 = " << kLog74 << '\n';
		PHYCAS_ASSERT(fabs(lnNumTrees - countedLnNumTrees) < 1.0e-6);	//TEMP
*/
        }
	else
		lnNumTrees = countDistancesUsingBryantSteel(firstFork, numLeaves);

	return lnNumTrees;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the probability distribution for the edge length defined by the supplied split `s'. If `s' is found in
|   `splitToEdgeLenDistMap' then that split-specific edge length reference distribution will be returned; otherwise,
|   `defEdgeLenDist' will be returned. Throws XProbDist exception if `defEdgeLenDist' is needed but not defined.
*/
ProbDistShPtr FocalTreeTopoProbCalculator::GetEdgeLenProbDistForSplit(const Split & s) const
    {
    std::map<Split, ProbDistShPtr>::const_iterator sIt = splitToEdgeLenDistMap.find(s);
    if (sIt == splitToEdgeLenDistMap.end())
        {
        if (defEdgeLenDist)
            return defEdgeLenDist;
        PHYCAS_ASSERT(false);
    	throw XProbDist("Split not found, and no default edge length distribution has been specified");
        }
    return sIt->second;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the log probability density of the edge length defined by the supplied split `s'. If `s' is found in
|   `splitToEdgeLenDistMap' then that split-specific edge length reference distribution will be used to compute the
|   probability density; otherwise, `defEdgeLenDist' will be used to compute the density. Throws XProbDist exception if
|   `defEdgeLenDist' is needed but not defined.
*/
double FocalTreeTopoProbCalculator::LnEdgeLenProbForSplit(const Split & s, const double b) const
    {
    ProbDistShPtr p = GetEdgeLenProbDistForSplit(s);
    double log_density = p->GetLnPDF(b);
    //std::cerr << boost::str(boost::format("print %s.getLnPDF(%g),%g # %s") % p->GetDistributionDescription() % b % log_density % s.CreatePatternRepresentation()) << std::endl;
    return log_density;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Computes the log of the probability of the supplied `testTree' using the reference distribution described in
|   Holder, MT, Lewis, PO, and Swofford, DL. 2013. Variable tree topology stepping-stone marginal likelihood estimation.
|   In Chen, M.-H., Kuo, L, and Lewis PO (eds.), Bayesian Phylogenetics: Methods, Algorithms, and Applications.
|   Chapman/Hall.
*/
std::pair<double, double> FocalTreeTopoProbCalculator::CalcTopologyLnProb(Tree & testTree, bool calcEdgeLenLnProb)
    {
    testTree.RecalcAllSplits(ntips);
    const TreeID & testTreeID = testTree.getTreeID();
    const TreeID::const_iterator ttIDIt = testTreeID.end();
    //scratchTree.RebuildTopologyFromMirror(*focalTree);

    //std::cerr << "testTree  = " << testTree.MakeNewick() << std::endl;
    //std::cerr << "focalTree = " << focalTree->MakeNewick() << std::endl;

    // Start with fnd (focal node) set to the internal node just above the root
    TreeNode * fnd = focalTree->GetFirstPreorder();
    fnd->SetIsSelected(false);
    fnd = fnd->GetNextPreorder();

    double lnNumerator = 0.0;       // this will be the numerator of the first value of the pair returned
    double lnDenominator = 0.0;     // this will be the denominator of the first value of the pair returned
    double lnEdgeLenProb = 0.0;     // this will be the second value of the pair returned

    // create a map for storing polytomies
    std::map<TreeNode *, std::vector<TreeNode *> > polytomyToCollapsed;

    // First, compute the log of the numerator of the probability of the tree topology defined by testTree
    while (fnd)
        {
        fnd->SetIsSelected(false);

        //std::string splitstr = fnd->GetSplitConst().CreatePatternRepresentation();
        //if (splitstr == "----*-*--*************--*********")
        //    {
        //    std::cerr << "target split: " << fnd->GetSplitConst().CreateIdRepresentation() << " (" << fnd->GetSplitConst().CreatePatternRepresentation() << ")" << std::endl;
        //    std::cerr << "testTreeID: " << std::endl;
        //    for (TreeID::const_iterator idIt = testTreeID.begin(); idIt != testTreeID.end(); idIt++)
        //        {
        //        std::cerr << "  " << idIt->CreateIdRepresentation() << " (" << idIt->CreatePatternRepresentation() << ")" << std::endl;
        //        }
        //    }

        if (!fnd->IsExternalEdge())
            {
            // See if the split defined by fnd is present in testTree
            const TreeID::const_iterator idIt = testTreeID.find(fnd->GetSplitConst());
            if (idIt == ttIDIt)
                {
                //std::cerr << "NOT FOUND: focal split = " << splitstr << std::endl;

                // split defined by fnd NOT found in testTree, so add log[1-p(s)] to numerator
                lnNumerator += log(1.0 - fnd->GetEdgeLen()); // could store log(1-p) in support and log(p) in edge_len to cut down on logs

                // select fnd because this node is in part of the tree that represents a resolution of a polytomy in Tu
                fnd->SetIsSelected(true);
                TreeNode * p = fnd->GetParent();
                if (p->IsSelected())
                    fnd->adhocTreeNodePtr = p->adhocTreeNodePtr; // using adhocTreeNodePtr to store the "deepest node" that will represent the polytomy.
                else
                    fnd->adhocTreeNodePtr = p;

                if (verboseMode)
                    std::cerr << "Adding node " << fnd->GetNodeNumber() << " to polytomy at " << fnd->adhocTreeNodePtr->GetNodeNumber() << '\n';

                std::vector<TreeNode *> & vc = polytomyToCollapsed[fnd->adhocTreeNodePtr];
                vc.push_back(fnd);
                }
            else
                {
                //std::cerr << "* FOUND *  focal split = " << splitstr << std::endl;

                // split defined by fnd WAS found in testTree, so add log[p(s)] to numerator
                lnNumerator += log(fnd->GetEdgeLen()); // could store log(1-p) in support and log(p) in edge_len to cut down on logs
                }
            }
        fnd = fnd->GetNextPreorder();
        }

    // If requested, compute lnEdgeLenProb (second value in returned pair)
    if (calcEdgeLenLnProb)
        {
        // Loop through all nodes in testTree and build up lnEdgeLenProb as we go
        preorder_iterator testNdIt = testTree.begin();
        ++testNdIt;
        for (; testNdIt != testTree.end(); ++testNdIt)
            {
            double currLnEdgeLenProb = LnEdgeLenProbForSplit(testNdIt->split, testNdIt->GetEdgeLen());
            lnEdgeLenProb += currLnEdgeLenProb;
            }
        }

    // Second, compute the log of the denominator of the probability of the tree topology defined by testTree
    for (std::map<TreeNode *, std::vector<TreeNode *> >::const_iterator pNdIt = polytomyToCollapsed.begin(); pNdIt != polytomyToCollapsed.end(); ++pNdIt)
        {
        TreeNode * polytomyNd = pNdIt->first;
        PHYCAS_ASSERT(polytomyNd->IsInternal());
        PHYCAS_ASSERT(!polytomyNd->IsTip());
        const std::vector<TreeNode *> & vc = pNdIt->second;

        // unselect all nodes adjacent to nodes in vc list
        for (std::vector<TreeNode *>::const_iterator ndIt = vc.begin(); ndIt != vc.end(); ++ndIt)
            {
            TreeNode * adj = (*ndIt)->GetLeftChild();
            if (adj)
                {
                adj->SetIsSelected(false);
                if (adj->rSib)
                    adj->GetRightSib()->SetIsSelected(false);
                }
            adj = (*ndIt)->GetParent();
            if (adj)
                adj->SetIsSelected(false);
            }

        // select all nodes in vc list
        for (std::vector<TreeNode *>::const_iterator ndIt = vc.begin(); ndIt != vc.end(); ++ndIt)
            (*ndIt)->SetIsSelected(true);

        const unsigned numLeaves = 3 + vc.size();
        PHYCAS_ASSERT(polytomyNd->IsInternal());
        PHYCAS_ASSERT(!polytomyNd->IsTip());
        lnDenominator += CalcLnNumTreesMaxDistFromTreeInSelectedRegion(polytomyNd, numLeaves);

        if (verboseMode)
            std::cerr << "lnDenominator = " << lnDenominator << "\t numLeaves = " << numLeaves << std::endl;

        // "unflag" edges that are collapsed in the consensus
        for (std::vector<TreeNode *>::const_iterator ndIt = vc.begin(); ndIt != vc.end(); ++ndIt)
            (*ndIt)->SetIsSelected(false);
        }

    //std::cerr << boost::str(boost::format("--> lnNumerator = %g, lnDenominator = %g, lnEdgeLenProb = %g, total = %g") % lnNumerator % lnDenominator % lnEdgeLenProb % (lnNumerator - lnDenominator + lnEdgeLenProb)) << std::endl;	//POLPOL
    focalTree->UnselectAllNodes();
    focalTree->clearAdHocTreeNodePointers();

    return std::pair<double, double>(lnNumerator - lnDenominator, lnEdgeLenProb);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `counts' vector for the supplied number of internal nodes (`n') using the method outlined by Joe
|	Felsenstein in his 2004 book and also in Felsenstein (1978) and Felsenstein (1981).
|
|	Felsenstein, J. 1978. The number of evolutionary trees. Syst. Zool. 27: 27-33.
|	Felsenstein, J. 1981. Syst. Zool. 30: 122.
|   Felsenstein, J. 2004. Inferring Phylogeny. Sinauer, Sunderland, Massachusetts.
|
|   Below I've reproduced the table from Felsenstein 2004 illustrating how to calculate the number of possible rooted
|   tree topologies for any number of internal nodes:
|>
|                      number of taxa in rooted tree            counts    nfactors
|             2       3       4       5        6        7         8
|	n     +-------+-------+-------+--------+--------+--------+---------+   +---+  This table shows results for 8
|	u  1  |   1   |   1   |   1   |   1    |    1   |    1   |    1    |   | 0 |  taxa and rooted trees. Columns
|	m     +-------+-------+-------+--------+--------+--------+---------+   +---+  corresponding to 2-7 taxa would have
|	b  2          |   3   |   10  |   25   |   56   |   119  |   246   |   | 0 |  been replaced by the column
|	e             +-------+-------+--------+--------+--------+---------+   +---+  corresponding to 8 taxa. This column
|	r  3                  |   15  |  105   |   490  |   1918 |   6825  |   | 0 |  is now stored in the `counts' vector.
|	                      +-------+--------+--------+--------+---------+   +---+  The scaling factor used in this case
|	o  4                          | 0.0105 | 0.1260 |  0.945 |  5.6980 |   | 1 |  was 10,000. The number of times this
|	f                             +--------+--------+--------+---------+   +---+  factor was applied is stored in the
|	   5                                   | 0.0945 | 1.7325 | 19.0575 |   | 1 |  nfactors vector. Thus, the "count"
|	n                                      +--------+--------+---------+   +---+  for 5 internal nodes is not really
|	o  6                                            | 1.0395 | 27.0270 |   | 1 |  19.0575, but instead 190575:
|	d                                               +--------+---------+   +---+  count[4]*(scaling_factor)^nfactors[4]
|	e  7                                                     | 13.5135 |   | 1 |
|	s                                                        +---------+   +---+
|>
|   The RecalcCountsAndPriorsImpl function works from left to right, calculating each column in turn. Within a column,
|   it works down. Each cell except the first in a given column requires knowledge of two cells from the column to the
|   left (the cell to its immediate left as well as the cell above the cell to its immediate left). This is because
|   in order to know how many tree topologies there are for N taxa, one needs to know how many places a taxon could be
|   added to a tree with one fewer taxa. The N-1 taxon tree could have the same number of internal nodes as the new
|   tree (the new taxon was added to an existing node, either creating a new polytomy or enlarging an existing one), or
|   the N-1 taxon tree could have one fewer internal nodes (in which case the new taxon inserts a new node).
*/
void PolytomyTopoPriorCalculator::RecalcCountsAndPriorsImpl(
  unsigned n) /**< is the number of internal nodes (equals ntax - 1 for rooted trees and ntax - 2 for unrooted trees) */
	{
    if (is_resolution_class_prior)
        counts_dirty = true;
    if (counts_dirty)
        {
        double scaling_factor = exp(log_scaling_factor);

	    counts.clear();
	    counts.push_back(1.0); // counts are always 1 for m = 1

        int last_factor = 0;
	    nfactors.clear();
	    nfactors.push_back(last_factor); // never need to scale this one

        // temporary variables
        double a, b, c;
        double max_log_double = log(DBL_MAX);
        double epsilon = scaling_factor/10.0; // value arbitrary, but must be larger than zero and less than scaling_factor

        // Compute the vector of counts for the number of internal nodes specified
        // This is the main loop over columns. z is the number of taxa minus 1 for rooted trees, and is the number
        // of taxa minus 2 for unrooted trees
	    for (unsigned z = 2; z <= n; ++z)
		    {
		    counts.push_back(0.0);  // this column is one element longer than the column to its left
            nfactors.push_back(last_factor);
            b = counts[0];  // counts[0] is always 1.0 because there is only one star tree topology

            // This is the loop over rows within the current column. m + 1 is the number of internal nodes.
            for (unsigned m = 1; m < z; ++m)
                {
                unsigned num_internal_nodes = m + 1;
                double diff = (double)(nfactors[m - 1] - nfactors[m]);
                double log_factor = diff*log_scaling_factor;
                if (log_factor >= max_log_double)
                    {
                    std::cerr << "Oops! log_factor >= max_log_double" << std::endl;
                    }
                PHYCAS_ASSERT(log_factor < max_log_double);
                a = b*exp(log_factor);
                b = counts[m];
                c = a*((double)(z + num_internal_nodes - 1));
                if (num_internal_nodes < z)
                    {
                    c += b*(double)num_internal_nodes;
                    }
                if (c > scaling_factor)
                    {
                    //unsigned prev_nfactors = nfactors[m];
                    double incr = floor(log(c - epsilon)/log_scaling_factor);
                    nfactors[m] += (unsigned)incr;
                    last_factor = nfactors[m];
                    counts[m] = exp(log(c) - incr*log_scaling_factor);
                    b = exp(log(b) - incr*log_scaling_factor);
                    }
                else
                    counts[m] = c;
                }
		    }

        // Now compute the log of the total number of tree topologies over all possible resolution classes
        // (i.e. number of internal nodes)
        // Begin by creating a vector of log counts and finding the largest value (this will be factored out
        // to avoid overflow)
        std::vector<double> v;
        unsigned sz = (unsigned)nfactors.size();
        PHYCAS_ASSERT(sz == counts.size());
        double max_log_count = 0.0;
        for (unsigned i = 0; i < sz; ++i)
            {
            double num_factors = (double)nfactors[i];
            double log_count = num_factors*log_scaling_factor + log(counts[i]);
            if (log_count > max_log_count)
                max_log_count = log_count;
            v.push_back(log_count);
            }

        // Compute log sum of counts by factoring out the largest count. Underflow will occur, but only for
        // counts that are so much smaller than the dominant counts that the underflow can be ignored for
        // all practical purposes
        double sum = 0.0;
        for (std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it)
            {
            double diff = (*it) - max_log_count;
            sum += exp(diff);
            }
        PHYCAS_ASSERT(sum > 0.0);
        log_total_count = log(sum) + max_log_count;
        counts_dirty = false;
        }
    else
        {
        nfactors.clear();
        counts.clear();
        counts_dirty = true;    // ensures that counts will be calcualated if someone asks for one, say by calling GetCount
        }

	// Recalculate the `topology_prior' vector too
	RecalcPriorsImpl();

	topo_priors_dirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `topology_prior' vector for the supplied number of internal nodes (`n'). The element `topology_prior'[m]
|	is the natural logarithm of the unnormalized prior probability of a (rooted/unrooted) tree with `ntax' taxa and m
|	internal nodes. The rooted/unrooted status is determined by the state of the data member `is_rooted'. The element
|	`topology_prior'[0] is the natural log of the normalizing constant. Normally, only unnormalized values are needed
|	because MCMC deals in prior ratios, but if for some reason the normalized prior probability is needed, it can be
|	computed as exp{`topology_prior'[m] - `topology_prior'[0]}. This function requires `counts' to be correct if using
|   the resolution class prior. Thus, never call RecalcPriorsImpl directly, only invoke it indirectly by calling the
|   function RecalcCountsAndPriorsImpl.
|
|	Consider an unrooted tree with ntax = 6. Such a tree has 4 internal nodes, and calling RecalcPriorsImpl
|	would yield the following if `is_resolution_class_prior' is false:
|>
|	     unnormalized
|	m    polytomy prior     C = 2
|	----------------------------------
|	1        C^3              8		     topology_prior[1] = ln(8)  = 2.079
|	2        C^2              4		     topology_prior[2] = ln(4)  = 1.386
|	3        C^1              2		     topology_prior[3] = ln(2)  = 0.693
|	4        C^0              1		     topology_prior[4] = ln(1)  = 0.000
|	----------------------------------
|	       C^4 - 1         16 - 1
|	       -------         ------ = 15   topology_prior[0] = ln(15) = 2.708
|	        C - 1           2 - 1
|>
|	If instead `is_resolution_class_prior' is true, we have:
|>
|	              unnormalized
|	              resolution
|	m    counts   class prior        C = 2
|	---------------------------------------------
|	1         1     (C^3)/1         8/1   = 8.000   topology_prior[1] = ln(8.000) =  2.079
|	2        25     (C^2)/25        4/25  = 0.160   topology_prior[2] = ln(0.160) = -1.833
|	3       105     (C^1)/105       2/105 = 0.019   topology_prior[3] = ln(0.019) = -3.963
|	4       105     (C^0)/105       1/105 = 0.010   topology_prior[4] = ln(0.010) = -4.605
|	---------------------------------------------
|	            (no easy formula)           8.189   topology_prior[0] = ln(8.189) =  2.103
|>
*/
void PolytomyTopoPriorCalculator::RecalcPriorsImpl()
	{
	topology_prior.clear();
	topology_prior.push_back(0.0);	// This will hold the normalizing constant in the end

	// Figure out the maximum possible value for m, the number of internal nodes
	unsigned maxm = ntax - (is_rooted ? 1 : 2);

    if (is_resolution_class_prior)
        {
	    // counts vector should have length equal to maxm if everything is ok
	    PHYCAS_ASSERT(maxm == (unsigned)counts.size());

	    double total = 0.0;
	    double logC = std::log(C);
	    for (unsigned m = 1; m <= maxm; ++m)
		    {
		    double logCterm = (double)(maxm - m)*logC;
		    double log_count_m = std::log(counts[m - 1]) + log_scaling_factor*(double)nfactors[m - 1];
		    double log_v = logCterm - log_count_m;
		    total += std::exp(log_v);
		    topology_prior.push_back(log_v);
		    }
	    topology_prior[0] = std::log(total);
        }
    else
        {
	    double total = 0.0;
	    double logC = std::log(C);
	    for (unsigned m = 1; m <= maxm; ++m)
		    {
		    double logCterm = (double)(maxm - m)*logC;
		    total += std::exp(logCterm);
		    topology_prior.push_back(logCterm);
		    }
	    topology_prior[0] = std::log(total);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the number of trees having `n' taxa and `m' internal nodes. Calls RecalcCountsAndPriors
|   function if `n' is not equal to `ntax'. Assumes `m' is greater than 0. If `is_rooted' is true, assumes `m' is less
|   than `ntax'. If `is_rooted' is false, assumes `m' less than `ntax' - 1.
*/
double PolytomyTopoPriorCalculator::GetLnCount(
  unsigned n,	/**< is the number of taxa */
  unsigned m)	/**< is the number of internal nodes */
	{
	PHYCAS_ASSERT((is_rooted && (m < n)) || (!is_rooted && (m < n - 1)));
	if (n != ntax)
		SetNTax(n);
	if (counts_dirty)
		Reset();
    double nf = (double)(nfactors[m - 1]);
    double log_count = nf*log_scaling_factor + log(counts[m - 1]);
	return log_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of saturated (i.e. fully-resolved and thus having as many internal nodes as possible) trees
|	of `n' taxa. Calls RecalcCountsAndPriors function if `n' is not equal to `ntax'.
*/
double PolytomyTopoPriorCalculator::GetLnSaturatedCount(
  unsigned n)	/**< is the number of taxa */
	{
	if (n != ntax)
		SetNTax(n);
	if (counts_dirty)
		Reset();
    unsigned last = (unsigned)(counts.size() - 1);
    double nf = (double)(nfactors[last]);
    double log_count = nf*log_scaling_factor + log(counts[last]);
	return log_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the total number of trees for `n' taxa, including all resolution classes from the star
|   tree to fully resolved (saturated) trees. Calls RecalcCountsAndPriors function if `n' is not equal to `ntax' or if
|   not using the resolution class prior (in which case counts have not been calculated).
*/
double PolytomyTopoPriorCalculator::GetLnTotalCount(
  unsigned n)	/**< is the number of taxa */
	{
	if (n != ntax)
		SetNTax(n);
	if (counts_dirty)
		Reset();
	return log_total_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a vector of realized resolution class priors from the values in the `topology_prior' vector. If
|	`topo_priors_dirty' is true, it recomputes the `topology_prior' vectors first. The mth element of the
|	returned vector is set to T_{n,m}*`topology_prior'[m] for m > 0. The 0th element of the returned vector holds the
|	normalization constant (sum of all other elements). This function is not efficient because it is intended only to
|	be used for providing information to the user on request. Table 2, p. 248, in the "Polytomies and Bayesian
|	Phylogenetic Inference" paper (Lewis, P. O., M. T. Holder and K. E. Holsinger. 2005. Systematic Biology 54(2):
|	241-253) presented (normalized) values from this vector.
*/
std::vector<double> PolytomyTopoPriorCalculator::GetRealizedResClassPriorsVect()
	{
	if (!is_resolution_class_prior)
        counts_dirty = true;
	if (topo_priors_dirty || counts_dirty)
		Reset();

	std::vector<double> v;
	v.reserve(topology_prior.size());
	v.push_back(0.0);

	unsigned sz = (unsigned)topology_prior.size();

    // First loop will be to determine largest value, which will be factored out
    // the second time through so that the total does not overflow
    double log_factored_out = 0.0;
	for (unsigned i = 1; i < sz; ++i)
		{
        double c = counts[i - 1];
        double nf = (double)nfactors[i - 1];
        double log_Tnm = log(c) + log_scaling_factor*nf;
		double log_prior = log_Tnm + topology_prior[i];
		v.push_back(log_prior);
        if (log_prior > log_factored_out)
            log_factored_out = log_prior;
		}

    // Now we can compute the total
    double total = 0.0;
    std::vector<double>::const_iterator it = v.begin();
    for (++it; it != v.end(); ++it)
		{
		total += exp((*it) - log_factored_out);
		}
	v[0] = log(total) + log_factored_out;

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a vector in which the element having index m (i = 0, 1, ..., max. num. internal nodes) represents the
|   natural logarithm of the number of tree topologies having m internal nodes. If `counts_dirty' is true, it recomputes
|   the `counts' vectors first. The 0th element of the returned vector holds the natural log of the total number of tree
|   topologies (log of the sum of all other elements).
*/
std::vector<double> PolytomyTopoPriorCalculator::GetLnCounts()
	{
    if (is_resolution_class_prior)
        counts_dirty = true;
	if (counts_dirty)
		Reset();

	std::vector<double> v;
	unsigned sz = ntax - (is_rooted ? 0 : 1);
	v.reserve(sz);
	v.push_back(log_total_count);

	//@POL could use a version of the transform algorithm here
	for (unsigned i = 1; i < sz; ++i)
		{
        double log_Tnm = log(counts[i - 1]) + log_scaling_factor*(double)(nfactors[i - 1]);
		v.push_back(log_Tnm);
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Samples a resolution class (i.e. number of internal nodes) from the realized resolution class distribution. This
|   function is not very efficient because it calls PolytomyTopoPriorCalculator::GetRealizedResClassPriorsVect, resulting in an
|   unnecessary vector copy operation.
*/
unsigned PolytomyTopoPriorCalculator::sample(
  LotShPtr rng) /**< is the random number generator to use for sampling */
	{
    std::vector<double> v = GetRealizedResClassPriorsVect();
#if 0
    // Geez, what was I thinking when I wrote this?!
    double u = rng->Uniform();
    double logu = (u > 0.0 ? log(u) : -DBL_MAX);
    double log_x = logu + v[0];
    double cum = 0.0;
    for (unsigned i = 1; i < v.size(); ++i)
        {
        cum += v[1];
        if (log_x <= cum)
            return i;
        }
#else
    double u = rng->Uniform();
    //std::cerr << "PolytomyTopoPriorCalculator::sample: seed = " << rng->GetSeed() << ", u = " << u << std::endl;
    //std::cerr << "v[0]  = " << v[0] << std::endl;
    double z = v[0];
    double cum = 0.0;
    for (unsigned i = 1; i < v.size(); ++i)
        {
        //std::cerr << "v[" << i << "]  = " << v[i] << std::endl;
        cum += exp(v[i] - z);
        if (u <= cum)
            return i;
        }
#endif
    PHYCAS_ASSERT(0);
    return (unsigned)(v.size() - 1);
    }


// methods below here were formerly in topo_prior_calculator.inl

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `is_rooted' to false, `is_resolution_class_prior' to true, `C' to 1.0, `ntax' to 4, and
|   `topo_priors_dirty' to true.
*/
PolytomyTopoPriorCalculator::PolytomyTopoPriorCalculator()
	{
	topo_priors_dirty			= true;
	is_rooted					= false;
	is_resolution_class_prior	= true;
	C							= 1.0;
	ntax						= 4;
    counts_dirty                = true;
    log_scaling_factor          = 10.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor clears the vectors `counts', `nfactors' and `topology_prior'.
*/
PolytomyTopoPriorCalculator::~PolytomyTopoPriorCalculator()
	{
	counts.clear();
	topology_prior.clear();
	nfactors.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `new_ntax' differs from `ntax', sets `ntax' to `new_ntax' and sets `topo_priors_dirty'. Returns immediately
|	without taking any action if `ntax' equals `new_ntax'. Assumes that `new_ntax' is greater than 1 if `is_rooted' is
|	true, or that `new_ntax' is greater than 2 if `is_rooted' is false.
*/
void PolytomyTopoPriorCalculator::SetNTax(
  unsigned new_ntax)	/**< is the new number of taxa */
	{
	if (ntax != new_ntax)
		{
		// Set ntax to the new value
		PHYCAS_ASSERT(new_ntax > (unsigned)(is_rooted ? 1 : 2));
		ntax = new_ntax;

		counts_dirty = true;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Forces recalculation of `polytomy_prior' if `is_resolution_class_prior' is false, and both `counts' and
|   `polytomy_prior' if `is_resolution_class_prior' is true (or if `counts_dirty' is true).
*/
void PolytomyTopoPriorCalculator::Reset()
	{
	unsigned num_internal_nodes = (is_rooted ? (ntax - 1) : (ntax - 2));
	RecalcCountsAndPriorsImpl(num_internal_nodes);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `ntax'.
*/
unsigned PolytomyTopoPriorCalculator::GetNTax() const
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_rooted' data member to true. There are more rooted than unrooted trees for the same value of `ntax', so
|	this setting is important when asking questions that require knowledge of the numbers of possible trees.
*/
void PolytomyTopoPriorCalculator::ChooseRooted()
	{
	if (!is_rooted)
		{
		is_rooted = true;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_rooted' data member to false. There are more rooted than unrooted trees for the same value of `ntax', so
|	this setting is important when asking questions that require knowledge of the numbers of possible trees.
*/
void PolytomyTopoPriorCalculator::ChooseUnrooted()
	{
	if (is_rooted)
		{
		is_rooted = false;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns copy of the `counts' vector, which contains in its (m-1)th element the number of tree topologies having
|   exactly m internal nodes. Note that you will need to also call GetNFactorsVect if there is any chance that some of
|   the counts are larger than exp(log_scaling_factor). In such cases, the actual log count is
|   log count = `nfactors[m]'*`log_scaling_factor' + log(`count[m - 1]')
*/
std::vector<double> PolytomyTopoPriorCalculator::GetCountsVect()
	{
	//@POL this function could be const were it not for lazy evaluation
	if (counts_dirty)
		Reset();
	return counts;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns copy of the `nfactors' vector, which contains in its (m-1)th element the number of times counts[m] has been
|   rescaled by dividing by the scaling factor (the log of which is `log_scaling_factor').
*/
std::vector<int> PolytomyTopoPriorCalculator::GetNFactorsVect()
	{
	//@POL this function could be const were it not for lazy evaluation
	if (counts_dirty)
		Reset();
	return nfactors;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_resolution_class_prior' data member to true.
*/
void PolytomyTopoPriorCalculator::ChooseResolutionClassPrior()
	{
	if (!is_resolution_class_prior)
		{
		is_resolution_class_prior = true;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `is_resolution_class_prior' data member to false.
*/
void PolytomyTopoPriorCalculator::ChoosePolytomyPrior()
	{
	if (is_resolution_class_prior)
		{
		is_resolution_class_prior = false;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `C' data member to `c'. Assumes `c' is greater than 0.0.
*/
void PolytomyTopoPriorCalculator::SetC(double c)
	{
	PHYCAS_ASSERT(c > 0.0);
	if (c != C)
		{
		C = c;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `C' data member.
*/
double PolytomyTopoPriorCalculator::GetC() const
	{
	return C;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `log_scaling_factor' data member to the supplied value `lnf'.
*/
void PolytomyTopoPriorCalculator::SetLnScalingFactor(double lnf)
	{
	PHYCAS_ASSERT(lnf > 0.0);
	if (lnf != log_scaling_factor)
		{
		log_scaling_factor = lnf;
		counts_dirty = true;
		topo_priors_dirty = true;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `log_scaling_factor' data member.
*/
double PolytomyTopoPriorCalculator::GetLnScalingFactor() const
	{
	return log_scaling_factor;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns copy of the `topology_prior' vector, which contains in its mth element the unnormalized prior for tree
|	topologies having exactly m internal nodes. The 0th element of `topology_prior' holds the normalizing constant.
*/
std::vector<double> PolytomyTopoPriorCalculator::GetTopoPriorVect()
	{
	//@POL this function could be const were it not for lazy evaluation
	if (topo_priors_dirty)
		Reset();
	return topology_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns result of call to GetLnTopologyPrior(m), where m is the number of internal nodes in `t'.
*/
double PolytomyTopoPriorCalculator::GetLnTopoProb(TreeShPtr t)
    {
    unsigned n = t->GetNTips();
	PHYCAS_ASSERT(n > 3);
    if (n != GetNTax())
        SetNTax(n);
    unsigned m = t->GetNInternals();
    double topo_prior = GetLnTopologyPrior(m);
    return topo_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the unnormalized topology prior. This represents the resolution class prior if
|	`is_resolution_class_prior' is true, otherwise it represents the polytomy prior. Assumes `m' is less than the
|	length of the `topology_prior' vector.
*/
double PolytomyTopoPriorCalculator::GetLnTopologyPrior(
  unsigned m)	/**< is the number of internal nodes */
	{
	if (topo_priors_dirty)
		Reset();
	PHYCAS_ASSERT(m < topology_prior.size());
	return topology_prior[m];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the normalized topology prior. This represents the resolution class prior if
|	`is_resolution_class_prior' is true, otherwise it represents the polytomy prior. Assumes `m' is less than the length
|	of the `topology_prior' vector. The log of the normalized topology prior is obtained as `topology_prior'[m] minus
|	`topology_prior'[0] (the 0th element of `topology_prior' holds the log of the normalization constant).
*/
double PolytomyTopoPriorCalculator::GetLnNormalizedTopologyPrior(
  unsigned m)	/**< is the number of internal nodes */
	{
	if (topo_priors_dirty)
		Reset();
	PHYCAS_ASSERT(m < topology_prior.size());
	return (topology_prior[m] - topology_prior[0]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the normalizing constant for the topology prior. This value is stored in
|	`topology_prior'[0].
*/
double PolytomyTopoPriorCalculator::GetLnNormConstant()
	{
	if (topo_priors_dirty)
		Reset();
	return topology_prior[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `is_resolution_class_prior'.
*/
bool PolytomyTopoPriorCalculator::IsResolutionClassPrior() const
	{
	return is_resolution_class_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `is_resolution_class_prior' is true, returns false; if `is_resolution_class_prior' is false, returns true.
*/
bool PolytomyTopoPriorCalculator::IsPolytomyPrior() const
	{
	return !is_resolution_class_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `is_rooted'.
*/
bool PolytomyTopoPriorCalculator::IsRooted() const
	{
	return is_rooted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `is_rooted' is true, returns false; if `is_rooted' is false, returns true.
*/
bool PolytomyTopoPriorCalculator::IsUnrooted() const
	{
	return !is_rooted;
	}

}	// namespace phycas
