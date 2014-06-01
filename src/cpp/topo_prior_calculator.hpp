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

#ifndef TOPO_PRIOR_CALCULATOR_HPP
#define TOPO_PRIOR_CALCULATOR_HPP

#include <boost/shared_ptr.hpp>

#include "split.hpp"
#include "basic_tree.hpp"

namespace phycas
{

class Tree;
typedef boost::shared_ptr<Tree> TreeShPtr;

class Lot;
typedef boost::shared_ptr<Lot> LotShPtr;

class ProbabilityDistribution;
typedef boost::shared_ptr<ProbabilityDistribution> ProbDistShPtr;

class TopoProbCalculator
    {
    public:
        TopoProbCalculator() {}
        virtual ~TopoProbCalculator() {}
        virtual double  GetLnTopoProb(TreeShPtr t);
    };

typedef boost::shared_ptr<TopoProbCalculator> TopoProbCalculatorShPtr;

class FocalTreeTopoProbCalculator: public TopoProbCalculator
    {
    public:
                                        FocalTreeTopoProbCalculator(TreeShPtr);
                                        ~FocalTreeTopoProbCalculator();

        void                            SampleTree(TreeShPtr, LotShPtr) const;
        std::pair<double, double>       CalcTopologyLnProb(Tree &, bool calcEdgeLenLnProb);
        virtual double                  GetLnTopoProb(TreeShPtr t);
        std::vector<double>             CalcTopologyLnProbVector(Tree &tree, bool calcEdgeLenLnProb);
        double                          CalcLnNumTreesMaxDistFromTreeInSelectedRegion(TreeNode *selectedFirstFork, unsigned numLeaves);
        void                            SetEdgeLenDist(const Split &, ProbDistShPtr);
        void                            SetDefaultEdgeLenDist(ProbDistShPtr d);
        ProbDistShPtr                   GetEdgeLenProbDistForSplit(const Split & s) const;
		double                          LnNumBinaryTrees();
        double                          CalcTopologyLnProb(Tree &) const;
        double                          CalcLnNumTreesMaxDistFromTreeInSelectedRegion(const TreeNode *, unsigned numLeaves) const;

    protected:
        TreeShPtr                       focalTree;              /**< the reference tree T^{\ast} in Holder et al. 2013 */
        std::map<Split, double>         splitToProbMap;
        unsigned                        ntips;
        ProbDistShPtr                   defEdgeLenDist;
        std::map<Split, ProbDistShPtr>  splitToEdgeLenDistMap;  /**< the edge length reference distribution for well-supported edges in T^{\ast} */

		mutable std::stack<TreeNode *>  postorderNodeStack;
        mutable Tree                    scratchTree;
        mutable std::set<TreeNode *>    omittedNodes;

        void                            ResolveToAvoidSharedSplits(TreeShPtr dest, TreeNode * destPolytomy, TreeNode * correspondingFocalNode, LotShPtr rng) const;
        void                            FindTabuSplitsFromSelectedNodes(TreeNode * correspondingFocalNode, std::set<Split> & splitSet) const;
        std::vector<TreeNode *>         RandomlyResolve(TreeShPtr dest, TreeNode * destPolytomy, const std::vector<TreeNode *> & polytomyChildren, LotShPtr rng) const;
        bool                            HasTabuSplit(TreeShPtr dest, TreeNode * destPolytomy, const std::vector<TreeNode *> & polytomyChildren, std::set<Split> & tabuSplits) const;
        double                          LnEdgeLenProbForSplit(const Split & s, const double b) const;
        void                            buildScratchTree();
        void                            DoTreeChecks(Tree & destTree, bool doRefreshesFirst, const char * tag) const;
		double                          countDistancesUsingBryantSteel(TreeNode *ff, unsigned n);
    };

typedef boost::shared_ptr<FocalTreeTopoProbCalculator> FocalTreeTopoProbCalculatorShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Computes topological priors used by BushMove to handle polytomous trees in MCMC analyses. Also provides several
|	utility functions for computing the number of tree topologies with varying degrees of resolution.
*/
class PolytomyTopoPriorCalculator : public TopoProbCalculator
	{
	public:
										PolytomyTopoPriorCalculator();
										~PolytomyTopoPriorCalculator();

		bool							IsResolutionClassPrior() const;
		bool							IsPolytomyPrior() const;

		bool							IsRooted() const;
		bool							IsUnrooted() const;

		void							SetNTax(unsigned n);
		unsigned						GetNTax() const;

		void							ChooseRooted();
		void							ChooseUnrooted();

		double							GetLnCount(unsigned n, unsigned m);
		double							GetLnSaturatedCount(unsigned n);
		double							GetLnTotalCount(unsigned n);
		std::vector<double>				GetLnCounts();

		std::vector<double>				GetCountsVect();
		std::vector<int>			    GetNFactorsVect();

		void							ChooseResolutionClassPrior();
		void							ChoosePolytomyPrior();

		void							SetC(double c);
		double							GetC() const;

		void							SetLnScalingFactor(double lnf);
		double							GetLnScalingFactor() const;

        virtual double                  GetLnTopoProb(TreeShPtr t);

		double							GetLnTopologyPrior(unsigned m);
		double							GetLnNormalizedTopologyPrior(unsigned m);
		double							GetLnNormConstant();

		std::vector<double>				GetTopoPriorVect();
		std::vector<double>				GetRealizedResClassPriorsVect();

        unsigned                        sample(LotShPtr rng);

		void							Reset();

	private:

		void							RecalcCountsAndPriorsImpl(unsigned n);
		void							RecalcPriorsImpl();

		unsigned						ntax;						/**< PolytomyTopoPriorCalculator is currently holding counts for this number of taxa */
		bool							is_rooted;					/**< If false, ntax is number of tips in an unrooted tree; if true, ntax is number of tips in a rooted tree */
		bool							is_resolution_class_prior;	/**< True for the resolution class prior, false for the polytomy prior */
		double							C;							/**< Determines strength of prior */

		bool							topo_priors_dirty;			/**< True if `ntax', `C', `is_rooted', or `is_resolution_class_prior' has changed since Reset was last called. Causes Reset to be called if any function is called that requires access to either `counts' or `polytomy_prior' */
		bool							counts_dirty;			    /**< True if counts need to be recalculated */

        double							log_scaling_factor;			/**< The natural log of the base scaling factor; the scaling factor raised to some power is removed from each count in `counts' for large numbers of taxa */
        std::vector<int>			    nfactors;				    /**< Holds, for each count in `counts', the number of times a factor `scaling_factor' needed to be removed from the count (count i is thus counts[i]*(scaling_factor^nfactors[i])) */
		std::vector<double> 			counts;						/**< Vector of length `ntax' holding counts for trees having all possible numbers of internal nodes for ntax taxa (note: counts might be scaled; see `scaling_factor' and the `nfactors' vector) */
        double							log_total_count;			/**< Holds the natural log of the sum of all counts (useful for normalizing) */
		std::vector<double>				topology_prior;				/**< Vector the mth element of which holds the unnormalized prior probability of a tree with m internal nodes; the 0th element holds the normalizing constant */
	};

typedef boost::shared_ptr<PolytomyTopoPriorCalculator> PolytomyTopoPriorCalculatorShPtr;

}	// namespace phycas

#endif
