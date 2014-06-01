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

#if ! defined(TREE_LIKELIHOOD_HPP)
#define TREE_LIKELIHOOD_HPP

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "ncl/nxscxxdiscretematrix.h"

#include "states_patterns.hpp"
#include "model.hpp"
#include "cond_likelihood.hpp"
#include "cond_likelihood_storage.hpp"
#include "underflow_manager.hpp"
#include "partition_model.hpp"

namespace phycas
{
class CharSuperMatrix;
unsigned ** getNodeSMat(TreeNode * nd, unsigned subsetIndex);

typedef const double * const * const * ConstPMatrices;
typedef std::vector<unsigned int> StateListPos;
class CondLikelihood;
class Tree;
class TreeLikelihood;
template<typename T> class GenericEdgeEndpoints;
typedef GenericEdgeEndpoints<TreeNode *> EdgeEndpoints;
typedef GenericEdgeEndpoints<const TreeNode *> ConstEdgeEndpoints;

class TipData;
class InternalData;

class SimData;
typedef boost::shared_ptr<SimData>	SimDataShPtr;

class Lot;
typedef boost::shared_ptr<Lot>	LotShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Used for computing the likelihood on a tree.
*/
class TreeLikelihood
	{
	public:

										TreeLikelihood(PartitionModelShPtr mod);
								virtual ~TreeLikelihood(); //needed as long as startTreeViewer is virtual

		// Accessors
		bool                            getNoData() const {return no_data;}
		unsigned						getNTaxa() const;
		unsigned						getNChar() const;
		unsigned						getNumPatterns() const;
		unsigned                        getNumSubsets() const
			{
			if (!partition_model)
				return 1;
			return partition_model->getNumSubsets();
			}

		PartitionModelShPtr				getPartitionModel() const;
		const double_vect_t &			getRateMeans(unsigned i) const;
		const double_vect_t &			getRateProbs(unsigned i) const;
		unsigned						getNRatesTotal(unsigned i) const;
		double_vect_t					getCategoryLowerBoundaries(unsigned i) const;
		unsigned						getNumStates(unsigned i) const;
		const state_list_vect_t &		getStateList() const;
		const state_list_pos_vect_t &	getStateListPos() const;
		void							replacePartitionModel(PartitionModelShPtr);
		const count_vect_t &			getPatternCounts() const;
		void							recalcRelativeRates();
		const std::vector<unsigned> &	getListOfAllMissingSites() const;
		const std::vector<double> &		getSiteLikelihoods() const;
		const std::vector<double> &		getSiteUF() const;
		bool							storingSiteLikelihoods() const;
		const std::vector<unsigned> &	getCharIndexToPatternIndex() const;

		// Modifiers
		void							setNoData();
		void							setHaveData();
		void							storeSiteLikelihoods(bool yes);

		// Utilities
		unsigned						sumPatternCounts() const;
		void							releasePartitionModel();
		double							calcLogLikeAtSubstitutionSaturation() const;

		void							addOrphanTip(TreeShPtr t, unsigned row, std::string name);
		void							addDecoratedInternalNode(TreeShPtr t, unsigned num = UINT_MAX);
		void							prepareForSimulation(TreeShPtr);
		void							prepareForLikelihood(TreeShPtr);
		void							prepareInternalNodeForLikelihood(TreeNode * nd);

		void							patternMapToVect(const pattern_map_vect_t & pattern_map_vect, pattern_to_sites_map_t & pattern_to_sites_map);
		void							copyDataFromDiscreteMatrix(const CharSuperMatrix *, const std::vector<unsigned> & partition_info);
		void							copyDataFromSimData(SimDataShPtr sim_data);

		bool							invalidateNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root);
		bool							invalidateBothEnds(TreeNode * ref_nd, TreeNode * unused = NULL);
		bool							invalidateBothEndsDiscardCache(TreeNode * ref_nd, TreeNode * unused = NULL);
		void							invalidateAwayFromNode(TreeNode & focalNode);

		bool							restoreFromCacheNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root);
		bool							restoreFromCacheBothEnds(TreeNode * ref_nd, TreeNode * unused = NULL);
		bool							restoreFromCacheParentalOnly(TreeNode * ref_nd, TreeNode * unused = NULL);
		bool							discardCacheBothEnds(TreeNode * ref_nd, TreeNode * unused = NULL);
		void							restoreFromCacheAwayFromNode(TreeNode & focalNode);
		void							discardCacheAwayFromNode(TreeNode & focalNode);

        void							debugCompressedDataInfo(std::string filename);
        void                            debugUncompressedDataInfo(std::string filename);

		const CondLikelihoodStorageShPtr	getCLAStorage() const;

		unsigned						bytesPerCLA() const;
		unsigned						numCLAsCreated() const;
		unsigned						numCLAsStored() const;

		TreeNode *						storeAllCLAs(TreeShPtr t);
		bool							debugCheckCLAsRemainInTree(TreeShPtr t) const;
		bool							debugCheckForUncachedCLAs(TreeShPtr t) const;
		bool							debugCheckCLAsRemainInNode(TreeNode * nd) const;

		bool							isValid(const TreeNode *focal, const TreeNode *avoidNd);
		void							refreshCLA(TreeNode & nd, const TreeNode * avoid);
		double							calcLnLFromNode(TreeNode & focal_node, TreeShPtr t);
		double							calcLnL(TreeShPtr);

		std::string						listPatterns(bool translate);
		std::string						getStateStr(unsigned i, state_code_t state) const;

		void							setDebug(bool on);

		void							simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);
		void							simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);

		void							swapInternalDataAndEdgeLen(TreeNode * nd1, TreeNode * nd2);
		void							addDataTo(SimData & other);

		unsigned						getNumLikelihoodEvals() const;
		void							incrementNumLikelihoodEvals();
		void							resetNumLikelihoodEvals();

		TreeNode *						getLikelihoodRoot();
		void							useAsLikelihoodRoot(TreeNode * nd);

		int								getLikelihoodRootNodeNum() const;
		void							debugSaveCLAs(TreeShPtr t, std::string fn, bool overwrite);
		virtual int						startTreeViewer(TreeShPtr, std::string, unsigned site = 0) const {return 0;}

		void							setUFNumEdges(unsigned nedges);

		void							calcPMatTranspose(unsigned i, double * * * transPMats, const uint_vect_t & stateListPosVec, double edgeLength);
		void							calcPMat(unsigned i, double * * * p, double edgeLength); //

		void							calcCLATwoTips(CondLikelihood & condLike, const TipData & leftTip, const TipData & rightTip);
		void							calcCLAOneTip(CondLikelihood & condLike, const TipData & leftChild, const InternalData & rightChild, const CondLikelihood & rightCondLike);
		void							calcCLANoTips(CondLikelihood & condLike, const InternalData & leftChild, const CondLikelihood & leftCondLike, const InternalData & rightChild, const CondLikelihood & rightCondLike);

		void							conditionOnAdditionalTip(CondLikelihood & condLike, const TipData & tipData);
		void							conditionOnAdditionalInternal(CondLikelihood & condLike, const InternalData & child, const CondLikelihood & childCondLike);

		double							harvestLnL(EdgeEndpoints & focalEdge, TreeShPtr t);
		double							harvestLnLFromValidEdge(ConstEdgeEndpoints & focalEdge);
		double							harvestLnLFromValidNode(TreeNode *focalNode);
        void                            debugWalkTreeShowCondLikes(TreeShPtr t);

		unsigned						buildConstantStatesVector();

		CondLikelihoodStorageShPtr		getCondLikelihoodStorage();

		void							slideNode(double fraction, TreeNode * slider, TreeNode * other);

	protected:

		UnderflowManager				underflow_manager;		/**< The object that takes care of underflow correction when computing likelihood for large trees */

		TreeNode *						likelihood_root;		/**< If not NULL< calcLnL will use this node as the likelihood root, then reset it to NULL before returning */
		CondLikelihoodStorageShPtr		cla_pool;

		bool							store_site_likes;		/**< If true, calcLnL always stores the site likelihoods in the `site_likelihood' data member; if false, the `site_likelihood' data member is not updated by calcLnL */
		bool							no_data;				/**< If true, calcLnL always returns 0.0 (useful for allowing MCMC to explore the prior) */

		PartitionModelShPtr				partition_model;		/**< The object that holds information about the model applied to each subset of the data partition */

		//@POL these next four should logically be inside the PartitionModel class
		double_vect_vect_t				rate_means;				/**< Vector of relative rate vectors (one rate vector for each partition subset) */
		double_vect_vect_t				rate_probs;				/**< Vector of relative rate probability vectors  (one probability vector for each partition subset) */
		state_list_vect_t				state_list;				/**< `state_list[i]' provides a vector of state code definitions for subset i */
		state_list_pos_vect_t			state_list_pos;			/**< `state_list_pos[i]' is a vector of positions of states in `state_list[i]' */

		bool							debugging_now;			/**< For debugging, indicates whether user wants to see debugging output */

	protected:

		// Utilities
		void							buildPatternReprVector(std::vector<std::string> &, TreeShPtr);
		TipData *						allocateTipData(unsigned);
		InternalData *					allocateInternalData();

		void	 						storePattern(pattern_map_t & pattern_map, pattern_to_sites_map_t & pattern_to_site_map, const std::vector<int8_t> & pattern, const unsigned pattern_index, const pattern_count_t weight, bool codon_model);
		unsigned						compressDataMatrix(const NxsCXXDiscreteMatrix &, const std::vector<unsigned> & partition_info);
		void							calcPMatCommon(unsigned i, double * * * pMatrices, double edgeLength);

		void							calcTMatForSim(unsigned i, TipData &, double);
		void							simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs);

	public:

        // Tree Length Prior
        void                            setTreeLengthPrior(TreeLengthDistributionShPtr tl_prior);
        TreeLengthDistributionShPtr     getTreeLengthPrior();
        void                            setTreeLengthRefDist(TreeLengthDistributionShPtr tl_ref_dist);
        TreeLengthDistributionShPtr     getTreeLengthRefDist();
        //void                            educateTreeLenWorkingPrior(TreeShPtr tree, double edge_length, unsigned edge_len_index);
        //void                            finalizeTreeLenWorkingPrior(LotShPtr lot);

	protected:

		unsigned						nevals;					/**< For debugging, records the number of times the likelihood is calculated */

	public: //@POL these should be protected rather than public

		unsigned						_num_taxa;                  /**< The total number of taxa */
		unsigned						_num_char;                  /**< The total number of characters (or sites), value assigned in copyDataFromDiscreteMatrix */

        double_vect_vect_t              tl_fitting_sample;          /**< storage used for fitting the tree length prior */
        TreeLengthDistributionShPtr     _tl_prior;                  /**< the Rannala-Zhu-Yang tree length prior (null if not being used) */
        TreeLengthDistributionShPtr     _tl_ref_dist;               /**< the distribution serving as the reference distribution when using a tree length prior */

		count_vect_t					pattern_counts;				/**< vector of pattern counts */
		uint_vect_t						subset_offset;				/**< `subset_offset'[i] holds the index into `pattern_vect' where the patterns from subset i begin. The number of elements is one greater than the number of subsets (the last element holds num_patterns to make it easy to find the end of any subset, including the last subset). */
		pattern_vect_t					pattern_vect;				/**< vector of patterns (all patterns for a given partition subset are contiguous, but within a subset, order may differ from data file order) */
		double_vect_t					site_likelihood;			/**< site_likelihood[pat] stores the site likelihood for pattern pat, but only if `store_site_likes' is true */
		pattern_to_sites_t				pattern_to_sites;			/**< vector of lists that provides a list of character indices for each pattern in `pattern_vect'. For example, if pattern j is found at sites 0, 15, and 167, then pattern_to_sites[j] is the list [0, 15, 167] */
		uint_vect_t						charIndexToPatternIndex; 	/**< maps original character index to the position of the corresponding element in `pattern_vect' */
		uint_vect_t						constant_states;			/**< keeps track of the states for potentially constant sites. See TreeLikelihood::buildConstantStatesVector for description of the structure of this vector. */
		uint_vect_t						all_missing;				/**< keeps track of sites excluded automatically because they have missing data for all taxa. */
		double_vect_t					site_uf;					/**< site_uf[pat] stores the underflow correction factor used for pattern pat, but only if `store_site_likes' is true */
	};

/// used to get access to a CLA to write it
CondLikelihoodShPtr getCondLikePtr(EdgeEndpoints edge);
CondLikelihoodShPtr getCondLikePtr(TreeNode * focal_nd, TreeNode * avoid);

/// read-only access to a CLA (must be valid already)
ConstCondLikelihoodShPtr getValidCondLikePtr(ConstEdgeEndpoints edge);
ConstCondLikelihoodShPtr getValidCondLikePtr(const TreeNode * focal_nd, const TreeNode * avoid);

} // namespace phycas

//#include "tree_likelihood.inl"

#endif
