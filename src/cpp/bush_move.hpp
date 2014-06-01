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

#if ! defined(BUSH_MOVE_HPP)
#define BUSH_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include "mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class ExponentialDistribution;
typedef boost::shared_ptr<ExponentialDistribution>	ExponentialDistributionShPtr;

class PolytomyTopoPriorCalculator;
typedef boost::shared_ptr<PolytomyTopoPriorCalculator>		PolytomyTopoPriorCalculatorShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

typedef std::map<unsigned, std::vector<double> > PolytomyDistrMap;
typedef std::vector<double> VecPolytomyDistr;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the add-edge and delete-edge proposal described in Lewis, P., M. T. Holder and K. E. Holsinger. 2005.
|   Polytomies and Bayesian phylogenetic inference. Systematic Biology 54(2):241-253.
*/
class BushMove : public MCMCUpdater
	{
	public:
                                            BushMove();
                                            virtual ~BushMove();

		// Accessors
		//
		bool                                addEdgeMoveProposed() const;
		virtual bool                        computesTopologyPrior() const;

		// Modifiers
		//
		void                                setEdgeLenDistMean(double mean);
		PolytomyTopoPriorCalculatorShPtr    getPolytomyTopoPriorCalculator();

		// Utilities
		//
		void                                finalize();

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual bool                        update();
		virtual double                      getLnHastingsRatio() const;
		virtual double                      getLnJacobian() const;
		virtual void                        proposeNewState();
		virtual void                        revert();
		virtual void                        accept();

	private:

		void                                reset();
		void                                proposeAddEdgeMove(TreeNode * nd);
		void                                proposeDeleteEdgeMove(TreeNode * nd);
		const VecPolytomyDistr &            computePolytomyDistribution(unsigned nspokes);
		void                                refreshPolytomies();
		BushMove &                          operator=(const BushMove &);	// never use - don't define

	private:

		unsigned                            num_taxa;							/**< The number of taxa */
		unsigned                            num_nodes_in_fully_resolved_tree;	/**< Total number of nodes in a fully resolved tree having the same number of taxa, N (= 2N - 2) */
		bool                                add_edge_move_proposed;				/**< True if last move proposed was an add-edge move; false otherwise */

		double                              edgelen_mean;						/**< The mean of the exponential distribution governing new edge lengths */
		ExponentialDistributionShPtr        edgelen_dist;						/**< The distribution from which new edge lengths are chosen */

		PolytomyDistrMap                    poly_prob;							/**< Map of probability distributions for x where x is the number of spokes split off in an add-edge move that divides the spokes of a polytomy of size n into two groups. A polytomy with n spokes is split into two nodes, with x spokes going to one node and n - x remaining with the other node. polyProb[n][x] is prob(x|n). */

		std::vector<TreeNode *>             polytomies;							/**< Vector containing pointers to TreeNode objects representing polytomies */

		double                              new_edgelen;						/**< Stores new edge length created by an add-edge move */
		unsigned                            polytomy_size;						/**< Stores number of spokes in polytomy broken by last proposed add-edge move */
		unsigned                            num_polytomies;						/**< Stores number of polytomies in tree */

		double                              orig_edgelen;						/**< Length of deleted edge saved (in case revert of delete-edge move is necessary) */
		TreeNode *                          orig_par;							/**< Parent of deleted node (in case revert of delete-edge move is necessary) */
		TreeNode *                          orig_lchild;						/**< Leftmost child of deleted node (in case revert of delete-edge move is necessary) */
		TreeNode *                          orig_rchild;						/**< Rightmost child of deleted node (in case revert of delete-edge move is necessary) */
		double                              ln_jacobian;						/**< The natural log of the Jacobian for the move last proposed */
		double                              ln_hastings;						/**< The natural log of the Hastings ratio for the move last proposed */

		CDF                                 cdf;								/**< CDF object needed for its LnGamma function */

	};

} // namespace phycas

#endif
