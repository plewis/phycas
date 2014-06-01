/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#if ! defined(MCMC_MOVE_HPP)
#define MCMC_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "phypy/phylogeny/tree_manip.hpp"			// for TreeManip

struct CIPRES_Matrix;

//class ProbabilityDistribution;

namespace phycas
{

//class TreeNode;

//class Tree;
//typedef boost::shared_ptr<Tree>					TreeShPtr;

//class Model;
//typedef boost::shared_ptr<Model>					ModelShPtr;

//class TreeLikelihood;
//typedef boost::shared_ptr<TreeLikelihood>			TreeLikeShPtr;

class TreeManip;
typedef boost::shared_ptr<TreeManip>				TreeManipShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

//class HKY;

/*----------------------------------------------------------------------------------------------------------------------
|	Base class for classes that encapsulate Metropolis proposals used for MCMC.
*/
class MCMCMove : public boost::enable_shared_from_this<MCMCMove>
	{
	public:
								MCMCMove();
								virtual	~MCMCMove()
									{
									//std::cerr << "MCMCMove dying..." << std::endl;
									}

		// Accessors
		const std::string &		getName() const;
		double					getLnLike() const;
		double					getLnPrior() const;

		// Modifiers
		virtual void			setName(const std::string & s);
		virtual void			setLot(LotShPtr p);
		virtual void			setTree(TreeShPtr p);
		virtual void			setTreeLikelihood(TreeLikeShPtr p);
		virtual void			setModel(ModelShPtr p);
		virtual void			setChainManager(ChainManagerWkPtr p);

		// Utilities
		virtual bool			update();

	protected:

		// Virtual functions that should be overridden in derived classes
		//
		virtual double			getLnHastingsRatio() const	= 0;
		virtual double			getLnJacobian() const		= 0;
		virtual void 			proposeNewState()			= 0;
		virtual void			accept()					= 0;
		virtual void			revert()					= 0;

	protected:

		std::string				name;				/**< Storage for the name of this parameter to be used in reporting to the user */
		TreeShPtr				tree;				/**< The tree on which the likelihood will be calculated */
		TreeManip				tree_manipulator;	/**< The object that facilitates topological rearrangements */
		ModelShPtr				model;				/**< The substitution model to be used in computing the likelihood */
		TreeLikeShPtr			likelihood;			/**< The object that knows how to calculate the likelihood */
		LotShPtr				rng;				/**< The pseudorandom number generator object used in updating parameter value */
		ChainManagerWkPtr		chain_mgr;			/**< The object that knows how to compute the joint log prior density */
		double					curr_ln_prior;		/**< The most recently calculated log prior (for this parameter only) */
		double					curr_ln_like;		/**< The most recently calculated log likelihood */
		double					ln_zero;			/**< The value to return for the posterior density if the parameter value is out of range */
	};

typedef boost::shared_ptr<MCMCMove>		MCMCMoveShPtr;
typedef std::vector<MCMCMoveShPtr>		MCMCMoveVect;
typedef MCMCMoveVect::iterator			MCMCMoveIter;

} // namespace phycas

#endif
