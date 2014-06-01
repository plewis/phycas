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

#if ! defined(LARGET_SIMON_MOVE_HPP)
#define LARGET_SIMON_MOVE_HPP

#include <vector>							// for std::vector
#include <boost/shared_ptr.hpp>				// for boost::shared_ptr
#include <boost/weak_ptr.hpp>				// for boost::weak_ptr
#include "mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>	ChainManagerWkPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the Larget-Simon local move. The LS move is the default move, but an exception must be made in the case
|	of a star tree, which does not have three contiguous edges required for the LS move. For star trees, the proposal
|	changes the length of just one randomly-chosen edge in the tree. An edge chosen at random is set to the value
|	Y = m*exp(`lambda'*(u - 0.5)), where m is the original length and u is a Uniform(0,1) random deviate. Under this
|	proposal scheme, the random variable Y has the following properties:
|>
|	density         f(Y) = 1/(lambda*Y)
|	cdf             F(Y) = 0.5 + (1/lambda) log(Y/m)
|	minimum         m*exp(-lambda/2)
|	maximum         m*exp(lambda/2)
|	mean            (m/lambda)[exp(lambda/2) - exp(-lambda/2)]
|	variance        (m/lambda)^2 [(lambda/2)(exp(lambda) - exp(-lambda)) - (exp(lambda/2) - exp(-lambda/2))^2]
|>
|	With a starting edge length of 1.0, the proposed edge lengths have increasing mean and variance with increasing
|	lambda, but values of lambda in the range 0.1 to 2.0 appear to be reasonable.
|>
|	lambda       mean        s.d.
|	-----------------------------
|	 0.1      1.00042     0.02888
|	 0.5      1.01045     0.14554
|	 1.0      1.04219     0.29840
|	 2.0      1.17520     0.65752
|	10.0     14.84064    29.68297
|>
*/
class LargetSimonMove : public MCMCUpdater
	{
	public:
						LargetSimonMove();
						virtual ~LargetSimonMove();

		unsigned		getWhichCase() const;
		void			setLambda(double x);
		double			getLambda() const;
		bool			topologyChanged() const;
		void			defaultProposeNewState();
		void			starTreeProposeNewState();
        TreeNode *      randomInternalAboveSubroot();
        TreeNode *      randomChild(TreeNode * nd);
        TreeNode *      chooseZ(TreeNode * middle);

		virtual double  sampleWorkingPrior() const;

		virtual bool    computesTopologyPrior() const {return true;}
		// These are virtual functions in the MCMCUpdater base class
		//
		double			recalcWorkingPrior() const;
		bool 			isPriorSteward() const;
        virtual void    setPosteriorTuningParam(double x);
        virtual void    setPriorTuningParam(double x);
		virtual void	setBoldness(double x);
		virtual bool	update();
		virtual double	getLnHastingsRatio() const;
		virtual double	getLnJacobian() const;
		virtual void	proposeNewState();
		virtual void	revert();
		virtual void	accept();

	private:

		double			    boldness;		            /**< Ranges from 0 to 100 and determines the boldness of the move */
        double              min_lambda;                 /**< The tuning parameter used for exploring the posterior distribution */
        double              max_lambda;                 /**< The tuning parameter used for exploring the prior distribution */
		double			    lambda;						/**< The tuning parameter for this move (the factor used in modifying backbone length) */

		//TreeNode *		    ndBase;						/**< Most ancestral node involved in the move, used as the center of the likelihood calcuation (and in revert) */
		TreeNode *		    ndX;						/**< Node at one end of segment involved in move; used by Revert to undo a move */
		TreeNode *		    ndY;						/**< One of two nodes in the middle of segment involved in move; used by Revert to undo a move */
		TreeNode *		    ndZ;						/**< Node at other end (from ndX) of segment involved in move; used by Revert to undo a move */
		double			    x;						    /**< Original length of ndX's branch; used by Revert to undo a move */
		double			    y;						    /**< Original length of ndX's branch; used by Revert to undo a move */
		double			    z;						    /**< Original length of ndX's branch; used by Revert to undo a move */
		TreeNode *		    swap1;						/**< First of the two nodes involved in an NNI swap; NULL if no swap was performed; used by Revert to undo a move */
		TreeNode *		    swap2;						/**< Second of the two nodes involved in an NNI swap; NULL if no swap was performed; used by Revert to undo a move */

		double			    m;							/**< Original 3-segment length; needed for computing Hastings ratio */
		double			    mstar;						/**< Modified 3-segment length; needed for computing Hastings ratio */

        double              expand_contract_factor;     /**< The factor by which the selected 3-edge segment is expanded or contracted */

		bool			    topol_changed;				/**< If true, last proposal changed topology */
		unsigned		    which_case;					/**< Which of the eight possible cases was tried last */

		std::vector<double> three_edgelens;			/**< workspace declared here to avoid unnecessary allocs/deallocs */

		void			    reset();					/**< Returns variables involved with reversing a proposed move to the state needed for the start of another proposal */

		bool			    star_tree_proposal;			/**< True if last proposed move was on a star tree (only one randomly-chosen edge changed); False if last proposed move was not on a star tree */

		// These are needed for the star tree exception
		double				orig_edge_len;	/**< Length of modified edge saved (in case revert is necessary) */
		TreeNode *			orig_node;		/**< Node owning the modified edge (in case revert is necessary) */
		std::vector<double>	one_edgelen;	/**< workspace declared here to avoid unnecessary allocs/deallocs */
	};

} // namespace phycas

#endif
