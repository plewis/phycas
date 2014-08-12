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

#if ! defined(EDGE_MOVE_HPP)
#define EDGE_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "mcmc_updater.hpp"		// for base class MCMCUpdater

//class ExponentialDistribution;
//typedef boost::shared_ptr<ExponentialDistribution>	ExponentialDistributionShPtr;

namespace phycas
{


class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

//typedef std::map<unsigned, std::vector<double> > PolytomyDistrMap;
//typedef std::vector<double> VecPolytomyDistr;

/*----------------------------------------------------------------------------------------------------------------------
|	An EdgeMove changes the length of just one randomly-chosen edge in the tree. An edge chosen at random is set to the
|	value Y = m*exp(`lambda'*(u - 0.5)), where m is the original length and u is a Uniform(0,1) random deviate. Under
|	this proposal scheme, the random variable Y has the following properties:
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
class EdgeMove : public MCMCUpdater
	{
	public:
									EdgeMove();
									virtual ~EdgeMove()
										{
										//std::cerr << "EdgeMove dying..." << std::endl;
										}


		// Accessors
		//
		double						getTuningParameter() const;
		//bool						addEdgeMoveProposed() const;

		// Modifiers
		//
		void						setTuningParameter(double x);
		//void						setEdgeLenDistMean(double mean);

		// Utilities
		//
		void						reset();
		//void						finalize();

		// These are virtual functions in the MCMCUpdater base class
		//
        //POLTMP2 virtual void                setPosteriorTuningParam(double x);
        //POLTMP2 virtual void                setPriorTuningParam(double x);
		//POLTMP2 virtual void	            setBoldness(double x);
		virtual bool				update();
		virtual double				getLnHastingsRatio() const;
		virtual double				getLnJacobian() const;
		virtual void				proposeNewState();
		virtual void				revert();
		virtual void				accept();

	private:

		EdgeMove &					operator=(const EdgeMove &);	// never use - don't define

	private:

		//POLTMP2 double			            boldness;		/**< Ranges from 0 to 100 and determines the boldness of the move */
        double                      lambda;         /**< the tuning parameter used for this move */
        //POLTMP2 double                      min_lambda;     /**< the tuning parameter used for exploring the posterior distribution */
        //POLTMP2 double                      max_lambda;     /**< the tuning parameter used for exploring the prior distribution */

		double						origEdgelen;	/**< Length of modified edge saved (in case revert is necessary) */
		TreeNode *					origNode;		/**< Node owning the modified edge (in case revert is necessary) */
		TreeNode *					likeRoot;		/**< Node to be used as the likelihood root (equals origNode if origNode is internal, otherwise equals origNode's parent) */

		std::vector<double>		    one_edgelen;						/**< Workspace declared here to avoid unnecessary allocs/deallocs */
	};

} // namespace phycas

#endif
