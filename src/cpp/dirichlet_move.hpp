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

#if ! defined(DIRICHLET_MOVE_HPP)
#define DIRICHLET_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include "mcmc_updater.hpp"		// for base class MCMCUpdater
#include "partition_model.hpp"		// for PartitionModelShPtr definition
#include "multivariate_probability_distribution.hpp"

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

typedef boost::shared_ptr<DirichletDistribution>    DirichletShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	A DirichletMove proposes new parameter values that are slightly different than the current parameter values by
|   sampling from a Dirichlet distribution with parameters equal to the current frequencies multiplied by a large
|   value (the tuning parameter 'psi').
*/
class DirichletMove : public MCMCUpdater
	{
	public:
									DirichletMove();
									virtual ~DirichletMove()
										{
										//std::cerr << "DirichletMove dying..." << std::endl;
										}

		// Accessors
		double						getTuningParameter() const;
		unsigned					getDimension() const;


		// Modifiers
		void						setTuningParameter(double x);
		void						setMaxPsi(double x);
		void						setMinPsi(double x);
		void						setDimension(unsigned d);

		// Utilities
		void						reset();
		virtual void				sendCurrValuesToModel(const double_vect_t & v) {}
		virtual void				getCurrValuesFromModel(double_vect_t & v) const {}
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams() {}
        virtual void                setParams(const std::vector<double> & v) {}

		// These are virtual functions in the MCMCUpdater base class
		virtual bool				update();
		virtual double				getLnHastingsRatio() const;
		virtual double				getLnJacobian() const;
		virtual void				proposeNewState();
		virtual void				revert();
		virtual void				accept();

	private:

		DirichletMove &				operator=(const DirichletMove &);	// never use - don't define

	protected:

		unsigned					dim;			/**< The number of parameters governed by this move */
		double						psi;			/**< Larger values result in changes of smaller magnitude */
		std::vector<double> 		new_params;	    /**< Proposed new parameter values */
		std::vector<double> 		orig_params;	/**< Saved parameter values (in case revert is necessary) */
		std::vector<double> 		c_forward;	    /**< Dirichlet parameter vector used to propose new frequencies */
		std::vector<double> 		c_reverse;	    /**< Dirichlet parameter vector used to propose original frequencies (only used to compute Hastings ratio) */
        DirichletShPtr              dir_forward;    /**< Points to an ad hoc Dirichlet distribution object used to assess the forward move density and to propose a new frequency vector */
        DirichletShPtr              dir_reverse;    /**< Points to an ad hoc Dirichlet distribution object used to assess the reverse move density */
	};

} // namespace phycas

#endif
