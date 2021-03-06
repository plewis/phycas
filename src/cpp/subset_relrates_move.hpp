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

#if ! defined(SUBSET_RELRATES_MOVE_HPP)
#define SUBSET_RELRATES_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include "mcmc_updater.hpp"		// for base class MCMCUpdater
#include "partition_model.hpp"		// for PartitionModelShPtr definition
#include "dirichlet_move.hpp"

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

typedef boost::shared_ptr<DirichletDistribution>    DirichletShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	A SubsetRelRatesMove proposes new subset relative rates that are slightly different than the current rates by
|   sampling from a Dirichlet distribution with parameters equal to the current rates multiplied by a large value (the
|   tuning parameter 'psi').
*/
class SubsetRelRatesMove : public DirichletMove
	{
	public:
									SubsetRelRatesMove();
                                    virtual ~SubsetRelRatesMove() {}

		virtual void				sendCurrValuesToModel(const double_vect_t & v);
		virtual void				getCurrValuesFromModel(double_vect_t & v) const;
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams();
        virtual void                setParams(const double_vect_t & v);
		void						setPartitionModel(PartitionModelShPtr m);
		//void						setSubsetProportions(const double_vect_t & proportions);
		virtual void				proposeNewState();
		virtual bool				update();
		virtual void				accept();
		virtual void 				revert();

        std::string                 debugShowTransformedRelativeRates(const std::string msg) const;

	private:

		double_vect_t 				new_relrates;	    /**< Proposed new relative rate values */
		double_vect_t 				orig_relrates;		/**< Saved relative rate values (in case revert is necessary) */
		//double_vect_t				p;					/**< Vector of proportion of sites in each subset */
		PartitionModelShPtr			partition_model;	/**< Pointer to the partition model */
		SubsetRelRatesMove &		operator=(const SubsetRelRatesMove &);	// never use - don't define
    };

} // namespace phycas

#endif
