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

#if ! defined(STATE_FREQ_MOVE_HPP)
#define STATE_FREQ_MOVE_HPP

#include <vector>                           // for std::vector
#include <boost/shared_ptr.hpp>             // for boost::shared_ptr
#include <boost/weak_ptr.hpp>               // for boost::weak_ptr
#include "mcmc_updater.hpp"		// for base class MCMCUpdater
#include "partition_model.hpp"   // for PartitionModelShPtr definition
#include "dirichlet_move.hpp"
#include "multivariate_probability_distribution.hpp"

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

typedef boost::shared_ptr<DirichletDistribution>    DirichletShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	A StateFreqMove proposes new state frequencies that are slightly different than the current frequencies by sampling
|   from a Dirichlet distribution with parameters equal to the current frequencies multiplied by a large value (the
|   tuning parameter 'psi').
*/
class StateFreqMove : public DirichletMove
	{
	public:
									StateFreqMove();
                                    virtual ~StateFreqMove() {}

		virtual void				sendCurrValuesToModel(const double_vect_t & v);
		virtual void				getCurrValuesFromModel(double_vect_t & v) const;
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams();
        virtual void                setParams(const std::vector<double> & v);

	private:

		StateFreqMove &				operator=(const StateFreqMove &);	// never use - don't define
    };

} // namespace phycas

#endif
