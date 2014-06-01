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

//#include "phycas/force_include.h"
#include "basic_tree.hpp"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "tip_data.hpp"
#include "cond_likelihood.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will only be used for simulating data and thus needs only to allocate the
|	`state' and `pMatrixTranspose' data members.
*/
TipData::TipData(
  PartitionModelShPtr           partition,      /**< is the PartitionModel object containing information about the number of states and rates for each subset */
  CondLikelihoodStorageShPtr    cla_storage)
	:
	state(-1),
	cla_pool(cla_storage)
	{
	const unsigned num_subsets = partition->getNumSubsets();

	pMatrixTranspose.resize(num_subsets);

	for (unsigned i = 0; i < num_subsets; ++i)
		{
		// allocate memory for the ScopedThreeDMatrix of transition probabilities for subset i
		const unsigned num_states		= partition->subset_num_states[i];
		const unsigned num_obs_states	= partition->subset_num_states[i];
		const unsigned num_rates		= partition->subset_num_rates[i];
		pMatrixTranspose[i].Initialize(num_rates, num_obs_states, num_states);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will be used in likelihood calculations and thus needs to store the observed
|	data for the tip as well as information about local state codes.
*/
TipData::TipData(
  PartitionModelShPtr							partition,			/**< is the PartitionModel object containing information about the number of states and rates for each subset */
  const state_list_pos_vect_t &					positions,			/**< is the vector of vectors of positions of each state into the `state_codes' vector */
  const state_list_vect_t &						states,				/**< is the vector of vectors of `state_codes' */
  CondLikelihoodStorageShPtr 					cla_storage)		/**< is the pool of available conditional likelihood arrays */
	:
	state(-1),
	cla_pool(cla_storage)
	{
	const unsigned num_subsets = partition->getNumSubsets();
	state_list_pos	= state_list_pos_vect_t(num_subsets);
	state_codes 	= state_list_vect_t(num_subsets);

	pMatrixTranspose.resize(num_subsets);

	for (unsigned i = 0; i < num_subsets; ++i)
		{
		state_list_t & subset_state_codes = state_codes[i];
		// copy state list positions for subset i
		state_list_pos[i].resize(positions[i].size());
		std::copy(positions[i].begin(), positions[i].end(), state_list_pos[i].begin());

		// copy state codes for subset i
		subset_state_codes.resize(states[i].size());
		std::copy(states[i].begin(), states[i].end(), subset_state_codes.begin());

		// allocate memory for the ScopedThreeDMatrix of transition probabilities for subset i
		const unsigned num_obs_states	= partition->subset_num_states[i] + 1 + positions[i].size();
		const unsigned num_rates		= partition->subset_num_rates[i];
		const unsigned num_states		= partition->subset_num_states[i];
		//const unsigned num_patterns = partition->subset_num_patterns[i];
		pMatrixTranspose[i].Initialize(num_rates, num_obs_states, num_states);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor ensures that all CLA structures are returned to `cla_pool'. The `ownedPMatrices' and `univents' data
|   members takes care of deleting themselves when they go out of scope.
*/
TipData::~TipData()
	{
	// Invalidate the parental CLAs if they exist
	if (parWorkingCLA)
		{
		cla_pool->putCondLikelihood(parWorkingCLA);
		parWorkingCLA.reset();
		}

	// Remove cached parental CLAs if they exist
	if (parCachedCLA)
		{
		cla_pool->putCondLikelihood(parCachedCLA);
		parCachedCLA.reset();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member. If `parWorkingCLA' does not currently point to anything, a CondLikelihood
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
CondLikelihoodShPtr TipData::getParentalCondLikePtr()
	{
	if (!parWorkingCLA)
		parWorkingCLA = cla_pool->getCondLikelihood();
	return parWorkingCLA;
	}

}	// namespace phycas
