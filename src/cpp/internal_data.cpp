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

#include "cond_likelihood.hpp"
#include "cond_likelihood_storage.hpp"
#include "internal_data.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor initializes the conditional likelihood vectors using `nRates' and `nStates', and initializes the
|	`pMatrices' data member using `pMat'. The parameter `managePMatrices' determines whether to allocate space for the
|	transition probability matrices.
*/
InternalData::InternalData(
  PartitionModelShPtr	partition,			/**< is the PartitionModel object containing information about the number of states and rates for each subset */
  CondLikelihoodStorageShPtr cla_storage)
	:
	state(-1),
	cla_pool(cla_storage)
	{
	const unsigned num_subsets = partition->getNumSubsets();

	pMatrices.resize(num_subsets);
	for (unsigned i = 0; i < num_subsets; ++i)
		{
		const unsigned num_rates	= partition->subset_num_rates[i];
		const unsigned num_states	= partition->subset_num_states[i];
		pMatrices[i].Initialize(num_rates, num_states, num_states);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor ensures that all CLA structures are returned to `cla_pool'. The `ownedPMatrices' take care of deleting
|   themselves when they go out of scope.
*/
InternalData::~InternalData()
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

	// Invalidate the filial CLAs if they exist
	if (childWorkingCLA)
		{
		cla_pool->putCondLikelihood(childWorkingCLA);
		childWorkingCLA.reset();
		}

    // Remove cached filial CLAs if they exist
	if (childCachedCLA)
		{
		cla_pool->putCondLikelihood(childCachedCLA);
		childCachedCLA.reset();
		}
	}


}	// namespace phycas
