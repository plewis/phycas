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

#include <numeric>
#include "cond_likelihood.hpp"
#include "cond_likelihood_storage.hpp"

namespace phycas
{

static unsigned next_cond_like_storage = 0;

/*----------------------------------------------------------------------------------------------------------------------
|	Ensures that the `cl_stack' contains at least `capacity' CondLikelihoodShPtr objects.
*/
void CondLikelihoodStorage::fillTo(unsigned capacity)
	{
	PHYCAS_ASSERT(std::accumulate(num_patterns.begin(), num_patterns.end(), 0) > 0);
	PHYCAS_ASSERT(std::accumulate(num_rates.begin(), num_rates.end(), 0) > 0);
	PHYCAS_ASSERT(std::accumulate(num_states.begin(), num_states.end(), 0) > 0);
	unsigned curr_sz = (unsigned)cl_stack.size();
	unsigned num_needed = (capacity > curr_sz ? capacity - curr_sz : 0);
	for (unsigned i = 0; i < num_needed; ++i)
		{
		cl_stack.push(CondLikelihoodShPtr(new CondLikelihood(num_patterns, num_rates, num_states)));
		num_created++;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `realloc_min' to 1 and `num_patterns', `num_rates' and `num_states' all to 0.
*/
CondLikelihoodStorage::CondLikelihoodStorage()
  :
  num_created(0),
  realloc_min(1)
	{
	which = next_cond_like_storage++;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor calls clearStack to destroy all CondLikelihoodShPtr objects currently stored.
*/
CondLikelihoodStorage::~CondLikelihoodStorage()
	{
	clearStack();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of elements in `cl_stack' (which equals the number of stored conditional likelihood arrays).
*/
unsigned CondLikelihoodStorage::getCLAPoolSize() const
	{
	return (unsigned)cl_stack.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `num_patterns', which holds the number of patterns used in determining the size
|	of newly-created CondLikelihood objects.
*/
unsigned CondLikelihoodStorage::getNumPatterns(
  unsigned i) const	/**< is the subset */
	{
	PHYCAS_ASSERT(num_patterns.size() > i);
	return num_patterns[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `num_rates', which holds the number of rates used in determining the size of
|	newly-created CondLikelihood objects.
*/
unsigned CondLikelihoodStorage::getNumRates(
  unsigned i) const	/**< is the subset */
	{
	PHYCAS_ASSERT(num_rates.size() > i);
	return num_rates[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `num_states', which holds the number of states used in determining the size of
|	newly-created CondLikelihood objects.
*/
unsigned CondLikelihoodStorage::getNumStates(
  unsigned i) const	/**< is the subset */
	{
	PHYCAS_ASSERT(num_states.size() > i);
	return num_states[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of bytes allocated for each CLA. This equals sizeof(LikeFltType) times the product of the number of
|	patterns, number of rates and number of states, summed over all partition subsets.
*/
unsigned CondLikelihoodStorage::bytesPerCLA() const
	{
	unsigned sz = num_patterns.size();
	PHYCAS_ASSERT(num_rates.size() == sz);
	PHYCAS_ASSERT(num_states.size() == sz);

	unsigned total = 0;
	for (unsigned i = 0; i < sz; ++i)
		{
		total += (unsigned)(num_patterns[i]*num_rates[i]*num_states[i]*sizeof(LikeFltType));
		}
	return total;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of the data member `num_created', which keeps track of the number of CondLikelihood objects created
|	since this object was constructed, or since the last call to clearStack, which resets `num_created' to zero.
*/
unsigned CondLikelihoodStorage::numCLAsCreated() const
	{
	return num_created;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current number of elements stored in the data member `cl_stack'. The total number of CLAs currently
|	checked out to the tree can be obtained as `num_created' minus the value returned by this function.
*/
unsigned CondLikelihoodStorage::numCLAsStored() const
	{
	return (unsigned)cl_stack.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `realloc_min' to 1, `cl_len' to `cond_like_len', and then calls the CondLikelihoodStorage::fillTo
|	member function to add `starting_size' new objects to the stack.
*/
//CondLikelihoodStorage::CondLikelihoodStorage(unsigned cond_like_len, unsigned starting_size)
// :
//  cl_len(cond_like_len),
//  realloc_min(1)
//	{
//	PHYCAS_ASSERT(cl_len > 0);
//	if (starting_size > 0)
//		fillTo(starting_size);
//	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `cl_stack' is empty, calls CondLikelihoodStorage::fillTo to add `realloc_min' more objects to the stack. After
|	ensuring that the `cl_stack' is not empty, pops a CondLikelihoodShPtr off and returns it.
*/
CondLikelihoodShPtr CondLikelihoodStorage::getCondLikelihood()
	{
	if (cl_stack.empty())
		fillTo(realloc_min);

	CondLikelihoodShPtr cl_ptr = cl_stack.top();
	cl_stack.pop();

	return cl_ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Pushes supplied CondLikelihoodShPtr `p' onto `cl_stack'.
*/
void CondLikelihoodStorage::putCondLikelihood(CondLikelihoodShPtr p)
	{
	cl_stack.push(p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `realloc_min', which determines how many CondLikelihood objects are to be created when the
|	`cl_stack' is discovered to be empty.
*/
void CondLikelihoodStorage::setReallocMin(unsigned sz)
	{
	PHYCAS_ASSERT(sz > 0);
	realloc_min = sz;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all CondLikelihoodShPtr objects currently in the `cl_stack'.
*/
void CondLikelihoodStorage::clearStack()
	{
	while (!cl_stack.empty())
		cl_stack.pop();
	num_created = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data members `num_patterns', `num_rates' and `num_states', which determine the dimensions of all
|	CondLikelihood objects stored. If the `cl_stack' is not currently empty and if the new conditional likelihood array
|	length is greater than the current length, all existing objects in `cl_stack' are deleted so that CLAs supplied to
|	the tree in the future will be at least the minimum length needed. Because it is critical that CLAs already checked
|	out to a tree be removed if the new length is longer than the old length (otherwise, CLAs that are too short will
|	be used in the future), this function also checks to make sure there are no CLAs currently checked out. If there are
|	checked out CLAs, an XLikelihood exception is thrown. This function has no effect unless the supplied arguments
|	imply that newly-created CLAs will be longer than the existing ones. It is somewhat wastefull to leave in CLAs that
|	are longer than they need to be, but perhaps more wasteful to continually recall and delete all existing CLAs just
|	to ensure that their length is exactly correct.
*/
void CondLikelihoodStorage::setCondLikeDimensions(const uint_vect_t & np, const uint_vect_t & nr, const uint_vect_t & ns)
	{
	unsigned sz = (unsigned)np.size();
	PHYCAS_ASSERT(nr.size() == sz);
	PHYCAS_ASSERT(ns.size() == sz);

	bool no_old = (num_patterns.size() == 0) || (num_rates.size() == 0) || (num_states.size() == 0);
	bool old = !no_old;

	unsigned newlen = 0;
	unsigned oldlen = 0;
	for (unsigned i = 0; i < sz; ++i)
		{
		PHYCAS_ASSERT(np[i] > 0);
		PHYCAS_ASSERT(nr[i] > 0);
		PHYCAS_ASSERT(ns[i] > 0);
		newlen += np[i]*nr[i]*ns[i];
		if (old)
			oldlen += num_patterns[i]*num_rates[i]*num_states[i];
		}

	if (newlen > oldlen)
		{
		clearStack();
		num_patterns.resize(sz);
		num_rates.resize(sz);
		num_states.resize(sz);
		std::copy(np.begin(), np.end(), num_patterns.begin());
		std::copy(nr.begin(), nr.end(), num_rates.begin());
		std::copy(ns.begin(), ns.end(), num_states.begin());
		}
	}

} // namespace phycas
