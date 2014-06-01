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

#if ! defined(COND_LIKELIHOOD_STORAGE_HPP)
#define COND_LIKELIHOOD_STORAGE_HPP

#include "boost/shared_ptr.hpp"
#include <stdexcept>
#include <vector>
#include <stack>

namespace phycas
{

class CondLikelihood;
typedef boost::shared_ptr<CondLikelihood> CondLikelihoodShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	A stack that stores CondLikelihood shared pointers. CondLikelihoodStorage object can be asked for a pointer to a
|	CondLikelihoodShPtr object when one needed for a likelihood calculation. If the stack is not empty, the pointer on
|	top of the stack is popped off and returned. If the stack is empty, several new CondLikelihood objects are created
|	(`realloc_min' to be exact) on the heap and a pointer to one of them is returned. The `realloc_min' data member is
|	by default 1 but can be modified to improved efficiency if such "stack faults" are expected to be common.
*/
class CondLikelihoodStorage
	{
	public:
										CondLikelihoodStorage();
										//CondLikelihoodStorage(unsigned cond_like_len, unsigned starting_size = 1);
										~CondLikelihoodStorage();

		unsigned                        getCLAPoolSize() const;

		CondLikelihoodShPtr				getCondLikelihood();
		void							putCondLikelihood(CondLikelihoodShPtr p);
		void							fillTo(unsigned capacity);

		unsigned						getNumPatterns(unsigned i) const;
		unsigned						getNumRates(unsigned i) const;
		unsigned						getNumStates(unsigned i) const;
		void							setCondLikeDimensions(const uint_vect_t & np, const uint_vect_t & nr, const uint_vect_t & ns);

		void							setReallocMin(unsigned sz);
		void							clearStack();

		unsigned						bytesPerCLA() const;
		unsigned						numCLAsCreated() const;
		unsigned						numCLAsStored() const;

	private:

		uint_vect_t						num_patterns;	/**< The number of data patterns vector (needed for the CondLikelihood constructor) */
		uint_vect_t						num_rates;		/**< The number of discrete rate categories vector (needed for the CondLikelihood constructor) */
		uint_vect_t						num_states;		/**< The number of states vector (needed for the CondLikelihood constructor) */
		unsigned						num_created;	/**< The total number of CondLikelihood objects created in the lifetime of this object */
		unsigned						realloc_min;	/**< When a request is made and `cl_stack' is empty, `realloc_min' new objects are created and added to the stack */
		std::stack<CondLikelihoodShPtr>	cl_stack;		/**< The stack of CondLikelihoodShPtr */

		unsigned						which;		//TEMP
	};

typedef boost::shared_ptr<CondLikelihoodStorage> CondLikelihoodStorageShPtr;

} //namespace phycas

#endif
