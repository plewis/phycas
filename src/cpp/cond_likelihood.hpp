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

#if ! defined(COND_LIKELIHOOD_HPP)
#define COND_LIKELIHOOD_HPP

#include <vector>
#include <stack>
#include <boost/shared_ptr.hpp>
#include "states_patterns.hpp"

namespace phycas
{
typedef double LikeFltType;
typedef long UnderflowType;

/*----------------------------------------------------------------------------------------------------------------------
|	Manages a conditional likelihood array for one end of an edge.
*/
class CondLikelihood
	{
	public:

									CondLikelihood(const uint_vect_t & npatterns, const uint_vect_t & nrates, const uint_vect_t & nstates);

		LikeFltType *				getCLA();
		LikeFltType *				getCLA() const;
		unsigned					getCLASize() const;

		UnderflowType *				getUF();
		UnderflowType const *		getUF() const;
		unsigned					getUFSize() const;

		void						zeroUF();
		std::string					debugShowUF() const;

		unsigned					getUnderflowNumEdges() const;
		void						setUnderflowNumEdges(unsigned n);

		static unsigned				calcCLALength(const uint_vect_t & npatterns, const uint_vect_t & nrates, const uint_vect_t & nstates);

	private:

		LikeFltType *				cla;								/**< Pointer to conditional likelihood array stored by `claVec'  */
		std::vector<LikeFltType>	claVec;								/**< Each element contains the likelihood conditional on a particular state, rate and pattern */

		UnderflowType *				uf;									/**< Pointer to the underflow correction array stored in `underflowExponVec'. Used if UnderflowManager is in effect */
		std::vector<UnderflowType>	underflowExponVec;					/**< Stores log of the underflow correction factor for each pattern. Used if UnderflowManager is in effect */

		unsigned 					numEdgesSinceUnderflowProtection;	/**< The number of edges traversed since the underflow protection factor was last updated */

		unsigned 					total_num_patterns;					/**< The number of patterns over all partition subsets */
	};

typedef boost::shared_ptr<CondLikelihood> CondLikelihoodShPtr;
typedef boost::shared_ptr<const CondLikelihood> ConstCondLikelihoodShPtr;

} //namespace phycas

#endif
