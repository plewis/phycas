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
#include <string>
#include <boost/format.hpp>
#include "cond_likelihood.hpp"

namespace phycas
{

unsigned CondLikelihood::calcCLALength(
  const uint_vect_t & npatterns, 	/**< is a vector containing the number of data patterns for each partition subset */
  const uint_vect_t & nrates, 		/**< is a vector containing the number of among-site relative rate categories for each partition subset */
  const uint_vect_t & nstates) 		/**< is a vector containing the number of states for each partition subset */
	{
	unsigned sz = (unsigned)npatterns.size();
	PHYCAS_ASSERT(nrates.size() == sz);
	PHYCAS_ASSERT(nstates.size() == sz);

	unsigned total = 0;
	for (unsigned i = 0; i < sz; ++i)
		{
		total += (npatterns[i]*nrates[i]*nstates[i]);
		}
	return total;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	CondLikelihood constructor. Allocates `npatterns'*`nrates'*`nstates' elements to the conditional likelihood vector
|	`claVec', allocates `npatterns' elements to the vector `underflowExpon', and sets `numEdgesSinceUnderflowProtection'
|	to UINT_MAX. Sets data member `cla' to point to the first element of the `claVec' vector and `uf' to point to the
|	first element of `underflowExponVec'. All three arguments shoudl be non-zero, and no error checking is done to
|	ensure this because CondLikelihood objects are managed exclusively by CondLikelihoodStorage class, which ensures
|	that the dimensions are valid.
*/
CondLikelihood::CondLikelihood(
  const uint_vect_t & npatterns,	/**< is a vector containing the number of data patterns for each partition subset */
  const uint_vect_t & nrates,		/**< is a vector containing the number of among-site relative rate categories for each partition subset */
  const uint_vect_t & nstates)		/**< is a vector containing the number of states for each partition subset */
  : numEdgesSinceUnderflowProtection(UINT_MAX), total_num_patterns(0)
	{
	unsigned cla_length = calcCLALength(npatterns, nrates, nstates);
	total_num_patterns = (unsigned)std::accumulate(npatterns.begin(), npatterns.end(), 0);
	PHYCAS_ASSERT(total_num_patterns > 0);
	underflowExponVec.resize(total_num_patterns);
	claVec.resize(cla_length);
	cla = &claVec[0];
	uf = &underflowExponVec[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the data member `numEdgesSinceUnderflowProtection'.
*/
unsigned	CondLikelihood::getUnderflowNumEdges() const
	{
	return numEdgesSinceUnderflowProtection;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the data member `numEdgesSinceUnderflowProtection' to `n'.
*/
void	CondLikelihood::setUnderflowNumEdges(
  unsigned n)	/**< is the number of edges traversed since underflow correction was last applied */
	{
	numEdgesSinceUnderflowProtection = n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the conditional likelihood array stored by the vector data member `claVec'.
*/
LikeFltType * CondLikelihood::getCLA()
	{
	return cla;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the conditional likelihood array stored by the vector data member `claVec'.
*/
LikeFltType * CondLikelihood::getCLA() const
	{
	return cla;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the underflow correction array stored by the vector data member `underflowExponVec'.
*/
UnderflowType * CondLikelihood::getUF()
	{
	return uf;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to the conditional likelihood array stored by the vector data member `underflowExponVec'.
*/
UnderflowType const * CondLikelihood::getUF() const
	{
	return uf;	//PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns length of the vector data member `claVec'.
*/
unsigned CondLikelihood::getCLASize() const
	{
	return (unsigned)claVec.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns length of the vector data member `underflowExponVec'.
*/
unsigned CondLikelihood::getUFSize() const
	{
	return (unsigned)underflowExponVec.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all elements in `underflowExponVec' to zero.
*/
void CondLikelihood::zeroUF()
	{
	underflowExponVec.assign(underflowExponVec.size(), 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing a space-separated list of all elements in `underflowExponVec'.
*/
std::string CondLikelihood::debugShowUF() const
	{
	std::string s;
	for (std::vector<UnderflowType>::const_iterator it = underflowExponVec.begin(); it != underflowExponVec.end(); ++it)
		s += str(boost::format("%d ") % (*it));
	return s;
	}

} //namespace phycas

