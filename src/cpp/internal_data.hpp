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

#if ! defined(INTERNAL_DATA_HPP)
#define INTERNAL_DATA_HPP

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "ncl/nxsallocatematrix.h"
#include "states_patterns.hpp"
#include "partition_model.hpp"

struct CIPRES_Matrix;

namespace CipresNative
{
class DiscreteMatrix;
}

namespace phycas
{
class CondLikelihood;
class CondLikelihoodStorage;

/*----------------------------------------------------------------------------------------------------------------------
|	The InternalData class stores the data structures needed for computing likelihoods on trees. This includes both
|	conditional likelihood arrays and transition probability matrices. A Conditional Likelihood Array (CLA) is stored
|	internally as a std::vector<double>, but the accessor function getCLA() returns a pointer to the first element so
|	the array can be efficiently traversed during update loops. The CLAs are laid out as follows (for 3 relative rate
|	categores and DNA data):
|>
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	|					 pattern 1					|					 pattern 2					|
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	|	  rate 1	|	  rate 2	|	  rate 3	|	  rate 1	|	  rate 2	|	  rate 3	|	  rate 1	|
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	| A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|					  ^
|>
|	The element above indicated by the symbol ^ would represent the probability of that part of data pattern 1 above
|	this node given that the relative rate was rate 2 and this node had state A. The `ownedPMatrices' data member stores
|	an array of two-dimensional transition matrices similar to the following (for 3 relative rate categories, DNA data).
|	All elements of all of these transition matrices are laid out to occupy one contiguous block of memory; the index of
|	each element in memory is indicated by the numbers inside the cells below. Note that this order is compatible with
|	the layout of the CLA to minimize cache faults as the CLAs are computed.
|>
|				to state			   to state				  to state
|			A	 C	  G	   T	   A	C	 G	  T		  A	   C	G	 T
|	f	 +----+----+----+----+	+----+----+----+----+  +----+----+----+----+
|	r  A |	0 |	 1 |  2 |  3 |	| 16 | 17 | 18 | 19 |  | 32 | 33 | 34 | 35 |
|	o	 +----+----+----+----+	+----+----+----+----+  +----+----+----+----+
|	m  C |	4 |	 5 |  6 |  7 |	| 20 | 21 | 22 | 23 |  | 36 | 37 | 38 | 39 |
|		 +----+----+----+----+	+----+----+----+----+  +----+----+----+----+
|	s  G |	8 |	 9 | 10 | 11 |	| 24 | 25 | 26 | 27 |  | 40 | 41 | 42 | 43 |
|	t	 +----+----+----+----+	+----+----+----+----+  +----+----+----+----+
|	a  T | 12 | 13 | 14 | 15 |	| 28 | 29 | 30 | 31 |  | 44 | 45 | 46 | 47 |
|	t	 +----+----+----+----+	+----+----+----+----+  +----+----+----+----+
|	e			 rate 1					rate 2				   rate 3
|>
*/
class InternalData
  : boost::noncopyable
	{
	friend class TreeLikelihood;
	friend class InternalDataTest;

	public:

													~InternalData();

		unsigned									getCLASize() const;
		double * * *								getPMatrices(unsigned i);
		double * * *								getMutablePMatrices(unsigned i) const;
		const double * const * const *				getConstPMatrices(unsigned i) const;

		CondLikelihoodShPtr							getChildCondLikePtr();
		CondLikelihoodShPtr							getParentalCondLikePtr();
		ConstCondLikelihoodShPtr					getValidChildCondLikePtr() const;
		ConstCondLikelihoodShPtr					getValidParentalCondLikePtr() const;

		bool										filialCLAValid() const;
		bool										filialCLACached() const;
		bool										parentalCLAValid() const;
		bool										parentalCLACached() const;

		CondLikelihoodStorageShPtr					getCondLikeStoragePool() const {return cla_pool;}
	private:
													InternalData(PartitionModelShPtr partition, CondLikelihoodStorageShPtr cla_storage);

		//CLA's for an edge from a node to its parent are stored in the node's InternalData (or TipData).
		//bool										parCLAValid;	/**< true if parWorkingCLA is valid */
		CondLikelihoodShPtr							parWorkingCLA;	/**< conditional likelihood array for the parent node and beyond (valid if it points to something, invalid otherwise) */
		CondLikelihoodShPtr							parCachedCLA;	/**< parental conditional likelihood array is stored here to make reverting MCMC moves cheap */
		//bool										childCLAValid;	/**< true if childWorkingCLA is valid. */
		CondLikelihoodShPtr							childWorkingCLA;/**< conditional likelihood array for this node and above (valid if it points to something, invalid otherwise) */
		CondLikelihoodShPtr							childCachedCLA; /**< filial conditional likelihood array is stored here to make reverting MCMC moves cheap */

		state_code_t								state;			/**< Used in simulation to temporarily store the state for one character */
		std::vector< ScopedThreeDMatrix<double> >	pMatrices;		/**< pMatrix[s][r] is the transition matrix for subset s and relative rate r */
		CondLikelihoodStorageShPtr					cla_pool;		/**< CondLikelihood object storage facility */
	};

typedef boost::shared_ptr<InternalData>			InternalDataShPtr;
typedef std::vector<InternalDataShPtr>			VecInternalDataShPtr;

// **************************************************************************************
// ***** InternalData inlines ***********************************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrices'.
*/
inline double * * * InternalData::getPMatrices(
  unsigned i)		/**< is the subset of the partition */
	{
	return pMatrices[i].ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrices'.
*/
inline double * * * InternalData::getMutablePMatrices(
  unsigned i) const		/**< is the subset of the partition */
	{
	return pMatrices[i].ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrices'.
*/
inline const double * const * const * InternalData::getConstPMatrices(
  unsigned i) const		/**< is the subset of the partition */
	{
	return pMatrices[i].ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `childWorkingCLA' data member. If `childWorkingCLA' does not currently point to anything, a CondLikelihood
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
inline CondLikelihoodShPtr InternalData::getChildCondLikePtr()
	{
	if (!childWorkingCLA)
		childWorkingCLA = cla_pool->getCondLikelihood();
	return childWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member. If `parWorkingCLA' does not currently point to anything, a CondLikelihood
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
inline CondLikelihoodShPtr InternalData::getParentalCondLikePtr()
	{
	if (!parWorkingCLA)
		parWorkingCLA = cla_pool->getCondLikelihood();
	return parWorkingCLA;
	}

inline ConstCondLikelihoodShPtr InternalData::getValidChildCondLikePtr() const
	{
	//@POL-NESCENT Mark, I don't understand this - why not just assert that childWorkingCLA actually
	// points to a CondLikelihood object, then return childWorkingCLA directly?
	//PHYCAS_ASSERT(childCLAValid);
	//InternalData * t = const_cast<InternalData *>(this);
	//return t->getChildCondLikePtr();
	PHYCAS_ASSERT(childWorkingCLA);
	return ConstCondLikelihoodShPtr(childWorkingCLA);
	}

inline ConstCondLikelihoodShPtr InternalData::getValidParentalCondLikePtr() const
	{
	//PHYCAS_ASSERT(parCLAValid);
	//InternalData * t = const_cast<InternalData *>(this);
	//return t->getParentalCondLikePtr();
	PHYCAS_ASSERT(parWorkingCLA);
	return ConstCondLikelihoodShPtr(parWorkingCLA);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `childWorkingCLA' data member actually points to something, which means that the parental CLA is
|	valid. Returns false if `childWorkingCLA' does not point to a CondLikelihood object.
*/
inline bool InternalData::filialCLAValid() const
	{
	return (bool)childWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `childCachedCLA' data member actually points to something, which means the working CLA has been
|	cached in case a revert of the current move is needed. Returns false if `childCachedCLA' does not point to a
|	CondLikelihood object.
*/
inline bool InternalData::filialCLACached() const
	{
	return (bool)childCachedCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parWorkingCLA' data member actually points to something, which means that the parental CLA is
|	valid. Returns false if `parWorkingCLA' does not point to a CondLikelihood object.
*/
inline bool InternalData::parentalCLAValid() const
	{
	return (bool)parWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parCachedCLA' data member actually points to something, which means the working CLA has been cached
|	in case a revert of the current move is needed. Returns false if `parCachedCLA' does not point to a CondLikelihood
|	object.
*/
inline bool InternalData::parentalCLACached() const
	{
	return (bool)parCachedCLA;
	}

} // namespace phycas

#endif
