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

#if ! defined(TIP_DATA_HPP)
#define TIP_DATA_HPP

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
class TreeLikelihood;
class CondLikelihood;
class CondLikelihoodStorage;
typedef boost::shared_ptr<CondLikelihood> CondLikelihoodShPtr;
typedef boost::shared_ptr<const CondLikelihood> ConstCondLikelihoodShPtr;
class Tree;

/*----------------------------------------------------------------------------------------------------------------------
|	The TipData class stores the data associated with each tip in the tree. The description below assumes the data was
|	read from the following NEXUS data file (ambig.nex):
|>
|	#nexus
|
|	begin data;
|	  dimensions ntax=4 nchar=5;
|	  format datatype=dna missing=? gap=-;
|	  matrix
|		taxon1 A ? T  G	   (CGT)
|		taxon2 A C - (AGT)	T
|		taxon3 A C -  N	   (AG)
|		taxon4 A C T  R		Y
|	  ;
|	end;
|>
|	Assume that a partition with 2 subsets is specified. The first subset comprises sites 1-3, and the second subset
|	comprises sites 4 and 5. The data in the above NEXUS file would be stored internally in the form of state codes:
|	states or state combinations are replaced with integers ranging from -1 to 9.
|>
|	taxon1	0 4	 3 2 6
|	taxon2	0 1 -1 7 3
|	taxon3	0 1 -1 5 8
|	taxon4	0 1	 3 8 9
|>
|	The global state list (TreeLikelihood::state_list) for both subsets begins by defining the four primary DNA states
|	(A, C, G, T) and the completely ambiguous state (?). Note that ? includes all four primary states as well as
|	the gap state. Each state is preceded by the number of elements in the state list required to define that state.
|	Thus, each primary state is preceded by 1 because only one integer element is required to define a primary state.
|	The ? state requires 5 elements - gap (-1), A (0), C (1), G (2), and T (3) - and is thus preceded by the number 5.
|	The state list for the second subset must define 4 partially-ambiguous states: (CGT), (AGT), R=(AG) and Y=(CT).
|>
|	state_list[0]:
|	   1 0 1 1 1 2 1 3 5 -1 0 1 2 3
|	     A   C   G   T    - A C G T
|
|	state_list[1]:
|	   1 0 1 1 1 2 1 3 5 -1 0 1 2 3 4 0 1 2 3 3 1 2 3 3 0 2 3 2 0 2 2 1 3
|	     A   C   G   T    - A C G T   A C G T   C G T   A G T   A G   C T
|>
|	The global state list position vector (TreeLikelihood::state_list_pos) specifies the offset into the global state
|	lists needed to access each state.
|>
|	  state_list_pos[0]: 0 2 4 6 8
|	  state_list_pos[1]: 0 2 4 6 8 14 19 23 27 30
|>
|	For large datasets with many different ambiguity combinations, the global state lists could grow quite large. Each
|	state code ends up being a row in the augmented transposed transition probability matrix used in computing
|	conditional likelihood arrays, and it could be quite inefficient if these augmented matrices had many unnecessary
|	rows. Thus, when a row in a global state list is copied to the tip of a tree, where it is used to compute likelihoods,
|	a translation is done into local state codes. The local codes are identical to the global codes for the primary states
|	and the state representing complete ambiguity, but may differ from the global codes for codes representing partial
|	ambiguities. To illustrate, consider taxon4 in the example. It has these states in the NEXUS data file:
|	  A	 C	T  R  Y
|	which translate to these global state codes:
|>
|	  TreeLikelihood::state_list[0]: 0 1 3
|	  TreeLikelihood::state_list[1]: 8 9
|>
|	When copied to a tip node in a tree, however, these codes become translated to
|>
|	  TipData::state_list[0]: 0 1 3
|	  TipData::state_list[1]: 5 6
|>
|	Global codes 8 and 9 have become local codes 5 and 6, respectively.
|
|	Note that the TipData constructors are private and thus TipData objects can only be created by the friend function
|	allocateTipData().
*/
class TipData
  : boost::noncopyable
	{
	friend class TreeLikelihood;
	friend class SimData;

	public:
													TipData(PartitionModelShPtr partition, const state_list_pos_vect_t & positions, const state_list_vect_t & states, CondLikelihoodStorageShPtr cla_storage);
													TipData(PartitionModelShPtr partition, CondLikelihoodStorageShPtr cla_storage);// this constructor only used for simulations
                                           	 		~TipData();

		const double * const * const *				getConstTransposedPMatrices(unsigned i) const;
		double * * *								getTransposedPMatrices(unsigned i);
		double * * *								getMutableTransposedPMatrices(unsigned i) const;
		const state_code_t *						getConstStateCodes(unsigned i) const;

		CondLikelihoodShPtr							getParentalCondLikePtr();
		ConstCondLikelihoodShPtr					getValidParentalCondLikePtr() const;

		bool										parentalCLAValid() const;
		bool										parentalCLACached() const;

        state_list_t & 								getTipStatesArray(unsigned i) {return state_codes.at(i);}
        const state_list_t &						getConstTipStatesArray(unsigned i) const {return state_codes.at(i);}
		const uint_vect_t &							getConstStateListPos(unsigned i) const;

		friend void									calcPMatTranspose(const TreeLikelihood & treeLikeInfo, const TipData & tipData, double edgeLength);

	private:

		//bool										parCLAValid;
		CondLikelihoodShPtr							parWorkingCLA;		/**< conditional likelihood array for parent and beyond (valid if it points to something, invalid otherwise) */
		CondLikelihoodShPtr							parCachedCLA;		/**< parental conditional likelihood array is stored here to make reverting MCMC moves cheap */


		state_code_t								state;				/**< Used in simulation to temporarily store the state for one character */
		state_list_pos_vect_t						state_list_pos;		/**< Vector of indices into the tip-specific `state_codes' array */
		state_list_vect_t							state_codes;		/**< Array of tip-specific state codes */
		std::vector< ScopedThreeDMatrix<double> >	pMatrixTranspose;	/**< pMatrixTranspose[s][r] is the transposed transition matrix for subset s and relative rate r */
		CondLikelihoodStorageShPtr					cla_pool;			/**< Source of CondLikelihood objects if needed */
	};

typedef boost::shared_ptr<TipData> TipDataShPtr;
typedef std::vector<TipDataShPtr> VecTipDataShPtr;

// *********************************************************************************
// ***** TipData inlines ***********************************************************
// *********************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member as a shared pointer to const CondLikelihood. If `parWorkingCLA' does not
|	currently point to anything, a CondLikelihood object is first retrieved from `cla_pool', so this function always
|	returns a shared pointer that actually points to something.
*/
inline ConstCondLikelihoodShPtr TipData::getValidParentalCondLikePtr() const
	{
	//PHYCAS_ASSERT(parCLAValid);
	//TipData * t = const_cast<TipData *>(this);
	//return t->getParentalCondLikePtr();
	PHYCAS_ASSERT(parWorkingCLA);
	return ConstCondLikelihoodShPtr(parWorkingCLA);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the element of `pMatrixTranspose' (const) corresponding to subset `i'.
*/
inline const double * const * const * TipData::getConstTransposedPMatrices(
  unsigned i) const		/**< is the subset */
	{
	return pMatrixTranspose[i].ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the element of `pMatrixTranspose' (non-const) corresponding to subset `i'.
*/
inline double * * * TipData::getTransposedPMatrices(
  unsigned i)		/**< is the subset */
	{
	return pMatrixTranspose[i].ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrixTranspose' (non-const).
*/
inline double * * * TipData::getMutableTransposedPMatrices(
  unsigned i)		/**< is the subset */
  const
	{
	return pMatrixTranspose[i].ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the address of the first element of the vector of state codes for subset i.
*/
inline const state_code_t * TipData::getConstStateCodes(
  unsigned i)		/**< is the subset */
  const
	{
	return &(state_codes[i][0]);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `state_list_pos' data member. The returned vector of unsigned values holds
|	the position of each particular coded state for this tip in the local `state_codes' vector for subset i.
*/
inline const uint_vect_t & TipData::getConstStateListPos(
  unsigned i) 		/**< is the subset */
  const
	{
	return state_list_pos[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parWorkingCLA' data member actually points to something, which means that the parental CLA is
|	valid. Returns false if `parWorkingCLA' does not point to a CondLikelihood object.
*/
inline bool TipData::parentalCLAValid() const
	{
	return parWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parCachedCLA' data member actually points to something, which means the working CLA has been cached
|	in case a revert of the current move is needed. Returns false if `parCachedCLA' does not point to a CondLikelihood
|	object.
*/
inline bool TipData::parentalCLACached() const
	{
	return parCachedCLA;
	}

} // namespace phycas

#endif
