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

#ifndef PHYCAS_SPLIT_H
#define PHYCAS_SPLIT_H

#include <vector>
#include <set>

namespace phycas
{

typedef unsigned long           split_t;    // Split::InvertSplit assumes this is an unsigned integer type
typedef std::vector<split_t>    SplitTVect;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the notion of a taxon bipartition, or split. Each split is stored as a collection of bits, with the 
|	bits that are set representing the taxa above an edge in the tree. The type split_t defines the number of bits in 
|	one unit. If there are more taxa than bits in one unit, the split must be represented by multiple units. The 
|   function CalcNUnits figures out the number of units that must be used.
*/
class Split
	{
    friend class Tree;
	
	public:

							    Split();
							    Split(const Split & other);
							    ~Split();

        void                    Copy(const Split & other);
		
		std::string 	        GetDimensionInfo();
		unsigned 		        GetNTaxa() const;
		void                    SetNTaxa(unsigned n);
		void 			        CalcNUnits(unsigned ntax);

		unsigned			    CalcComplexity() const;

        std::string             Write() const;
        void                    Read(const std::string s);

        void                    CreateFromPattern(std::string s);
		void	 			    CreateAndAppendPatternRepresentation(std::string &) const;
		std::string			    CreateIdRepresentation() const;
		std::string			    CreatePatternRepresentation() const;
        std::string             CreateNewickRepresentation(bool zero_based = false) const;

        unsigned			    CountOnBits() const;
		unsigned			    CountOffBits() const;
		int 				    Cmp(const Split & other) const;
		bool 				    Equals(const Split & other) const;
		bool 				    IsBitSet(unsigned t) const;
		bool				    IsCompatible(const Split & other) const;
		bool 				    IsLessThan(const Split & other) const;
		bool 				    SubsumedIn(const Split & other, unsigned startUnit = 0) const;

        std::vector<unsigned>   GetOnList() const;
        std::vector<unsigned>   GetOffList() const;
        void                    SetExcluded(const std::vector<unsigned> excl);
        std::vector<unsigned>   GetExcludedList() const;

        void                    SetOnSymbol(char c);
        void                    SetOffSymbol(char c);
        void                    SetExcludedSymbol(char c);

        char                    GetOnSymbol() const;
        char                    GetOffSymbol() const;
        char                    GetExcludedSymbol() const;

		bool				    operator!=(const Split & other) const;
		Split 				    operator&(const Split & other) const;
		Split &                 operator&=(const Split & other);
		Split 				    operator^(const Split & other) const;
		Split &                 operator^=(const Split & other);
		Split 				    operator|(const Split & other) const;
		Split &                 operator|=(const Split & other);
		bool				    operator<(const Split & other) const;
		bool				    operator==(const Split & other) const;
		Split &                 operator=(const Split & other);

		friend std::istream &   operator>>(std::istream & in, Split & s);

        //@POL 13-Nov-2007 it is currently a mystery why uncommenting the line below leads to hundreds of VC error C2679 messages: 
        // binary '<<' : no operator found which takes a right-hand operand of type <whatever> (or there is no acceptable conversion)
        //friend std::string &    operator<<(std::string & out, const Split & s);

		void 				    Clear();
        void                    Reset();

		void 				    CombineWith(const Split & other);
		void 				    IntersectWith(const Split & other);

		void 				    SetBit(unsigned t);
        void 				    SetBits(const std::vector<unsigned> bits_to_set);
		void 				    UnsetBit(unsigned t);
        void 				    UnsetBits(const std::vector<unsigned> bits_to_unset);

		void 				    InvertSplit();
			
	private:

		void 				    Resize();

        void                    GetOnListImpl(std::vector<unsigned> & v) const;
        void                    GetOffListImpl(std::vector<unsigned> & v) const;

        void                    WriteImpl(std::string & out) const;
        void                    ReadImpl(std::istream & in);
		
    public:

		SplitTVect              unit;			/**< is the vector of split units */
		const unsigned		    bits_per_unit;	/**< is the number of bits in a variable of type split_t */
		const split_t		    split_unity;	/**< is a split_t variable with only the least significant bit set */
		unsigned		        split_ntax;		/**< is the number of taxa currently under consideration */
		unsigned		        nunits;			/**< is the length of the array necessary to represent a split */
		char			        on_symbol;		/**< is the symbol used to represent bits that are set (i.e. "on") */
		char			        off_symbol;		/**< is the symbol used to represent bits that have been cleared (i.e. "off") */
		char			        excl_symbol;	/**< is the symbol used to represent bits that have been excluded (i.e. should not be considered off or on) */
        std::vector<unsigned>   excl_bits;      /**< is a sorted vector containing bit positions corresponding to excluded taxa (lowest possible value is 0) */
	};

std::istream & operator>>(std::istream & in, Split & s);

typedef std::set<Split> SplitSet;

unsigned FindSplitsAbsentInTestTree(const SplitSet & ref_tree, const SplitSet & test_tree, SplitSet & missing);

}   //namespace phycas

#endif

