/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2010 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
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

#if ! defined(CHAR_SUPER_MATRIX_HPP)
#define CHAR_SUPER_MATRIX_HPP

#include <boost/shared_ptr.hpp>
#include "ncl/nxspublicblocks.h"
#include "nxs_file_path.hpp"
#include "ncl/nxscxxdiscretematrix.h"

class PhycasNexusReader;
namespace phycas
{
class CharSuperMatrix
	{
	public:
		CharSuperMatrix(){}
		CharSuperMatrix(std::vector<NxsCXXDiscreteMatrix *> mats) : allSubmatrices(mats) {}

		unsigned GetNumMatrices() const {return allSubmatrices.size();}
		unsigned getNTax() const {return allSubmatrices[0]->getNTax();}
		unsigned getNChar() const {return allSubmatrices[0]->getNChar();}
		NxsCXXDiscreteMatrix * GetMatrix(unsigned i) const { return allSubmatrices.at(i);}

	public:

		std::vector<NxsCXXDiscreteMatrix *> allSubmatrices;
	};


}	// namespace phycas

phycas::CharSuperMatrix * GetLastDiscreteMatrix(PhycasNexusReader & nexusReader, bool convertGapsToMissing);
phycas::CharSuperMatrix * createNativeDiscreteMatrix(PhycasNexusReader & nexusReader, NxsTaxaBlock * taxaBlockPtr, unsigned int charBlockIndex, bool convertGapsToMissing);

#endif
