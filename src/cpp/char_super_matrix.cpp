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

#include "char_super_matrix.hpp"
#include "phycas_nexus_reader.hpp"

// using namespace phycas;
phycas::CharSuperMatrix * GetLastDiscreteMatrix(PhycasNexusReader & nexusReader, bool convertGapsToMissing = false)
	{
	return createNativeDiscreteMatrix(nexusReader, 0L, UINT_MAX, convertGapsToMissing);
	}

phycas::CharSuperMatrix * createNativeDiscreteMatrix(PhycasNexusReader & nexusReader, NxsTaxaBlock * taxaBlockPtr, unsigned int charBlockIndex, bool convertGapsToMissing)
	{
	NxsCharactersBlock * cb;
	if (charBlockIndex == UINT_MAX || taxaBlockPtr == NULL)
		cb = nexusReader.getActiveCharactersBlock();
	else
		cb = nexusReader.GetCharactersBlock(taxaBlockPtr, charBlockIndex);
	if (cb == 0L)
		return 0L;
	std::vector<NxsCXXDiscreteMatrix *> vmats;
	std::vector<const NxsDiscreteDatatypeMapper *> mappers = cb->GetAllDatatypeMappers();
	if (mappers.empty() || mappers[0] == NULL)
		{
		std::string m("Characters block does not contain a matrix with a valid mapping data structure. Please report this error to the developers of Phycas.");
		nexusReader.NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
		return 0L;
		}
	if (mappers.size() > 1)
		{
		std::vector<const NxsDiscreteDatatypeMapper *> o;
		std::map<const NxsDiscreteDatatypeMapper *, NxsUnsignedSet> map2Inds;

		for (unsigned charIndex = 0; charIndex < cb->GetNChar(); ++charIndex)
			{
			const NxsDiscreteDatatypeMapper * currMapper = cb->GetDatatypeMapperForChar(charIndex);
			if (map2Inds.find(currMapper) == map2Inds.end())
				o.push_back(currMapper);
			map2Inds[currMapper].insert(charIndex);
			}

		for (std::vector<const NxsDiscreteDatatypeMapper *>::const_iterator oIt = o.begin(); oIt != o.end(); ++oIt)
			{
			std::map<const NxsDiscreteDatatypeMapper *, NxsUnsignedSet>::const_iterator cm = map2Inds.find(*oIt);
			PHYCAS_ASSERT(cm != map2Inds.end());
			const NxsUnsignedSet & unsSet = cm->second;
			vmats.push_back(new NxsCXXDiscreteMatrix(*cb, convertGapsToMissing, &unsSet));
			}


		}
	else
		vmats.push_back(new NxsCXXDiscreteMatrix(*cb, convertGapsToMissing));
	return new phycas::CharSuperMatrix(vmats);
	}
