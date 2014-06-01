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

#include <iostream>
#include <iterator>
#include <fstream>
#include "ncl/nxsexception.h"
#include "ncl/nxstreesblock.h"
#include "ncl/nxscxxdiscretematrix.h"
#include "phycas_nexus_reader.hpp"
#include "char_super_matrix.hpp"

PhycasNexusReader::~PhycasNexusReader()
	{
	//std::cerr << "\n\n>>>>> PhycasNexusReader dying..." << std::endl;
	// Unlike the PublicNexusReader, the Phycas Nexus Reader owns the blocks, and thus deletes them on destruction
	Clear();
	}

void PhycasNexusReader::Clear()
	{
	PublicNexusReader::DeleteBlocksFromFactories();
	}

// void PhycasNexusReader::NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult , NxsBlock* )
// 	{
// 	errorMsg = msg;
// 	filePositionOfError = pos;
// 	lineOfError = line;
// 	columnOfError = col;
// 	}

void PhycasNexusReader::ReadNxsFilePath(NxsInFilePath & filePath)
	{
	errorMsg.clear();
	const std::string fn = filePath.GetFullName();
	if (!filePath.Exists())
		{
		errorMsg = std::string("The file \'");
		errorMsg.append(fn);
		errorMsg.append("\' does not exist.");
		throw NxsException(errorMsg.c_str());
		}
	if (filePath.IsDirectory())
		{
		errorMsg  = std::string ("\'");
		errorMsg.append(fn);
		errorMsg.append("\' is a directory, not a file.");
		throw NxsException(errorMsg.c_str());
		}
	std::ifstream inStream;
	if (!filePath.Open(inStream))
		{
		errorMsg = std::string("The file \'");
		errorMsg.append(fn);
		errorMsg.append("\' could not be opened.");
		throw NxsException(errorMsg.c_str());
		}
	std::string msg;
	NxsToken token(inStream);
	Execute(token);
	const unsigned nTaxaBlocks = GetNumTaxaBlocks();
	if (nTaxaBlocks == 0)
		{
		std::string m("No taxa blocks encountered in ");
		m.append(fn);
		m.append(" No information from the file was stored.");
		NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
		activeTaxaBlockIndex = UINT_MAX;
		activeCharactersBlockIndex = UINT_MAX;
		return;
		}
	else
		{
		activeTaxaBlockIndex = nTaxaBlocks - 1;
		if (nTaxaBlocks > 1)
			{
			std::string m("More than one taxa block encountered, only the last Taxa block is active");
			NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
			}
		trees.clear();
		}
	NxsTaxaBlock * activeTB = GetTaxaBlock(activeTaxaBlockIndex);
	PHYCAS_ASSERT(activeTB);

	const unsigned nCharBlocks = GetNumCharactersBlocks(activeTB);
	if (nCharBlocks == 0)
		activeCharactersBlockIndex = UINT_MAX;
	else
		{
		activeCharactersBlockIndex = nCharBlocks - 1;
		if (nCharBlocks > 1)
			{
			std::string m("More than one Characters block encountered, only the last Characters block is active");
			NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
			}
		}

	const unsigned nTreesBlocks = GetNumTreesBlocks(activeTB);
	for (unsigned tInd = 0; tInd < nTreesBlocks; ++tInd)
		{
		NxsTreesBlock * tb = GetTreesBlock(activeTB, tInd);
		std::vector<NxsFullTreeDescription> & currBlockTrees = tb->GetProcessedTrees();
		/* PELIGROSO -- we are stealing the trees from the first trees block */
		if (trees.empty())
			trees.swap(currBlockTrees);
		else
			std::copy(currBlockTrees.begin(), currBlockTrees.end(), std::back_inserter(trees));
		tb->Reset();
		}
	}

NxsTaxaBlock * PhycasNexusReader::getActiveTaxaBlock() const
	{
	if (activeTaxaBlockIndex == UINT_MAX)
		return NULL;
	return GetTaxaBlock(activeTaxaBlockIndex);
	}

NxsCharactersBlock * PhycasNexusReader::getActiveCharactersBlock() const
	{
	NxsTaxaBlock * activeTB = getActiveTaxaBlock();
	if (activeTB == 0L)
		return 0L;
	if (activeCharactersBlockIndex == UINT_MAX)
		return 0L;
	return GetCharactersBlock(activeTB, activeCharactersBlockIndex);
	}

unsigned PhycasNexusReader::GetNChar() const
	{
	NxsCharactersBlock * cb = getActiveCharactersBlock();
	if (cb == 0L)
		return 0;
	return cb->GetNumChar();
	}

const std::vector<NxsFullTreeDescription> & PhycasNexusReader::GetTrees() const
	{
	return trees;
	}


std::vector<std::string> PhycasNexusReader::GetTaxLabels() const
	{
	NxsTaxaBlock * activeTB = getActiveTaxaBlock();
	if (activeTB == 0L)
		return std::vector<std::string>();
	return activeTB->GetAllLabels();
	}
