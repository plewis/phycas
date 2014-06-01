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
#include <fstream>
#include "ncl/nxsexception.h"
#include "nxs_file_path.hpp"

#if defined(NXS_USE_MAC_DIR_STYLE)
	const char nativeSeparator = ':';
#elif defined(NXS_USE_WINDOWS_DIR_STYLE)
	const char nativeSeparator = '\\';
#elif defined(NXS_USE_UNIX_DIR_STYLE)
	const char nativeSeparator = '/';
#else
#	error must define a filesystem
#endif

void NxsFilePath::ShortenDirStack(std::list<std::string> &dirStack)
	{
	std::list<std::string>::iterator currIt = dirStack.begin();
#	if defined(SUPPORT_ABSOLUTE_PATHS) && !defined(NXS_USE_WINDOWS_DIR_STYLE) //POL, can you do the windows side of ABSOLUTE_PATHS
		if (currIt != dirStack.end() && currIt->length() == 0)
			++currIt;
#	endif
	while (currIt != dirStack.end())
		{
		if (currIt->length() < 1 || (currIt->length() == 1 && (*currIt)[0] == '.')) // remove empty elements and this (.) directories
			currIt = dirStack.erase(currIt);
		else
			{
			if (IsUpDirNotation(*currIt))
				{
				if (currIt != dirStack.begin())
					{
					// shorten  "d/c/../"  to "d/"
					--currIt;
					if (IsUpDirNotation(*currIt))
						++currIt; // we can't erase ../../
					else
						{
						currIt = dirStack.erase(currIt);	//erase the previous directory name
						currIt = dirStack.erase(currIt);	// erase the .. directory name
						}
					}
				}
			++currIt;
			}
		}
	}

void NxsFilePath::TranslateDirStack(const std::list<std::string> &dirStack) const
	{
	std::list<std::string>::const_iterator dIt = dirStack.begin();
	const std::list<std::string>::const_iterator endIt = dirStack.end();

#	if defined(NXS_USE_UNIX_DIR_STYLE)
 		nativePath.clear();
 		if (dirStack.empty())
 			return;
 		if (dirStack.size() == 1 )
 			{
 			nativePath = (IsUpDirNotation(*dIt) ? "../" : *dIt);
 			return;
 			}
 		bool firstWord = true;

 		if (dIt->empty())  // first character was /
 			{
 			// absolute path
 			firstWord = false;
 			}
 		bool lastWasUp = false;
 		for (; dIt != endIt; ++dIt)
 			{
 			if (!firstWord)
 				nativePath.append("/");
 			firstWord = false;
 			lastWasUp = IsUpDirNotation(*dIt);
 			if (lastWasUp)
 				nativePath.append("..");
 			else
 				nativePath.append(*dIt);
 			}
 		if (lastWasUp)
 			nativePath.append("/");
#	elif defined(NXS_USE_MAC_DIR_STYLE)
 		nativePath.clear();
 		bool wordWritten = false;
 		if (dirStack.size() == 1 && !IsUpDirNotation(*dIt))
 			{
 			nativePath = *dIt;
 			return;
 			}
 		for (; dIt != endIt; ++dIt)
 			{
			nativePath.append(":");
 			if (!IsUpDirNotation(*dIt))
 				{
 				wordWritten = true;
 				nativePath.append(*dIt);
 				}
 			}
 		if (!wordWritten)
 			nativePath.append(":");
#	elif defined(NXS_USE_WINDOWS_DIR_STYLE)
		//@POL Mark, submitting a path to _stat to see if it is a directory will fail if the path
		// ends in a backslash ('\\') character. What I'm not sure of is whether elsewhere you depend
		// on having a slash on the end (e.g. when combining with a file name)
		//
 		nativePath.clear();
 		for (; dIt != endIt; ++dIt)
 			{
 			if (IsUpDirNotation(*dIt))
 				nativePath.append("..");
 			else
 				{
                if (dIt != dirStack.begin())
                    {
                    // Only add slash if in the middle of the path
 				    nativePath.append("\\");
                    }
 				nativePath.append(*dIt);
 				}
 			}
#	endif
	}

void NxsFilePath::CreateNative() const
	{
	nativePath = path;
	if (nativePath.empty())
		return;
	if (!nativePath.empty() && nativePath[0] != nativeSeparator)
		{
		std::list<std::string> dirStack;
		split<char, std::string>(nativePath, '/', &dirStack);
		ShortenDirStack(dirStack);
		TranslateDirStack(dirStack);
		}
	}

bool NxsFilePath::OpenInput(std::ifstream &inFStream) const
	{
	inFStream.close();
	inFStream.clear();
	inFStream.open(GetNative().c_str());
	if (inFStream.good())
		return true;
	inFStream.close();
	inFStream.clear();
	return false;
	}

bool NxsFilePath::IsDirectory() const
	{
	if (!this->Exists())
		return false;
	const char *const natDirName = this->GetNative().c_str();
	struct NCL_STAT_STRUCT dStat;
	if (NCL_STAT(natDirName, &dStat) != 0)
		{
		std::string msg;
		msg += "stat call failed for ";
		msg += natDirName;
		throw NxsException(msg);
		}
	bool retCode = NCL_S_ISDIR(dStat.st_mode);
	//std::cerr << "IsDirectory is returning "<< retCode << " for " << natDirName <<std::endl;
	return retCode;
	}

NxsFilePath::NxsFilePath(const std::string &s, bool shouldBeDir )
	:sourceOfName(kProgrammer),
	isAbsolute(false),
	isDirty(true),
	path(),
	pathAsEntered(),
	nativePath(),
	queryUserOnError(true)
	{
	PHYCAS_ASSERT(!shouldBeDir);	// file paths without filenames aren't supported yet
	if (ReadNewPathString(s))
		CreateNative();
	}
bool NxsFilePath::Exists() const
	{
	std::ifstream testIn;
	bool retVal = OpenInput(testIn);
	testIn.close();
	return retVal;
	}

bool NxsFilePath::ReadNewPathString(const std::string &s)
	{
	pathAsEntered = path = s;
	nativePath.clear();
	isDirty = true;
	return true;
	}

bool NxsInFilePath::Open(std::ifstream &outF, std::string purposeOfFile)
	{
	return OpenInput(outF);
	}

