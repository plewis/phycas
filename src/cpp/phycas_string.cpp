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

#include <cstdarg>
#include <boost/scoped_array.hpp>
#include "phycas_string.hpp"
#include "ncl/nxsexception.h"
using std::string;
using std::vector;
using std::va_list;
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::strncpy;
		using std::strcpy;
		using std::strchr;
		using std::isgraph;
		using std::isalpha;
		using std::isupper;
		using std::toupper;
		using std::tolower;
		using std::size_t;
#	endif

std::string	Join(const std::vector<std::string> &v, const std::string & separator)
	{
	std::vector<std::string>::const_iterator vIt = v.begin();
	const std::vector<std::string>::const_iterator endIt = v.end();
	std::string retString;
	if (vIt != endIt)
		{
		retString = *vIt++;
		const bool hasSeparator = !separator.empty();
		for (; vIt != endIt; vIt++)
			{
			if (hasSeparator)
				retString << separator;
			retString << *vIt;
			}
		}
	return retString;
	}

bool IsALegalNexusLabelForObjectN(const std::string &s, unsigned n, bool *isANumber)
	{
	unsigned u = 0;
	*isANumber = IsAnUnsigned(s, &u);
	if (*isANumber)
		return u == (n+1);
	return IsLegalNexusWord(s);
	}

void DblFormatter::RecreateFormatCode()
	{
	//@ perhaps we should be using %.g here?
	formatCode.clear();
	if (fieldWidth == UINT_MAX)
		{
		if (digitsAfterDecimal == UINT_MAX)
			formatCode = "%lf";
		else
			StrPrintF(formatCode, "%%.%dlf", digitsAfterDecimal);
		}
	else
		{
		if (digitsAfterDecimal == UINT_MAX)
			StrPrintF(formatCode, "%%%dlf", fieldWidth);
		else
			StrPrintF(formatCode, "%%%d.%dlf", fieldWidth, digitsAfterDecimal);
		}
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Written to make it easy to initialize a vector of strings. Similar to the perl split function. Converts a string like
| this -- "A|bro|ken strin|g" -- to a vector of strings with four elements:  "A", "bro", "ken string", and "g".
*/
vector<string> SplitString(
  const string &str, char separator)     /* the string submitted for splitting */
	{
	vector<string> retVec;
	if (str.empty())
		return retVec;
	string ss;
	string::const_iterator p = str.begin();
	const string::const_iterator endIt = str.end();
	for (;;)
		{
		const bool done = (p == endIt);
		if (done || (*p == separator))
			{
			retVec.push_back(ss);
			ss.clear();
			if (done)
				return retVec;
			++p;
			if (p == str.end())
				return retVec;
			}
		ss << *p++;
		}
	}


size_t CopyToCStr(const string & s, char *buffer, const size_t bufferLen)
	{
	PHYCAS_ASSERT(buffer != NULL);
	PHYCAS_ASSERT(bufferLen > 0);
	size_t len = s.length();
	if (len == 0)
		{
		buffer[0] = '\0';
		return 0;
		}
	if (len >= bufferLen)
		{
		strncpy(buffer, s.c_str(), bufferLen-1);
		buffer[bufferLen-1] = '\0';
		return  (bufferLen-1);
		}
	else
		{
		strcpy(buffer, s.c_str());
		return len;
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns a single-quoted version of the string.
*/
string GetSingleQuoted(const string & f)
	{
	string withQuotes;
	withQuotes.reserve(f.length() + 4);
	withQuotes = f;
	return ConvertToNexusSingleQuoted(withQuotes);
	}

unsigned replace_all_substr(std::string &target, const std::string & searchStr, const std::string & replaceStr)
	{
	unsigned nReplacements = 0;
	size_t i = (unsigned)target.find(searchStr, 0);
	while (i != string::npos)
		{
		++nReplacements;
		target.replace(i, searchStr.length(), replaceStr);
		i = (unsigned)target.find(searchStr, i + (unsigned)replaceStr.length());
		}
	return i;
	}

string & ConvertToNexusSingleQuoted(string & f)
	{
	const string orig(f);
	f = '\'';
	for (std::string::const_iterator oIt = orig.begin(); oIt != orig.end(); ++oIt)
		{
		f << *oIt;
		if (*oIt == '\'')
			f << '\'';
		}
	f << '\'';
	return f;
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if the string is a abbreviation (or complete copy) of the argument `s'.
*/
bool IsStdAbbreviation(
  const std::string &f, ///test string
  const string & s,   /* the string for which the stored string is potentially an abbreviation */
  bool 	respectCase)	     /* if true, comparison will be case-sensitive */
	{
	if (f.empty())
		return false;
		// s is the unabbreviated comparison string
	const unsigned slen = static_cast<unsigned long>(s.size());
	const unsigned flen = static_cast<unsigned long>(f.size());
	if (flen > slen)
		return false;

		// Examine each character in t and return false (meaning "not an abbreviation")
		// if at any point the corresponding character in s is different
	for (unsigned k = 0; k < flen; ++k)
		{
		if (respectCase)
			{
			if (f[k] != s[k])
				return false;
			}
		else if (toupper(f[k]) != toupper(s[k]))
			return false;
		}
	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if `f' is a case-insensitive abbreviation (or complete copy) of `s' and the stored string
| has all of the characters that are in the initial capitalized portion of `s'. For example if `s' is "KAPpa" then
| "kappa", "kapp", or "kap" (with any capitalization pattern) will return true and all other strings will return false.
| Always returns false if the stored string has length of zero.
*/
bool IsCapAbbreviation(
  const string & f,  /* the string to check */
  const string & s)   /* the abbreviation template to compare to */
	{
	if (f.empty())
		return false;

		// s is the unabbreviated comparison string
	const unsigned slen = static_cast<unsigned>(s.size());
	const unsigned flen = static_cast<unsigned>(f.size());
		// If the stored string is longer than s then it cannot be an abbreviation of s
	if (flen > slen)
		return false;

	unsigned k = 0;
	for (; k < slen; ++k)
		{
		if (isupper(s[k]) || !isalpha(s[k]))
			{	// Still in the mandatory portion of the abbreviation
				// Note that non-alphabetic characters in the
				// that are preceded by no lower-case letters are treated as mandatory
			if (k >= flen || toupper(f[k]) != s[k])
				return false;
			}
		else // Get here if we are no longer in the upper case portion of s
			break;
		}

	// Check the lower case portion of s and any corresponding characters in t for mismatches
	// Even though the abbreviation is valid up to this point, it will become invalid if
	// any mismatches are found beyond the upper case portion of s
	//
	for (; k < flen; ++k)
		{
		if (toupper(f[k]) != toupper(s[k]))
			return false;
		}

	return true;
	}
/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if the string needs to be surrounded by single-quotes to make it a single nexus token.
*/
NexusQuotingReq QuotesNeeded(const std::string & s)
	{
	NexusQuotingReq nrq = kNoQuotesNeededForNexus;
	for (std::string::const_iterator sIt = s.begin(); sIt != s.end(); ++sIt)
		{
		if (!isgraph(*sIt))
			{
			if (*sIt != ' ')
				return kSingleQuotesNeededForNexus;
			nrq  = kUnderscoresSufficeForNexus;
			}
		else if (strchr("(){}\"-]/\\,;:=*`+<>", *sIt) != NULL)
			{
			// Get here if c is any NEXUS punctuation mark except left square bracket ([) or apostrophe (').
			// [ and ' never get returned as punctuation by NxsToken,
			// so we should never encounter them here.
			//
			return (s.length() > 1 ? kSingleQuotesNeededForNexus : kNoQuotesNeededForNexus);
			}
		else if (strchr("\'[_", *sIt) != NULL)
			{
			// Get here if c is either an apostrophe or left square bracket. Quotes are needed if one of these
			// characters is all there is to this string
			//
			return kNoQuotesNeededForNexus;
			}
		}
	return nrq;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
bool IsAnUnsigned(const std::string &s, unsigned * u)
	{
	long l;
	if (IsALong(s, &l))
		{
		if (l >=0 && l < (long)UINT_MAX)
			{
			*u = (unsigned) l;
			return true;
			}
		}
	return false;
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if the stored string is a non-case-sensitive copy of the argument `s'. Note: will return true if both the
| stored string and `s' are empty strings.
*/
bool EqualsCaseInsensitive(
  const string &f,
  const string &s)   /* the comparison string */
	{
	unsigned k;
	unsigned slen = (unsigned)s.size();
	unsigned tlen = (unsigned)f.size();
	if (slen != tlen)
		return false;
	for (k = 0; k < tlen; ++k)
		{
		if ((char)toupper(f[k]) != (char)toupper(s[k]))
			return false;
		}
	return true;
	}



/*--------------------------------------------------------------------------------------------------------------------------
| Checks to see if the stored string begins with upper case letters and, if so, returns all of the contiguous capitalized
| prefix. If the stored string begins with lower case letters, an empty string is returned.
*/
string UpperCasePrefix(const string &s)
	{
	string x;
	unsigned i = 0;
	while (i < s.length() && isupper(s[i]))
		x << s[i++];
	return x;
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Transforms the vector of string objects by making them all lower case and then capitalizing the first portion of
| them so that the capitalized portion is enough to uniquely specify each. Returns true if the strings are long enough
| to uniquely specify each. Horrendously bad algorithm, but shouldn't be called often.
*/
bool SetToShortestAbbreviation(
  std::vector<std::string>  & strVec,		/* vector of string objects */
  bool			allowTooShort)  /* */
	{
	std::vector<std::string> upperCasePortion;
	unsigned long i;
	for (i = 0; i < strVec.size(); ++i)
		{
		// Change the next string to lower case
		//
		ToLower(strVec[i]);

		unsigned prefLen = 0;
		string pref;

		if (prefLen >= strVec[i].size())
			return false;
		pref << (char) toupper(strVec[i][prefLen++]);
		bool moreChars = true;

		// Keep adding letters from the current string until pref is unique.
		// Then add this pref to upperCasePortion (vector of previous prefs)
		//
		for (;moreChars;)
			{
			unsigned long prevInd = 0;
			for (; prevInd < upperCasePortion.size(); ++prevInd)
				{
				if (pref == upperCasePortion[prevInd])
					{
					//      Conflict  - both abbreviations need to grow
					//
					if (prefLen >= strVec[i].size())
						{
						if (allowTooShort)
							{
							if (prefLen < strVec[prevInd].size())
								upperCasePortion[prevInd] << (char) toupper(strVec[prevInd][prefLen]);
							moreChars = false;
							break;
							}
						else
							return false;
						}
					pref << (char) toupper(strVec[i][prefLen]);
					if (prefLen >= strVec[prevInd].size())
						{
						if (allowTooShort)
							{
							prevInd = 0;
							++prefLen;
							break;
							}
						else
							return false;
						}
					upperCasePortion[prevInd] << (char) toupper(strVec[prevInd][prefLen++]);
					prevInd = 0;
					break;
					}
				else
					{
					unsigned j;
					for (j = 0; j < prefLen; ++j)
						{
						if (pref[j] != upperCasePortion[prevInd][j])
							break;
						}
					if (j == prefLen)
						{
						//      pref agrees with the first part of another abbreviation, lengthen it.
						//
						if (prefLen >= strVec[i].size())
							{
							if (allowTooShort)
								{
								moreChars = false;
								break;
								}
							else
								return false;
							}
						pref << (char) toupper(strVec[i][prefLen++]);
						break;
						}
					}
				}
			if (prevInd == upperCasePortion.size() || !moreChars)
				{
				// Made it all the way through with no problems, add this
				// prefix as command i's upper case portion
				//
				upperCasePortion.push_back(pref);
				break;
				}
			}
		}

	for (i = 0; i < strVec.size(); ++i)
		{
		for (unsigned j = 0; j < upperCasePortion[i].size(); ++j)
			strVec[i][j] = upperCasePortion[i][j];
		}

	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Returns a vector of string objects that match the entire `testStr'.
*/
std::vector<std::string> GetVecOfPossibleAbbrevMatches(
  const string  & testStr,	       /* string to match */
  const std::vector<std::string> & possMatches)   /* vector of possible matches */
	{
	std::vector<std::string> matches;
	for (unsigned i = 0; i < possMatches.size(); ++i)
		{
		if (Abbreviates(testStr, possMatches[i], ncl::kStringNoRespectCase))
			matches.push_back(possMatches[i]);
		}
	return matches;
	}

void FillVectorWithNumbers(std::vector<std::string> &v, unsigned fromN, unsigned endN)
	{
	v.reserve(v.size() + endN-fromN);
	string s;
	for (; fromN < endN; ++fromN)
		{
		s.clear();
		s << fromN;
		v.push_back(s);
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if the character is a space, tab, newline or < 7.
*/
inline bool IsNexusWhitespace(char ch)
	{
	return (ch == ' ' || ch == '\n' || ch == '\t' || ch < 7);
	}

const char *gIllegalChars = "()[]{}/\\,;:=*\'\"`<>^";

/*--------------------------------------------------------------------------------------------------------------------------
|	Function used when a single character (e.g. Gap character in the characters block) is being assigned by the user
|	and it cannot be whitespace or any one of the following ()[]{}/\\,;:=*\'\"`<>^
|	any of these characters will generate an NxsException.
*/
bool IsLegalNexusChar(
  char testC) 	/* character that the user has chosen */
	{
	return !(IsNexusWhitespace(testC) || strchr(gIllegalChars, testC) != NULL);
	}

void ThrowIllegalNexusChar(
  char 			/*testC*/, 	/* character that the user has chosen */
  const char  * name) /* name of the field that character is being assigned to (e.g. match, missing, gap) used in the NxsException if the character is illegal */
	{
	string s;
	s << name << " character cannot be whitespace (including _) or any of the following:  " << gIllegalChars;
	throw NxsException(s);
	}


	// MakeStrPrintF and StrPrintF are almost identical (MakeStrPrintF would ideally be a wrapper around StrPrintF)
	// but the fact that they do multiple calls to va_start and va_end make it difficult to make one call the other
	//	The solution (for now at least) is to define the body of the function as a huge macro NCL_STR_PRINT_F_MACRO_HACK
	//
	// Comments for the NCL_STR_PRINT_F_MACRO_HACK
	// see  http://www.alphazed.co.uk/programming/vsnprintf.php for a discussion of problems with
	// the return value of vsnprintf.  Basically  we can trust that there was no truncation
	// iff -1 < nAdded < kInitialBufferSize - 1
	//
	//	Comment on the terminating for(;;)
	//  given that the last vsnprintf failed, we either double the buffer size or
	//  increase the buffer size to the previous outsize (with padding for \0 and avoiding
	//  return code ambiguity.

#if defined(_MSC_VER)
#   define VSN_PRINT_F_NAME _vsnprintf
#elif defined (__MWERKS__)
#   define VSN_PRINT_F_NAME  std::vsnprintf
#else
#   define VSN_PRINT_F_NAME vsnprintf
#endif

#define NCL_STR_PRINT_F_MACRO_HACK  	\
	const int kInitialBufferSize = 512; \
	char buf[kInitialBufferSize]; \
	va_list argList; \
	va_start(argList, formatStr); \
	int nAdded = VSN_PRINT_F_NAME(buf, kInitialBufferSize, formatStr, argList); \
	va_end(argList); \
	if (nAdded > -1 && nAdded < kInitialBufferSize - 1) \
		{ \
		str.append(buf); \
		return str; \
		} \
	int nextSizeGuess = kInitialBufferSize; \
	int outsize = nAdded; \
	typedef boost::scoped_array<char> ScopedCharArrPtr; \
	for (;;)  \
		{  \
		nextSizeGuess = (outsize < nextSizeGuess + 1 ? nextSizeGuess * 2 : outsize + 2); \
		ScopedCharArrPtr dynBuffer(new char[nextSizeGuess]); \
		va_start(argList, formatStr); \
		outsize = VSN_PRINT_F_NAME(dynBuffer.get(), (unsigned long) nextSizeGuess, formatStr, argList);\
		va_end(argList); \
		if (outsize > -1 && outsize < nextSizeGuess - 1) \
			{ \
			str.append(dynBuffer.get()); \
			return str; \
			} \
		}

std::string MakeStrPrintF(const char * formatStr, ...)
	{
	std::string str;
	NCL_STR_PRINT_F_MACRO_HACK
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Appends a printf-style formatted string onto the end of this string and returns the number of characters added to the
| string. For example, the following code would result in the string s being set to "ts-tv rate ratio = 4.56789":
|>
| double kappa = 4.56789;
| string s;
| StrPrintF(s, "ts-tv rate ratio = %.5f", kappa);
|>
*/
std::string & StrPrintF(string & str,
  const char * formatStr,	/* the printf-style format string */
  ...)					  /* other arguments referred to by the format string */
	{
	NCL_STR_PRINT_F_MACRO_HACK
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Appends a string version of the supplied unsigned value `v' to the string `s' and returns a reference to `s'.
*/
std::string &append_unsigned(std::string &s, unsigned v)
	{
	static char tmp[128];
	sprintf(tmp, "%d", v);
	s << tmp;
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Utility function useful for reporting strings in exceptions. If supplied string `s' is longer than `max_chars', it
|	is shortened to `max_chars' characters and an ellipsis ("...") is appended to the returned value. If `s' is less
|	than `max_chars', it does not need to be abbreviated and is thus returned in its entirety.
*/
std::string abbreviate(const std::string &s, unsigned max_chars)
	{
	if (s.length() < max_chars)
		return s;
	else
		{
		std::string ss(s.substr(0, max_chars));
		ss << "...";
		return ss;
		}
	}

