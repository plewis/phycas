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

#if !defined(PYPHY_STRING_HPP)
#define PYPHY_STRING_HPP

//POL 29-Nov-2005
//
// This file was originally intended to be a replacement for string_extensions.hpp and a workspace
// for trying new things (e.g. boost::format) without tampering with working code that depends on
// string_extensions.hpp. Unfortunately, there are too many places in the phycas code where I needed
// to include headers from ncl, which eventually leads to a conflict between string_extensions.hpp
// and phycas_string.hpp. So, now this file is largely useless and should be eliminated. I still like
// my DoubleFormatter function, however. It does the same job as string_extensions' DblFormatter,
// but uses boost::format (which is slower than sprintf but has the benefit of being typesafe).
// Thus, I have kept DoubleFormatter but now include string_extensions.hpp to get the other stuff.

#include <string>
#include <limits>
#include "boost/format.hpp"
#include <climits>
#include <string>
#include <vector>
#if defined(__APPLE__) && defined(__GNUC__)
#	include <math.h>
#	include <cfloat>
#	include <ctype.h>
#else
#	include <cmath>
#	include <cctype>
#endif
#include <cstring>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cctype>
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::strlen;
#endif
#include <limits>
#include <ncl/nxsstring.h>
class NxsException;

#if defined(_MSC_VER)
#	if _MSC_VER >= 1500
#		include <cstdio>
#		define sprintf sprintf_s
#		if !defined(vsnprintf)
#			define vsnprintf _vsnprintf_s
#		endif
#	endif
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Class used with the << operators to format doubles by providing the arguments for the sprintf function that is
|	used to append double.  
|	The first argument to the constructor is the minimum width of the double,
|	The second argument is the number of digits of precision to use
*/
class DblFormatter
	{
	public:
		DblFormatter(unsigned widthOfField, unsigned nDigitsAfterDecimal)
			:fieldWidth(widthOfField),
			digitsAfterDecimal(UINT_MAX)
			{
			SetPrecision(nDigitsAfterDecimal); //crops precision if needed and initializes formatcode
			}
		unsigned GetWidth() const
			{
			return fieldWidth;
			}
		unsigned GetPrecision() const
			{
			return digitsAfterDecimal;
			}
		std::string & FormatDouble(std::string & s, const double d) const;		
		void SetWidth(unsigned int widthOfField)
			{
			fieldWidth = widthOfField;
			RecreateFormatCode();
			}
		
		void SetPrecision(unsigned int nDigitsPrecision)
			{
			digitsAfterDecimal =  (unsigned int) (nDigitsPrecision > kMaxDoublePrintPrecision ? kMaxDoublePrintPrecision : nDigitsPrecision);
			RecreateFormatCode();
			}
		
		unsigned fieldWidth;
		unsigned digitsAfterDecimal;
	private:
		std::string formatCode;
		STATIC_CONST const unsigned int kMaxDoublePrintPrecision = std::numeric_limits<double>::digits10; 
		void RecreateFormatCode();
	};
class StringDblFormatterRef
	{
	public:
		StringDblFormatterRef(std::string &s, DblFormatter & d)
			:strRef(s),
			dfRef(d)
			{}
		std::string &strRef;
		const DblFormatter dfRef;
	};
namespace ncl
{
enum StringCmpEnum /* enum that is used to specify the mode for string comparisons/equality tests */
	{
	kStringRespectCase , 
	kStringNoRespectCase, 
	kStringCapAbbrev
	};

enum	IntegerConversion /*argument type to control how forgiving ConvertToUnsigned ConvertToInt and ConvertToLong are */
	{
	kStringIsInt,		/* throws NxsX_NotANumber if the string (in its entirety) isn't an integer */
	kStringIsDblRound,	/* throws NxsX_NotANumber if the string (in its entirety) isn't a number.  Rounds in the case of doubles */
	kStringIsDblFloor	/* throws NxsX_NotANumber if the string (in its entirety) isn't a number.  returns floor in the case of doubles */
	};
} //namespace ncl

	
/*--------------------------------------------------------------------------------------------------------------------------
| Written to make it easy to initialize a vector of strings. Similar to the perl split function. Converts a string like
| this -- "A|bro|ken strin|g" -- to a vector of strings with four elements:  "A", "bro", "ken string", and "g".
*/
std::string	  & AppendNumberThenWord(unsigned i, std::string s);
std::string	  & Capitalize(std::string &);
void			CharToUpper(char & c);
void			CharToLower(char & c);
std::size_t		CopyToCStr(const std::string &, char * buffer, const std::size_t bufferLen) ;
bool			EqualsLong(const std::string &, long) ;
bool			EqualsCaseInsensitive(const std::string &, const std::string &s) ;
std::string		GetCapitalized(const std::string &) ;
std::string		GetOrderString(unsigned ind);
bool			IsADouble(const std::string &, double *) ;
bool 			IsALong(const std::string &, long *) ;
bool	 		IsAnUnsigned(const std::string &, unsigned *) ;
bool			IsInVector(const std::vector<std::string> &s, ncl::StringCmpEnum mode = ncl::kStringRespectCase) ;
bool			IsSpace(char c);
bool			IsStdAbbreviation(const std::string &,  const std::string &s, bool respectCase) ;
std::string		Join(const std::vector<std::string> &toConcatenate, const std::string &separator);
std::string	  	MakeStrPrintF(const char * formatStr, ...);
std::string		MakeRightJustifiedString(const std::string &s, unsigned w);
std::string		MakeRightJustifiedLong(long x, unsigned w);
std::string		MakeRightJustifiedDbl(double x, unsigned w, unsigned p);
unsigned char * p_str(unsigned char *) ;
unsigned		replace_all_substr(std::string &target, const std::string & searchStr, const std::string & replaceStr); // mimics stl, but takes string not iterator
std::string	  & RightJustifyString(std::string &buff, const std::string &s, unsigned w);
std::string	  & RightJustifyLong(std::string &s, long x, unsigned w);
std::string	  & RightJustifyDbl(std::string &s, double x, unsigned w);
std::vector<std::string>		SplitString(const std::string & strList, char separator);
bool			StrEquals(const std::string &, const std::string & s, ncl::StringCmpEnum mode);
std::string	  &	StrPrintF(std::string & str, const char * formatStr, ...);
std::string	  & ToLower(std::string &);
std::string	  & ToUpper(std::string &);
template<typename T>
std::string  & AppendNumber(std::string &, const T);

//typedef std::pair<std::string &, DblFormatter &> StringDblFormatterRef;
//typedef std::pair<std::string &, DblFormatter> StringDblFormatter;
StringDblFormatterRef operator<<(std::string &, DblFormatter & d); // logically the DblFormatter is const, but we need an object (want to disallow construction of a temporary)
//StringDblFormatter operator<<(std::string &, const DblFormatter & d);
std::string operator<<(StringDblFormatterRef &, const double d);
//std::string operator<<(StringDblFormatter &, const double d);


std::string & operator<<(std::string &, int i);
std::string & operator<<(std::string &, unsigned i);
std::string & operator<<(std::string &, long l);
std::string & operator<<(std::string &, unsigned long l);
std::string & operator<<(std::string &, double d);
std::string & operator<<(std::string &, const char *c);
std::string & operator<<(std::string &, char c);
std::string & operator<<(std::string &, const std::string &s);
std::string & operator<<(std::string &, const std::string &s);
bool operator==(const std::string & s, const std::string & r);
bool operator==(const std::string & s, const char * r);
bool operator==(const std::string & s, const char & r);
bool operator!=(const std::string & s, const char & r);
const std::string operator+(const std::string & s, const char * r);
const std::string operator+(const std::string & s, const std::string & r);

/// Nexus related functions
bool			Abbreviates(const std::string &, const std::string &s, ncl::StringCmpEnum mode);
std::string   & ConvertBlanksToUnderscores(std::string & s);
std::string	  & ConvertToNexusToken(std::string &);
std::string	  & ConvertToNexusSingleQuoted(std::string &);
std::string	  & BlanksToUnderscores(std::string &);
void			CharBlankToUnderscore(char & c);
void			CharUnderscoreToBlank(char & c);
std::string		GetAsNexusToken(const std::string &);
std::string		GetBlanksToUnderscores(const std::string & s);
std::string		GetSingleQuoted(const std::string &) ;
std::vector<std::string>		GetVecOfPossibleAbbrevMatches(const std::string & testStr,const std::vector<std::string> & possMatches);
bool			IsCapAbbreviation(const std::string &testStr, const std::string &capAbbrev);
bool			IsLegalNexusChar(char testC);
bool			IsLegalNexusWord(const std::string &) ;
bool 			IsALegalNexusLabelForObjectN(const std::string &s, unsigned n, bool *isANumber);
bool			IsNexusPunctuation(const char c);
enum NexusQuotingReq
	{
	kNoQuotesNeededForNexus = 0, /// this enum value is kept equivalent to false
	kSingleQuotesNeededForNexus, /// punctuation or non-space whitespace characters present
	kUnderscoresSufficeForNexus  /// No nexus token-breakers
	};
NexusQuotingReq	QuotesNeeded(const std::string &) ;
bool			SetToShortestAbbreviation(std::vector<std::string> & strVec, bool allowTooShort = false);
void			ThrowIllegalNexusChar(char, const char *contextOfChar) X_SPEC_THROW(NxsException);
std::string	  & UnderscoresToBlanks(std::string &);
std::string		UpperCasePrefix(const std::string &) ;

class NxsIndexSet;
//class DblFormatter;
#if defined (__MSC_VER) 
#	define TO_UPPER toupper
#	define TO_LOWER tolower
#else
#	define TO_LOWER std::tolower
#	define TO_UPPER std::toupper
#endif
class NxsX_NotANumber {};       /* exception thrown if attempt to convert string to a number fails */
class NxsX_NumberIsTooLarge : public NxsX_NotANumber {}; //exceeds the bounds of the requested type 
class NxsX_NumberIsTooSmall : public NxsX_NotANumber {}; //smaller than the lower bound of the requested type

/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline bool IsALong(const std::string &s, long * l)
	{
	const char * c = s.c_str();
	return c && NxsString::to_long(c, l);
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline bool IsADouble(const std::string &s, double * d) 
	{
	const char * c = s.c_str();
	return c && NxsString::to_double(c, d);
	}

inline std::string  & ConvertBlanksToUnderscores(std::string & s)
	{
	std::replace(s.begin(), s.end(), ' ', '_');
	return s;
	}
	
inline std::string GetBlanksToUnderscores(const std::string & s)
	{
	std::string c(s);
	return ConvertBlanksToUnderscores(c);
	}

inline std::string & ConvertToNexusToken(std::string & s)
	{
	const NexusQuotingReq nqr = QuotesNeeded(s);
	if (nqr == kNoQuotesNeededForNexus)
		return s;
	return (nqr == kSingleQuotesNeededForNexus ? ConvertToNexusSingleQuoted(s) : ConvertBlanksToUnderscores(s));
	}
	
inline std::string GetAsNexusToken(const std::string & s)
	{
	const NexusQuotingReq nqr = QuotesNeeded(s);
	if (nqr == kNoQuotesNeededForNexus)
		return s;
	return (nqr == kSingleQuotesNeededForNexus ? GetSingleQuoted(s) : GetBlanksToUnderscores(s));
	}

inline bool operator==(const std::string & s, const std::string & r)
	{
	return StrEquals(s, r, ncl::kStringRespectCase);
	}
	
inline bool operator==(const std::string & s, const char * r)
	{
#	if defined(C_FUNCS_IN_STD_NMSPACE)
		using std::strlen;
#	endif
	return StrEquals(s, r, ncl::kStringRespectCase);
	}
	
inline bool operator==(const std::string & s, const char & r)
	{
	return (s.length() == 1 && s[0] == r);
	}

inline bool operator!=(const std::string & s, const char & r)
	{
	return (s.empty() || s[0] != r || s.length() > 1);
	}

inline const std::string operator+(const std::string & s, const char * r)
	{
	std::string result(s);
	result.append(r);
	return result;
	}

inline const std::string operator+(const std::string & s, const std::string & r)
	{
	std::string result(s);
	result.append(r);
	return result;
	}

inline std::string  & RightJustifyString(std::string &buff, const std::string &s, unsigned w )
	{
	if (s.length() >= w)
		buff.append(s);
	else
		buff.append(std::string(w - s.length(), ' ') + s);
	return buff;
	}
	
inline std::string  & RightJustifyLong(std::string &s, long x, unsigned w)
	{
	return StrPrintF(s, "%*ld", w, x);
	}

inline std::string  & RightJustifyDbl(std::string &s, double x, unsigned w, unsigned p)
	{
	return StrPrintF(s, "%*.*g", w, p, x);
	}

inline std::string MakeRightJustifiedString(const std::string &s, unsigned w)
	{
	if (s.length() >= w)
		return s;
	return std::string(w-s.length(), ' ') + s;
	}
	
inline std::string MakeRightJustifiedLong(long x, unsigned w)
	{
	std::string s;
	return RightJustifyLong(s, x, w);
	}
	
inline std::string MakeRightJustifiedDbl(double x, unsigned w, unsigned p)
	{
	std::string s;
	return RightJustifyDbl(s, x, w, p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Less than functor that can be used as comparison function when keeping case-insensitive NxsStrings in sorted 
|	containers (e.g. sets or maps)
*/
class NxsStringNoCaseLess 
  : public std::less<std::string>
	{
	public:
		bool operator()(const std::string &lOp, const std::string &rOp) const 
		{
		std::string::const_iterator lOpIt = lOp.begin();
		std::string::const_iterator rOpIt = rOp.begin();
		while (lOpIt != lOp.end() && rOpIt != rOp.end())
			{
			const char l = (char) TO_UPPER(*lOpIt++); //POL VC says toupper not in std
			const char r = (char) TO_UPPER(*rOpIt++); //POL VC says toupper not in std
			if (l != r)
				{
				if (l < r)
					return true;
				return false;
				}
			}
		return (lOpIt == lOp.end() && rOpIt != rOp.end());
		}
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	Function object (Unary Predicate functor) that stores one string. The ()(const std::string &) operator then returns 
|	true if the stored string is an abbreviation of the argument (following the rules of IsCapAbbreviation)
*/
class NStrCapAbbrevEquals 
	{
	public :
					NStrCapAbbrevEquals(const std::string &s);
		bool	    operator()(const std::string &s) const;

	protected :
		std::string compStr;
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	Provides a unary predicate functor that calls GetName() and then perform a case-insensitive StrEquals test
|	the Obj must contain the function:
|		std::string GetName() const
*/
template <class Obj> class NamedObjEqualsStored
	{
	protected :
		std::string compStr;
	public :
		NamedObjEqualsStored(const std::string &s)	
			:compStr(GetCapitalized(s))	
			{
			}
		bool operator()(const Obj &o) const
			{
			return (compStr == GetCapitalized(o.GetName()));
			}

	};

/*--------------------------------------------------------------------------------------------------------------------------
|	Function object that compares pointers to any class that provides a function with the prototype
|	std::string GetName() const
|	bool sent to the constructor controls whether comparison is case-sensitive (if true is sent) or not (false is the default)
*/	
template <class T> class NamedObjLessThan : public std::less<T *>
	{
		bool caseSensitive;
	public :
		NamedObjLessThan(bool caseSens = false) : caseSensitive(caseSens) {}
		bool operator() (const T &lOp, const T &rOp) const 
			{
			std::string lStr = lOp->GetName();
			std::string rStr = rOp->GetName();
			if (!caseSensitive)
				{
				Capitalize(lStr);
				Capitalize(rStr);
				}
			return (lStr < rStr);
			}
	};


/*--------------------------------------------------------------------------------------------------------------------------
|	return true if the std::string stored in NStrCapAbbrevEquals is an abbreviation of the argument.
*/
inline bool NStrCapAbbrevEquals::operator()(
  const std::string &s) const/* the string to be compared */
	{
	return IsCapAbbreviation(compStr, s);
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Returns a Capitalized version of the std::string.  The calling object is not altered.
|	Written for ease of use.  Simply copies the std::string, calls Capitalize() to the copy, and then returns it.
*/
inline std::string GetCapitalized(const std::string & f) 
	{
	std::string s(f);
	return ToUpper(s);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	returns true if the StrEquals comparison function is true for this and any element in the vector s.
*/
inline bool IsInVector(
  std::string & s,
  const std::vector<std::string> & v,
  ncl::StringCmpEnum mode)/* the argument of the StrEquals(std::string, ncl::StringCmpEnum) function called for every element in the vector */
	{
	for (std::vector<std::string>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt)
		{
		if (StrEquals(s, *vIt, mode))
			return true;
		}
	return false;
	}

template<>
inline std::string & AppendNumber<int>(std::string & s, const int i)
	{
	char tmp[81];
	sprintf(tmp, "%d", i);
	s.append(tmp);
	return s;
	}

template<>
inline std::string & AppendNumber<unsigned>(std::string & s, const unsigned i)
	{
	char tmp[81];
	sprintf(tmp, "%u", i);
	s.append(tmp);
	return s;
	}
template<>
inline std::string & AppendNumber<long>(std::string & s, const long l)
	{
	char tmp[81];
	sprintf(tmp, "%ld", l);
	s.append(tmp);
	return s;
	}

template<>
inline std::string & AppendNumber<unsigned long>(std::string & s, const unsigned long l)
	{
	char tmp[81];
	sprintf(tmp, "%lu", l);
	s.append(tmp);
	return s;
	}
	
template<>
inline std::string & AppendNumber<double>(std::string & s, const double d)
	{
	return StrPrintF(s, "%g", d);
	}

template<>
inline std::string & AppendNumber<float>(std::string & s, const float d)
	{
	return StrPrintF(s, "%g", d);
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	Capitalizes all lower case letters in the stored string by calling ToUpper.
*/
inline std::string & Capitalize(std::string & s)
	{
	ToUpper(s);
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the stored string is an abbreviation (or complete copy) of the supplied string `s'.
*/
inline bool Abbreviates(
  const std::string & f,
  const std::string & s,	     /* the full comparison string */
  ncl::StringCmpEnum    mode)   /* if equal to abbrev, a non-case-sensitive comparison will be made, otherwise comparison will respect case */
  	{
	if (mode == ncl::kStringCapAbbrev)
		return IsCapAbbreviation(f, s);
	return IsStdAbbreviation(f, s, mode == ncl::kStringRespectCase);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Uses the mode argument to call (and return the result of) the correct string comparison function. 
*/
inline bool StrEquals(
  const std::string & f,
  const std::string & s,	   /* the string to which *this is compared */
  ncl::StringCmpEnum mode)      /* should be one of these three: respect_case, no_respect_case or abbrev */
  	{
	switch (mode) {
		case ncl::kStringRespectCase :
			return (std::strcmp(f.c_str(), s.c_str()) == 0);
		case ncl::kStringNoRespectCase :
			return EqualsCaseInsensitive(f, s);
		case ncl::kStringCapAbbrev :
			return IsCapAbbreviation(f, s);
		default :
			PHYCAS_ASSERT(0);// incorrect setting for mode
		}
	return false;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Uses the mode argument to call (and return the result of) the correct string comparison function. 
*/
inline bool EqualsLong(
  const std::string & s,
  long l)      
  	{
	std::string f;
	AppendNumber<long>(f, l);
	return EqualsCaseInsensitive(s, f);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if `c' is any Nexus punctuation character:
|>
|	()[]{}/\,;:=*'"`-+<>
|>
*/
inline bool IsNexusPunctuation(
  const char c) /* the character in question */
	{
	return (std::strchr("()[]{}/\\,;:=*\'\"`-+<>", c) != NULL);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Creates a new string (and returns a reference to the new string) composed of the integer `i' followed by a space and
|	then the string `s'. If `i' is not 1, then an 's' character is appended to make `s' plural. For example, if `i' were 0,
|	1, or 2, and `s' is "character", then the returned string would be "0 characters", "1 character" or "2 characters", 
|	respectively. Obviously this only works if adding an 's' to the supplied string makes it plural.
*/
inline std::string &AppendNumberThenWord(
  std::string & f,
  unsigned i,		   /* the number */
  const std::string s)    /* the string needing to be pluralized */
	{
	f << i << ' ' << s;
	if (i != 1)
		f << 's';
	return f;
	}
	
inline std::string & DblFormatter::FormatDouble(std::string & s, const double d) const
	{
	return StrPrintF(s, formatCode.c_str(), d);
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber(written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  int i)	/* the integer to append */
	{
	return AppendNumber<int>(s, i);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  unsigned i)   /* the unsigned integer to append */
	{
	return AppendNumber<unsigned>(s, i);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  long l)       /* the long integer to append */
	{
	return AppendNumber<long>(s, l);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  unsigned long l)      /* the unsigned long integer to append */
	{
	return AppendNumber<unsigned long>(s, l);
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  double d)     /* the double floating point value to append */
	{
	return AppendNumber<double>(s, d);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the append operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string & operator<<(
  std::string & s,
  const char * c)	/* the C-string to append */
	{
	s.append(c);
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the += operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string & operator<<(
  std::string & s,
  char c)       /* the char to append */
	{
	s.push_back(c); //std namespace's operator+=
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the += operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & f,
  const std::string &s)   /* the std::string to append */
	{
	f.append(s); //std namespace's operator+=
	return f;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns string as a Pascal string (array of unsigned characters with the length in the first byte).
*/
inline unsigned char *p_str(
  const std::string & s,
  unsigned char * buffer)	/* buffer to receive current string in Pascal form (i.e. length in first byte) */
  	{
	std::memmove(buffer + 1, s.c_str(), s.length());
	buffer[0] = (unsigned char) s.length();
	return buffer;
	}

// ############################# start of standalone functions ##########################

/*--------------------------------------------------------------------------------------------------------------------------
|	Writes the string `s' to the ostream `out'.
*/
inline std::ostream & operator<<(
  std::ostream & out,		    /* the stream to which the string `s' is to be written */
  const std::string & s)   /* the string to write */
	{
	out << s.c_str();
	return out;
	}

inline std::string GetOrderString(
  unsigned ind)
	{
	std::string s;
	s << ind; 
	const unsigned i = ind % 10;
	if (i == 1)
		s << "st";
	else if (i == 2)
		s << "nd";
	else if (i == 3)
		s << "rd";
	else 
		s << "th";
	return s;
	}

inline void	CharToUpper(char & c)
	{
	c = (char) TO_UPPER(c); //POL VC says toupper not in std
	}

inline void	CharToLower(char & c)
	{
	c = (char) TO_LOWER(c); //POL VC says tolower not in std
	}

inline void	CharBlankToUnderscore(char & c)
	{
	if (c == ' ')
		c = '_';
	}

inline void	CharUnderscoreToBlank(char &c)
	{
	if (c == '_')
		c = ' '; 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts every character in the stored string to its lower case equivalent.
*/
inline std::string & ToLower(std::string & s)
	{
	std::for_each(s.begin(), s.end(), CharToLower);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts any blank spaces found in the stored string to the underscore character.
*/
inline std::string & BlanksToUnderscores(std::string & s)
	{
	std::for_each(s.begin(), s.end(), CharBlankToUnderscore);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts any underscore characters found in the stored string to blank spaces.
*/
inline std::string & UnderscoresToBlanks(std::string & s)
	{
	std::for_each(s.begin(), s.end(), CharUnderscoreToBlank);
	return s;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Converts every character in the stored string to its lower case equivalent.
*/
inline std::string & ToUpper(std::string & s)
	{
	std::for_each(s.begin(), s.end(), CharToUpper);
	return s;
	}

inline StringDblFormatterRef operator<<(std::string &s, DblFormatter & d)
	{
	return StringDblFormatterRef(s, d);
	}
	
inline std::string operator<<(const StringDblFormatterRef & sdf, const double d)
	{
	return sdf.dfRef.FormatDouble(sdf.strRef, d);
	}
	
inline bool IsLegalNexusWord(const std::string & s)
	{
	//@POL check for non-ascii characters here?
	return (s.length() > 1 || !IsNexusPunctuation(s[0]));
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Converts to the templated type or returns the max size (DBL_MAX, LONG_MAX, ULONG_MAX) if the entire string isn't a number
*/
template<typename T>
T ConvertEntireStr(const std::string &s);

//template<typename OUT, typename T>
//inline T ErrReportAndReturn(OUT &outStream, const char *msg, T toReturn);

//@ inlines here to prevent link errors in metrowerks
template<>
inline double ConvertEntireStr<double>(const std::string &s)
	{
	if (s.empty())	
		return DBL_MAX;
	char *endPtr;
	const char *const startPtr = s.c_str();
	double d = std::strtod(startPtr, &endPtr);
	const unsigned charRead = (unsigned) (endPtr - startPtr);
	if (charRead != s.length() || d == HUGE_VAL || d == -HUGE_VAL)
		return DBL_MAX;
	return d;
	}

template<>
inline long ConvertEntireStr<long>(const std::string &s)
	{
	if (s.empty())	
		return LONG_MAX;
	char *endPtr;
	const char *const startPtr = s.c_str();
	long l = std::strtol(startPtr, &endPtr, 10);
	const unsigned charRead = (unsigned) (endPtr - startPtr);
	if (charRead != s.length() || l == LONG_MAX || l == LONG_MIN)
		return LONG_MAX;
	return l;
	}

template<>
inline unsigned long ConvertEntireStr<unsigned long>(const std::string &s)
	{
	if (s.empty() || s[0] == '-')	
		return ULONG_MAX;
	char *endPtr;
	const char *const startPtr = s.c_str();
	unsigned long ul = std::strtoul(startPtr, &endPtr, 10);
	const unsigned charRead = (unsigned) (endPtr - startPtr);
	if (charRead != s.length() || ul == ULONG_MAX)
		return ULONG_MAX;
	return ul;
	}

template<>
inline unsigned ConvertEntireStr<unsigned>(const std::string &s)
	{
	unsigned long ul = ConvertEntireStr<unsigned long>(s);
	if (ul == ULONG_MAX || ul > (unsigned long) UINT_MAX)
		return UINT_MAX;
	return (unsigned) ul;
	}

inline bool StartsWithANumber(const std::string &s)
	{
#   if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::isdigit;
#   endif
	return (!s.empty() && std::isdigit(s[0]));
	}
	
inline bool IsSpace(char c)	//mac gcc isn't liking isspace as third arg to find_if //@ need to move
	{
#	if defined(__MWERKS__)
		return (std::isspace(c) != 0);
#	else
		return (isspace(c) != 0);
#	endif
	}
	
	

/*----------------------------------------------------------------------------------------------------------------------
|	Functor that takes two std::string arguments (`a' and `b') and returns true if and only if the length of string `a'
|	is less than the length of string `b'. Can be used with algorithms std::min_element (or std::max_element) to find
|	the shortest (or longest) string in a vector, for example.
*/
class StringLengthLess
	{
	public:
		bool operator()(const std::string & a, const std::string & b) const
			{
			return (a.length() < b.length());
			}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Class used with operator<<() to format doubles by providing the arguments needed by sprintf() to append a double 
|	value to the supplied std::string `s'. The first constructor parameter (`widthOfField') is the minimum width of the 
|	double; the second parameter (`nDigitsAfterDecimal') is the number of digits of precision to use.
*/
class DoubleFormatter
	{
	public:

						DoubleFormatter(unsigned widthOfField, unsigned nDigitsAfterDecimal);

		// Accessors
		//
		unsigned		GetWidth() const;
		unsigned		GetPrecision() const;

		// Modifiers
		//
		void			SetWidth(unsigned int widthOfField);
		void			SetPrecision(unsigned int nDigitsPrecision);

		// Utilites
		//
		std::string &	FormatDouble(std::string &s, double d);		

	protected:

		unsigned fieldWidth;					/**< is the number of characters to use when displaying the floating point value */
		unsigned digitsAfterDecimal;			/**< is the number of decimal places to show in the floating point value */
		boost::format formatter;				/**< is the formatting object */

	private:

		static const unsigned int	kMaxDoublePrintPrecision = std::numeric_limits<double>::digits10; 
		void						RecreateFormatCode();
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor that returns value of data member `fieldWidth'.
*/
inline unsigned DoubleFormatter::GetWidth() const
	{
	return fieldWidth;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor that returns value of data member `digitsAfterDecimal'.
*/
inline unsigned DoubleFormatter::GetPrecision() const
	{
	return digitsAfterDecimal;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets the value of data member `fieldWidth' to the supplied value `widthOfField', then calls 
|	RecreateFormatCode().
*/
inline void DoubleFormatter::SetWidth(
  unsigned widthOfField)	/**< is the new field width */
	{
	fieldWidth = widthOfField;
	RecreateFormatCode();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier that sets the value of data member `digitsAfterDecimal' to the supplied value `nDigitsPrecision', then
|	calls RecreateFormatCode(). Note: this function ensures that `digitsAfterDecimal' will be no larger than
|	`kMaxDoublePrintPrecision'.
*/
inline void DoubleFormatter::SetPrecision(
  unsigned nDigitsPrecision)	/**< is the new precision */
	{
	PHYCAS_ASSERT(nDigitsPrecision < fieldWidth - 1);
	digitsAfterDecimal =  (unsigned) (nDigitsPrecision > kMaxDoublePrintPrecision ? kMaxDoublePrintPrecision : nDigitsPrecision);
	RecreateFormatCode();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The first constructor parameter (`widthOfField') is the minimum width of the double; the second parameter
|	(`nDigitsAfterDecimal') is the number of digits of precision to use.
*/
inline DoubleFormatter::DoubleFormatter(
  unsigned widthOfField,		/**< is the field width */
  unsigned nDigitsAfterDecimal)	/**< is the precision */
	: fieldWidth(widthOfField), digitsAfterDecimal(UINT_MAX)
	{
	SetPrecision(nDigitsAfterDecimal); //crops precision if needed and calls RecreateFormatCode()
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the string `formatCode' sent to sprintf() function based on the values of `fieldWidth' and 
|	`digitsAfterDecimal'.
*/
inline void DoubleFormatter::RecreateFormatCode()
	{
	//@ perhaps we should be using %.g here?
	std::string formatStr;
	if (fieldWidth == UINT_MAX)
		{
		if (digitsAfterDecimal == UINT_MAX)
			formatStr = "%f"; //strcpy(formatCode, "%lf");
		else
			formatStr = str(boost::format("%%.%||f") % digitsAfterDecimal); //sprintf(formatCode, "%%.%dlf", digitsAfterDecimal);
		}
	else
		{
		if (digitsAfterDecimal == UINT_MAX)
			formatStr = str(boost::format("%%%||f") % fieldWidth); //sprintf(formatCode, "%%%dlf", fieldWidth);
		else
			formatStr = str(boost::format("%%%||.%||f") % fieldWidth % digitsAfterDecimal); //sprintf(formatCode, "%%%d.%dlf", fieldWidth, digitsAfterDecimal);
		}
	formatter = boost::format(formatStr);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline std::string &DoubleFormatter::FormatDouble(std::string &s, double d)
	{
	formatter % d;
	s.append(str(formatter));
	return s;
	}

std::string	&			append_unsigned(std::string &s, unsigned v);
std::string				abbreviate(const std::string &s, unsigned max_chars);

#if 0

template<typename T>
std::string &AppendNumber(std::string &s, const T i);

//class StringDblFormatterRef
//	{
//	public:
//		StringDblFormatterRef(std::string &s, DoubleFormatter &d) : strRef(s), dfRef(d) {}
//		std::string			&	strRef;
//		const DoubleFormatter		dfRef;
//	};

//StringDblFormatterRef	operator<<(std::string &, DoubleFormatter &);
//std::string				operator<<(StringDblFormatterRef &, const double d);

std::string &			operator<<(std::string &, int i);
std::string &			operator<<(std::string &, unsigned i);
std::string &			operator<<(std::string &, long l);
std::string &			operator<<(std::string &, unsigned long l);
std::string &			operator<<(std::string &, double d);
std::string &			operator<<(std::string &, const char *c);
std::string &			operator<<(std::string &, char c);
std::string &			operator<<(std::string &, const std::string &s);
bool					operator==(const std::string &s, const std::string &r);
bool					operator==(const std::string &s, const char *r);
bool					operator==(const std::string &s, const char &r);
bool					operator!=(const std::string &s, const char &r);
const std::string		operator+(const std::string &s, const char *r);
const std::string		operator+(const std::string &s, const std::string &r);

// moved up, out of #if 0
//std::string	&			append_unsigned(std::string &s, unsigned v);
//std::string				abbreviate(const std::string &s, unsigned max_chars);

template<>
inline std::string &AppendNumber<int>(std::string &s, const int i)
	{
	char tmp[81];
	std::sprintf(tmp, "%d", i);
	s.append(tmp);
	return s;
	}

template<>
inline std::string &AppendNumber<unsigned>(std::string & s, const unsigned i)
	{
	char tmp[81];
	std::sprintf(tmp, "%u", i);
	s.append(tmp);
	return s;
	}

template<>
inline std::string &AppendNumber<long>(std::string & s, const long l)
	{
	char tmp[81];
	std::sprintf(tmp, "%ld", l);
	s.append(tmp);
	return s;
	}

template<>
inline std::string &AppendNumber<unsigned long>(std::string & s, const unsigned long l)
	{
	char tmp[81];
	std::sprintf(tmp, "%lu", l);
	s.append(tmp);
	return s;
	}
	
template<>
inline std::string &AppendNumber<double>(std::string & s, const double d)
	{
	s = str(boost::format("%g") % d);
	return s;
	}

template<>
inline std::string & AppendNumber<float>(std::string & s, const float d)
	{
	s = str(boost::format("%g") % d);
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber(written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,	/**< is the string */
  int i)			/**< is the integer to append */
	{
	return AppendNumber<int>(s, i);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  unsigned i)   /* the unsigned integer to append */
	{
	return AppendNumber<unsigned>(s, i);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  long l)       /* the long integer to append */
	{
	return AppendNumber<long>(s, l);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  unsigned long l)      /* the unsigned long integer to append */
	{
	return AppendNumber<unsigned long>(s, l);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the AppendNumber (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & s,
  double d)     /* the double floating point value to append */
	{
	return AppendNumber<double>(s, d);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the append operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string & operator<<(
  std::string & s,
  const char * c)	/* the C-string to append */
	{
	s.append(c);
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the += operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string & operator<<(
  std::string & s,
  char c)       /* the char to append */
	{
	s.push_back(c); //std namespace's operator+=
	return s;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Another way to call the += operator (written to make it possible to use a std::string like an ostream)
*/
inline std::string &operator<<(
  std::string & f,
  const std::string &s)   /* the std::string to append */
	{
	f.append(s); //std namespace's operator+=
	return f;
	}

#endif

#endif
