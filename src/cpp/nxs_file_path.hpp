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
#if ! defined(NXS_FILE_PATH_HPP)
#define NXS_FILE_PATH_HPP



#include <algorithm>

/*----------------------------------------------------------------------------------------------------------------------
|	Similar to the Perl split function.  Copies the original container of T objects into the list outList as
|	smaller containers broken at the specified value (which is omitted).
|	e.g. if the incoming string "usr/home/mine" and the splitAtVal was the char '/', then the resulting list would
|	have three strings "usr" "home" and "mine"
*/
template <typename T, class ORIG_CONTAINER>
void split(const ORIG_CONTAINER &origContainer, T splitAtVal, std::list<ORIG_CONTAINER> *outList)
	{
	typedef typename ORIG_CONTAINER::const_iterator OrigCIt;
	OrigCIt begIt = origContainer.begin();
	const OrigCIt endIt = origContainer.end();
	if (begIt == endIt)
		return;	//empty container sent
	OrigCIt copyToIt;
	do	{
		copyToIt = std::find(begIt, endIt, splitAtVal);
		if (begIt == copyToIt)
			outList->push_back(ORIG_CONTAINER());
		else
			{
			outList->push_back(ORIG_CONTAINER(begIt,copyToIt));
			begIt = copyToIt;
			}
		if (begIt != endIt)	//POL-31Dec2008
			++begIt;
		}
	while (copyToIt != endIt);
	}

#if defined(_MSC_VER)
#	define NCL_STAT				_stat
#	define NCL_STAT_STRUCT		_stat
#	define NCL_DIRMODE			0
#	define NCL_MODE_T			unsigned
#	define NCL_MKDIR(a,b)		_mkdir(a)
#	define NCL_S_ISDIR(a)		((_S_IFDIR & (a)) != 0)
#	define NCL_TIME_T			time_t
#	include <direct.h>
#else
#	define NCL_STAT				stat
#	define NCL_STAT_STRUCT		stat
#	define NCL_DIRMODE			(S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP| S_IROTH | S_IXOTH)
#	define NCL_MODE_T			mode_t
#	define NCL_MKDIR(a,b)		mkdir(a,b)
#	define NCL_S_ISDIR(a)		S_ISDIR(a)
#	define NCL_TIME_T			std::time_t
#endif


// enum CmdResult
// 	{
// 	kCmdSucceeded,
// 	kCmdFailedSilent,
// 	kCmdFailedGenerateMessage
// 	};
//
enum DirDescriptionFormat  // names for the 3 systems of interpretting strings as paths
	{
	kUNIXDirStyle,
	kMacDirStyle,
	kDOSDirStyle
	};

class NxsFilePath
	{
	public:
		NxsFilePath(const std::string &, bool isDir = false);

		bool				Exists() const;
		std::string 		GetFileName() const;
		std::string 		GetFullName() const;
		const std::string & GetNative() const;
		NxsFilePath			GetParent() const;
		bool 				IsDirectory() const;
		bool				IsFile() const;
		void				SetQueryUserOnError(bool v) {queryUserOnError = v;}

		NxsFilePath & operator += (const NxsFilePath &);
		NxsFilePath & operator += (const std::string &);

		enum SourceOfFileName
			{
			kProgrammer = 0x01,		/* file path is not being used as command option */
			kDefaultOpt  = 0x02,		/* file path is controlled by a command option and the user has never entered anything */
			kUserSupplied = 0x03		/* file path was set by the user */
			};
		SourceOfFileName sourceOfName;
		bool				WasSetByUser() const {return (sourceOfName == kUserSupplied);}

		enum QueryResult
			{
			kCouldNotQuery,
			kProblemWasFixed,
			kUserCancelled
			};

		STATIC_DATA DirDescriptionFormat dirStringInterpretStyle; /* specifies the format that will be used to interpret strings as path descriptions */
	protected:
		enum OpenOutputReturnCode
			{
			kBothReplaceAndAppend,
			kOpenAndReplacing,
			kOpenAndAppending,
			kNeitherReplaceOrAppend,
			kOpenedNewFile,
			kCouldNotOpen
			};
		void			 	CreateNative() const;
		bool			 	OpenInput(std::ifstream &outF) const;
		bool				ReadNewPathString(const std::string &);
		void 				TranslateDirStack(const std::list<std::string> &dirStack) const;

		STATELESS_FUNC void			ShortenDirStack(std::list<std::string> &dirStack);

		bool				isAbsolute;
		bool				isDirty;
		std::string			path;
		std::string			pathAsEntered;
		mutable std::string	nativePath;
		bool				queryUserOnError;

	};


class NxsInFilePath : public NxsFilePath
	{
	public :
		NxsInFilePath() : NxsFilePath(0) {}
		NxsInFilePath(const std::string &s) : NxsFilePath(s, 0) {}

		bool		Open(std::ifstream &outF, std::string purposeOfFile = std::string());

	};
bool IsUpDirNotation(const std::string d);

inline const std::string &NxsFilePath::GetNative() const
	{
	if (isDirty)
		CreateNative();
	return nativePath;
	}
inline std::string NxsFilePath::GetFullName() const
	{
	return GetNative();
	}

inline bool IsUpDirNotation(const std::string d)
	{
	return (d.length() == 2 && d[0] == '.' && d[1] == '.');
	}
#endif
