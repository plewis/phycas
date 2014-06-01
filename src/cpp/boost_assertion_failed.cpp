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

#if !defined(NDEBUG)
#include <iostream>
#include <string>
#include <boost/assert.hpp>

#define DEBUG_INPUT_ON_ASSERT
/*    #define CRASH_ON_ASSERT
*/
#if defined(BOOST_ENABLE_ASSERT_HANDLER)
void boost::assertion_failed(char const * expr, char const * function, char const * file, long line)
	{
	std::cerr << "\nBoost assertion failed:";
	std::cerr << "\n  expr: " << expr;
	std::cerr << "\n  func: " << function;
	std::cerr << "\n  file: " << file;
	std::cerr << "\n  line: " << line;
	std::cerr << std::endl;
#	if defined (DEBUG_INPUT_ON_ASSERT)
		std::cerr << "Enter a key and return to exit" << std::endl;
		std::string s;
		std::cin >> s;
#	endif 
#	if defined (CRASH_ON_ASSERT)
		// this hack works to stop gdb on an assert (even in cases in which adding 
		//	a breakpoint seems to be difficult).
		char * c = 0L;
		*c = 'q';
#	endif 
	std::exit(1);
	}
#endif 

#endif
