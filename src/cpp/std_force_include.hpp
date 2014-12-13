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

#ifndef STD_FORCE_INCLUDE_HPP
#define STD_FORCE_INCLUDE_HPP

#if defined(__cplusplus)

// This file (or one like it) should be specified as the force include file (a header file that is
// included first in every source code file in the project).

#define PYTHON_ONLY	// for code that only makes sense when exported to Python (most code can be used to construct a pure C++ program)
#define NO_IDL_TYPES
#define POLPY_NEWWAY  1
#define POLPY_OLDWAY  !(POLPY_NEWWAY)

#if defined(_MSC_VER)
// The macros below allows you to insert #pragma TODO("fix this") in C++ code and have it show up in the Visual Studio IDE
// see http://mail.python.org/pipermail/python-list/2002-January/121346.html
#   define __STRINGIZE__(x) #x
#   define TODO(x) message(__FILE__"("__STRINGIZE__(__LINE__)") : " "TODO: "#x)
#endif

#if ! defined(NDEBUG)
#	define BOOST_DEBUG_PYTHON
#	define BOOST_ALL_NO_LIB
#	define BOOST_ENABLE_ASSERT_HANDLER
#endif

#undef DECLARING_ALL_EXCEPTIONS
#if defined DECLARING_ALL_EXCEPTIONS
#	define X_SPEC_THROW(x) throw(x)
#else
#	define X_SPEC_THROW(x)
#endif

#if defined(_MSC_VER)
#	define _CRT_SECURE_NO_WARNINGS	// avoid warnings about unsafe functions (e.g. sprintf)
//#	pragma warning(disable: 4996)	// warning about loss of data when converting size_t to int
#endif

#include "phycas_config.h"

#endif	/* defined(__cplusplus) */
#endif	/* STD_FORCE_INCLUDE_HPP */
