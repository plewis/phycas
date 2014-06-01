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

#ifndef PHYCAS_CONFIG_H
#define PHYCAS_CONFIG_H

#if defined(PYTHON_ONLY)
#	include <Python.h>
//#	include <Python/Python.h>	// this is probably better for OS X, but we'll just use an appropriate search path within the framework for now
#endif

#if defined(NDEBUG)
// no asserts in release version
#	define PHYCAS_ASSERT(expr)
#else
// debug versions use boost's assert, making it easy to set a breakpoint in
// boost_assertion_failed.cpp and inspect the call stack
#	include <boost/assert.hpp>
#	define PHYCAS_ASSERT(expr)  if (!(expr)) boost::assertion_failed((const char *)#expr, (const char *)__FUNCTION__, __FILE__, __LINE__)
#	define IGNORE_NXS_ASSERT	// used in nxs_defs.hpp to cause NXS_ASSERT to expand to nothing
#endif

#define NXS_SUPPRESS_OUTPUT

#if defined (__APPLE__)
#   define NXS_USE_UNIX_DIR_STYLE
#elif defined(_WIN32)
#   define NXS_USE_WINDOWS_DIR_STYLE
#else
#   define NXS_USE_UNIX_DIR_STYLE
#endif

///////////////////////////////// just replace these with their definition ///////////////////////////////
#define STATIC_CONST static
#define STATIC_DATA static
#define STATELESS_FUNC static

///////////////////////////////// everything below here should go ///////////////////////////////
#if 0

#if defined(__MWERKS__) && defined(__APPLE__)
#	include <MSLCarbonPrefix.h>
#endif

/* to make cipres services supplied by phycas reentrant, MTH is on a crusade against statics
The following defines will help me wade through harmless, class-level functions and
  and static data and functions that hinder thread-safety.
*/
 /* STATIC_DATA - for use with static data members and local variables that are declared static */
 /* STATIC_DATA_FUNC -  for use to declare static functions that rely on static data */
#define STATIC_DATA_FUNC static
 /* STATELESS_FUNC - for functions that are part of a class's namespace, but are static AND do not
   access in static data */

  /* STATIC_SINGLETON -  used for singletons that are declared as statics: */
#define STATIC_SINGLETON static
  /* STATIC_CONST_ACCESSOR - used to declare functions that refer to static data that is constant
  these functions will be thread-safe despite using static data */
#define STATIC_CONST_ACCESSOR static
  /* STATIC_CONST declares data as const and static */
#define STATIC_CONST static


#if defined (CORBA_CLIENT_PHYCAS) || defined(CORBA_SERVER_PHYCAS)
#	define CORBA_PHYCAS
#endif

	// if this is just a CIPRES module, then we don't read PHYCAS blocks
#if !defined (CIPRES_USING_PHYCAS_CODE_BASE)
#	define READ_PHYCAS_BLOCK
#endif

#if defined (__APPLE__)
#	define MAC_PHOREST
#	define MACH_O_PHOREST
#elif defined(_WIN32)
#	define WIN_PHOREST
#	define _WIN32_WINNT   0x0501	// assuming target OS is Windows XP or later, see phycas/misc/stopwatch.hpp
#else
#	define LINUX_PHOREST	//POL added
#endif

	// try to read _run.nex
#define DEBUGGING_FROM_RUN_NEX
#define STOP_DEBUGGER_ON_TRIPPED_ASSERTS 1

#if defined(__MWERKS__)
#	define C_FUNCS_IN_STD_NAMESPACE
#	define HAVE_PRAGMA_UNUSED
#	if defined(MAC_PHOREST)
#		define BOOST_DISABLE_THREADS
#	endif
#	if __ide_target("d_phycas_mac_console") || __ide_target("d_PhycasWinConsole") || __ide_target("d_phycas_mac_socket")|| __ide_target("mth_dev_mac_console")|| __ide_target("dummy")
#		undef NDEBUG
#	else
#		define NDEBUG
#	endif
#	if __ide_target("d_phycas_mac_memchk")
#		define MONITORING_ALLOCATION	// Special target for instantiating memory checker (it really slows things down, and so it shouldn't be on for all debugging
#	endif
#endif	/* defined(__MWERKS__) */

#if defined (WIN_PHOREST)
	//@@ this misspecification of the directory style is a nxsfilepath bug that needs to be fixed
	//#	define NXS_USE_WINDOWS_DIR_STYLE
#	define NXS_USE_WINDOWS_DIR_STYLE
#elif defined (LINUX_PHOREST)
#	define NXS_USE_UNIX_DIR_STYLE
#elif defined (MAC_PHOREST)
#	if defined (MACH_O_PHOREST)
#		define NXS_USE_UNIX_DIR_STYLE
#	else
#		define NXS_USE_MAC_DIR_STYLE
#	endif
#else
#	error need to define a platform
#endif

class Phylodiversity;
class DistributionManager;
class PhoSampler;
class PhoMCMC;
class PWD;
class LongOperationManager;


#include <vector>
typedef std::vector<double> VecDbl;
#include "phycas.h"

#endif

#endif	/* __PHOREST_CONFIG_H */
