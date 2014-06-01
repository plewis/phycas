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

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include <iostream>
#include <fstream>

#include <boost/python.hpp>

#include "ncl/nxsexception.h"
#include "ncl/nxstreesblock.h"

#include "ncl/nxscxxdiscretematrix.h"

#include "char_super_matrix.hpp"

#if defined(USING_NUMARRAY)
// These next three lines required, otherwise get link error "unresolved external symbol _PyArrayHandle"
// (at least in VC7)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif

using namespace boost::python;
//using namespace phycas;
void translate(const NxsException & e)
	{
	// Use the Python 'C' API to set up an exception object
	PyErr_SetString(PyExc_Exception, e.what());
	}

#include "ncl/nxspublicblocks.h"
#include "phycas_nexus_reader.hpp"

BOOST_PYTHON_MODULE(_ReadNexusExt)
{
	class_<NxsFullTreeDescription>("FullTreeDescription", init<NxsFullTreeDescription>())
		.def("getName", &NxsFullTreeDescription::GetName, return_value_policy<copy_const_reference>())
		.def("getNewick", &NxsFullTreeDescription::GetNewick, return_value_policy<copy_const_reference>())
		.def("isRooted", &NxsFullTreeDescription::IsRooted)
		;

	class_<std::vector<NxsFullTreeDescription> >("VecFullTreeDescription", no_init)
		.def("__iter__",  iterator<std::vector<NxsFullTreeDescription> >())
		.def("size", &std::vector<NxsFullTreeDescription>::size)
		;

	class_<PhycasNexusReader, boost::noncopyable>("NexusReaderBase", init<int>())
		.def("readFile", &PhycasNexusReader::ReadFilePath)
		.def("getNChar", &PhycasNexusReader::GetNChar)
		.def("getErrorMessage", &PhycasNexusReader::GetErrorMessage)
		.def("getTrees", &PhycasNexusReader::GetTrees, return_value_policy<copy_const_reference>())
		.def("getTaxLabels", &PhycasNexusReader::GetTaxLabels)
		.def("clear", &PhycasNexusReader::Clear)
		;

	class_<phycas::CharSuperMatrix, boost::noncopyable>("CharSuperMatrix", no_init)
		.def("getNumMatrices", &phycas::CharSuperMatrix::GetNumMatrices)
		.def("getMatrix", &phycas::CharSuperMatrix::GetMatrix, return_value_policy<manage_new_object>())
		;

	def("getLastNexusDiscreteMatrix", GetLastDiscreteMatrix, return_value_policy<manage_new_object>());

	register_exception_translator<NxsException>(&translate);
}
