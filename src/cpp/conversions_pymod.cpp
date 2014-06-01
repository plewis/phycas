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

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include <boost/python.hpp>

//#include "phycas/force_include.h"
#include "states_patterns.hpp"	// for int8_t definition
#include "basic_tree.hpp"
//#include <boost/python/tuple.hpp>
//#include <boost/python/numeric.hpp>
#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
#include "thirdparty/pyconversions.h"	// from HippoDraw
using namespace boost::python;

BOOST_PYTHON_MODULE(_ConversionsExt)
{
#if defined(USING_NUMARRAY)
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");
	//numeric::array::set_module_and_type();
#endif

	// these lines taken from HippoGraph
	std_vector_to_tuple<unsigned>();
	std_vector_to_tuple<int8_t>();
	std_vector_to_tuple<int>();
	std_vector_to_tuple<double>();
	std_vector_to_tuple<std::string>();
	std_vector_to_tuple<phycas::TreeNode *>();

	from_python_sequence<std::vector<unsigned>, variable_capacity_policy>();
    from_python_sequence<std::vector<int8_t>, variable_capacity_policy>();
    from_python_sequence<std::vector<int>, variable_capacity_policy>();
    from_python_sequence<std::vector<double>, variable_capacity_policy>();
    from_python_sequence<std::vector<std::string>, variable_capacity_policy>();
}
