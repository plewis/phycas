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

#if ! defined(STATES_PATTERNS_HPP)
#define STATES_PATTERNS_HPP

#include <vector>
#include <map>
#include <set>
#include <list>
#include <string>

// some basic types derived from STL containers
typedef std::vector< char >				char_vect_t;
typedef std::vector< int >				int_vect_t;
typedef std::vector< unsigned int >		uint_vect_t;
typedef std::vector< float >			float_vect_t;
typedef std::vector< double >			double_vect_t;
typedef std::vector< double_vect_t >	double_vect_vect_t;
typedef std::set< int >					int_set_t;
typedef std::set< unsigned int >		uint_set_t;
typedef std::list< unsigned int >		uint_list_t;
typedef std::vector< std::string >		string_vect_t;

#if HAVE_INTTYPES_H
#	include <inttypes.h>
#elif HAVE_STDINT_H
#	include <stdint.h>
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#	include <basetsd.h>
	typedef   INT8 int8_t;
	typedef  UINT8 uint8_t;
	typedef  INT64 int64_t;
	typedef UINT64 uint64_t;
#elif defined(_MSC_VER)
	typedef signed char int8_t;
	typedef unsigned char uint8_t;
	typedef long long int64_t;
	typedef unsigned long long uint64_t;
#elif defined(_WIN32)
#	include <stdint.h>
#endif

typedef std::vector< int8_t >						int8_vect_t;
typedef std::vector< int8_t >						state_list_t;
typedef std::vector< state_list_t >					state_list_vect_t;
typedef  std::vector< unsigned int >				state_list_pos_t;
typedef std::vector< state_list_pos_t >				state_list_pos_vect_t;
typedef int8_t 										state_code_t;
typedef std::vector< int8_t >						pattern_t;
typedef double										pattern_count_t;
typedef	std::map<int8_vect_t, double>				pattern_map_t;
typedef std::vector< pattern_map_t >                pattern_map_vect_t;
typedef	std::vector<int8_vect_t>					pattern_vect_t;
typedef std::vector< uint_vect_t > 					pattern_to_sites_t;//UINT_LIST
typedef std::map<int8_vect_t, uint_vect_t>			pattern_to_sites_map_t;//UINT_LIST
typedef	std::vector<pattern_count_t>				count_vect_t;
typedef std::vector< int8_t >                       StateMapping;

#endif


