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

// NOTE: This file is currently (19-Nov-2007) not being used in Phycas. All PDF generation is done in
// pure Python. It is still here, however, because it is anticipated that at some point we will want to
// be able to generate PDF files from the C++ side directly, at which point the Python code currently
// in _PDFGenerator.py will need to be translated to C++. With the stubs pdfgen.cpp, pdfgen.hpp and
// pdfgen_pymod.cpp in place, along with (currently-disabled) code in Jamroot and the (currently-unused)
// VC project named pdfgen, it will not be (as) hard to get C++ PDF generation working.

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include <boost/python.hpp>

#include "pdfgen.hpp"
using namespace boost::python;
using namespace phycas;

BOOST_PYTHON_MODULE(_PDFGen)
{
	class_<PDFGen>("PDFGenBase")
		.def("createDocument", &PDFGen::createDocument)
		;
}
