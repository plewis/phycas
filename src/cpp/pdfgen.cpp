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

#include "pdfgen.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets ...
*/
PDFGen::PDFGen()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor clears ...
*/
PDFGen::~PDFGen()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
void PDFGen::createDocument()
	{
    std::ofstream outf("doof.pdf");
    outf << "%PDF-1. 4\n";
    outf << "1 0 obj\n";
    outf << "<< /Type /Catalog\n";
    outf << "/Outlines 2 0 R\n";
    outf << "/Pages 3 0 R\n";
    outf << ">>\n";
    outf << "endobj\n";
    outf << "2 0 obj\n";
    outf << "<< /Type /Outlines\n";
    outf << "/Count 0\n";
    outf << ">>\n";
    outf << "endobj\n";
    outf << "3 0 obj\n";
    outf << "<< /Type /Pages\n";
    outf << "/Kids [ 4 0 R ]\n";
    outf << "/Count 1\n";
    outf << ">>\n";
    outf << "endobj\n";
    outf << "4 0 obj\n";
    outf << "<< /Type /Page\n";
    outf << "/Parent 3 0 R\n";
    outf << "/MediaBox [ 0 0 612 792 ]\n";
    outf << "/Contents 5 0 R\n";
    outf << "/Resources << /ProcSet 6 0 R\n";
    outf << "/Font << /F1 7 0 R >>\n";
    outf << ">>\n";
    outf << ">>\n";
    outf << "endobj\n";
    outf << "5 0 obj\n";
    outf << "<< /Length 73 >>\n";
    outf << "stream\n";
    outf << "BT\n";
    outf << "/F1 24 Tf\n";
    outf << "100 100 Td\n";
    outf << "( Hello World ) Tj\n";
    outf << "ET\n";
    outf << "endstream\n";
    outf << "endobj\n";
    outf << "6 0 obj\n";
    outf << "[ /PDF /Text ]\n";
    outf << "endobj\n";
    outf << "7 0 obj\n";
    outf << "<< /Type /Font\n";
    outf << "/Subtype /Type1\n";
    outf << "/Name /F1\n";
    outf << "/BaseFont /Helvetica\n";
    outf << "/Encoding /MacRomanEncoding\n";
    outf << ">>\n";
    outf << "endobj\n";
    outf << "xref\n";
    outf << "0 8\n";
    outf << "0000000000 65535 f\n";
    outf << "0000000009 00000 n\n";
    outf << "0000000074 00000 n\n";
    outf << "0000000120 00000 n\n";
    outf << "0000000179 00000 n\n";
    outf << "0000000364 00000 n\n";
    outf << "0000000466 00000 n\n";
    outf << "0000000496 00000 n\n";
    outf << "trailer\n";
    outf << "<< /Size 8\n";
    outf << "/Root 1 0 R\n";
    outf << ">>\n";
    outf << "startxref\n";
    outf << "625\n";
    outf << "%%EOF\n";
    outf.close();
	}
}	// namespace phycas
