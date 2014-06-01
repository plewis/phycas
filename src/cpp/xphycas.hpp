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

#ifndef PHYCAS_XPYPHY_H 
#define PHYCAS_XPYPHY_H

#include <string>
#include <exception>

// General exception class for transmitting exceptions from C++ to Python.
// The attached message should be meaningful if caught and displayed in
// a Python environment.
//
class XPyPhy : public std::exception
	{
	public:
		XPyPhy() throw() {}
		XPyPhy(const std::string s) throw()
			:msg() 
			{
			try {
				msg = s;
				}
			catch (...)
				{
				}
			}
		virtual ~XPyPhy() throw() {}
		std::string	msg;	/**< the error to display to the user */
		const char * what () const throw ()
			{
			return msg.empty() ? "Unknown exception" : msg.c_str();
			}
	}; 

#endif
