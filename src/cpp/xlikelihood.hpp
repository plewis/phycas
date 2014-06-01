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

#ifndef PYPHY_XLIKELIHOOD_HPP 
#define PYPHY_XLIKELIHOOD_HPP

#include <string>
#include <exception>

/*----------------------------------------------------------------------------------------------------------------------
|	Exception class for transmitting exceptions originating within the Likelihood module from C++ to Python. The 
|	attached message should be meaningful if caught and displayed in a Python environment.
*/
class XLikelihood : public std::exception
	{
	public:
							XLikelihood()						throw();
							XLikelihood(const std::string s)	throw();
		virtual				~XLikelihood()						throw();

		const char		*	what() const						throw ();

		std::string			msg;	/**< the error to display to the user */
	}; 

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor.
*/
inline XLikelihood::XLikelihood() throw()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor that sets string assigned to the message delivered by this exception (data member `msg').
*/
inline XLikelihood::XLikelihood(const std::string s) throw()
  :msg()
	{
	try {
		msg = s;
		}
	catch (...)
		{
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Virtual destructor.
*/
inline XLikelihood::~XLikelihood() throw()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `msg' is empty, returns string "Likelihood module"; otherwise returns string stored by `msg'. Overrides the 
|	virtual function in the base class (std::exception).
*/
inline const char *XLikelihood::what () const throw ()
	{
	return msg.empty() ? "Likelihood module" : msg.c_str();
	}

#endif
