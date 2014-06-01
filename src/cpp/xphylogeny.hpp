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

#ifndef PYPHY_XPHYLOGENY_HPP 
#define PYPHY_XPHYLOGENY_HPP

#include <string>
#include <exception>

/*----------------------------------------------------------------------------------------------------------------------
|	Exception class for transmitting exceptions originating within the Phylogeny module from C++ to Python. The attached
|	message should be meaningful if caught and displayed in a Python environment.
*/
class XPhylogeny : public std::exception
	{
	public:
							XPhylogeny()					throw();
							XPhylogeny(const std::string s)	throw();
		virtual				~XPhylogeny()					throw();

		const char		*	what() const					throw ();

		std::string			msg;	/**< the error to display to the user */
	}; 

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor.
*/
inline XPhylogeny::XPhylogeny() throw()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor that sets string assigned to the message delivered by this exception (data member `msg').
*/
inline XPhylogeny::XPhylogeny(const std::string s) throw()
  :msg() /*@POL Mark, did you add this section? If so, why not just assign msg here?
  			No, I didn't add this section (as far as I can remember. cvs indicates that I haven't touched this file before
  			now (though if it was moved I might have changed it elsewhere).
  			string assignment can throw a memory related exception, so that might be why the assignment is done in the 
  			try-catch below */
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
inline XPhylogeny::~XPhylogeny() throw()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `msg' is empty, returns string "Phylogeny module"; otherwise returns string stored by `msg'. Overrides the virtual
|	function in the base class (std::exception).
*/
inline const char *XPhylogeny::what () const throw ()
	{
	return msg.empty() ? "Phylogeny module" : msg.c_str();
	}

#endif
