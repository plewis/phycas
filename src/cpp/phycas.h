/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#if !defined(PHYCAS_H)
#define PHYCAS_H

#include <cfloat>

#include "ncl/nxsdefs.h"
	//USING_VARIANCE_FOR_RATEHET as opposed to using the shape parameter of the gamma distribution for rate het.
#define USING_VARIANCE_FOR_RATEHET
	//	If USING_NODE_NUMBERS_ONLY is defined, tree Nodes do not have a std::string data member that stores their name
	//	an external class must be used to translate from the taxon number to a name.
#define USING_NODE_NUMBERS_ONLY
typedef char NStateInt;	// to be used in code where speed is critical, but we need a integer type big enough to hold the max num states
						// signed because - numbers are used to indicate ambiguity
#define NSTATEINT_MAX CHAR_MAX

	//	USING_PROGRAMMER_SPECIFIC_DEFINES added so that we can define 
	//	macros in our projects without touching phorest.h (I keep forgetting that I've tweaked
	//	this file and committing changes).
	// 	the !defined(USING_PROGRAMMER_SPECIFIC_DEFINES) section below provides the 
	//	stable defaults.
	//
#if !defined(USING_PROGRAMMER_SPECIFIC_DEFINES)
#		define PHYLODIVERSITY_MODULE
#endif

#define SCM_MODULE

#endif
