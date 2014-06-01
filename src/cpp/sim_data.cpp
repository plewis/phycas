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

#include <cmath>
#include <fstream>
#include <algorithm>
#include <functional>
#include <boost/format.hpp>
#include "basic_tree.hpp"
#include "cond_likelihood.hpp"
#include "cond_likelihood_storage.hpp"
#include "tip_data.hpp"
#include "sim_data.hpp"
#include "phycas_string.hpp"

const int8_t phycas::SimData::missing_state = -1;	// GCC 3.2.3 (Red Hat Linux 3.2.3-20) requires initializing here rather than in header
using std::exp;

class LookupStateSymbol : public std::unary_function<int8_t, char>
	{
	public:
		LookupStateSymbol(const std::vector<std::string> & state_symbols) : symbols(state_symbols) {}

		char operator()(int8_t i)
			{
			PHYCAS_ASSERT((unsigned)i < (unsigned)symbols.size());
			return symbols[i][0];
			}

	private:
		const std::vector<std::string> & symbols;
	};

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the vector of counts (`binv') that is created in the calctBinned function. The table below shows the
|   identity of each bin for the case of 4 states. The length of `binv' is 2k - 1, where k is the number of states.
|   The number of states can thus be obtained as (sz + 1)/2, where sz is the length of `binv'.
|>
|	Count   Bin description
|   ----------------------------------------
|	n_0     all patterns containing only A
|	n_1     all patterns containing only C
|	n_2     all patterns containing only G
|	n_3     all patterns containing only T
|	n_4     patterns containing any 2 states
|	n_5     patterns containing any 3 states
|	n_6     patterns containing any 4 states
|   ----------------------------------------
|>
*/
std::vector<double> SimData::getBinnedCounts()
	{
    if (binv.empty())
        {
        std::cerr << "Doh!" << std::endl;
        }
    PHYCAS_ASSERT(!binv.empty());

    return binv;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Builds up the `binv' vector by classifying stored site patterns into the following categories (bins):
|>
|	Count   Bin description
|   ----------------------------------------
|	n_0     all patterns containing only A
|	n_1     all patterns containing only C
|	n_2     all patterns containing only G
|	n_3     all patterns containing only T
|	n_4     patterns containing any 2 states
|	n_5     patterns containing any 3 states
|	n_6     patterns containing any 4 states
|   ----------------------------------------
|>
|	Warning: this function currently assumes no missing data! A missing state is treated as if it were one of the
|	other states in the pattern. For example, the count for the pattern 'AAACA?GAA' would be stuffed into the 3-state
|	bin, with the implicit assumption being that the ? equals either A, C or G.
*/
void SimData::buildBinVector(
  //unsigned nstates,   /**< is the number of states */ //POLOLD
  unsigned nstates,   /**< is the number of states */
  bool minbins)     /**< if true, minimum number of bins will be used; if false, a bin will be created for every permutation of every possible subset of states */
	{
	PHYCAS_ASSERT(nstates == 4);    // Currently assumes DNA or RNA sequence data
	unsigned nbins = 0;
    if (minbins)
        nbins = 7;
    else
        nbins = 15;

    binv.clear();
	binv.resize(nbins, 0.0);
	std::set<int8_t> state_set;

    std::vector<bool> seen(4, false);

    // pattern_map_t associates int8_vect_t keys (pattern) with double values (count)
	for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		const int8_vect_t & p = it->first;

        seen.assign(4, false);

        // For this pattern, count distinct states by creating a set
		state_set.clear();
		unsigned last_state = UINT_MAX;
		for (int8_vect_t::const_iterator pit = p.begin(); pit != p.end(); ++pit)
			{
			int8_t curr_state = *pit;
			int cs = (int)curr_state;
			if (cs >= 0 && cs < (int)nstates)
				{
				// state is not a gap (-1), missing (nstates), or an ambiguity code (> nstates), so add to set
				state_set.insert(curr_state);
				last_state = cs;

                seen[cs] = true;
				}
			}
		unsigned sz = (unsigned)state_set.size();
		PHYCAS_ASSERT(sz > 0);
		PHYCAS_ASSERT(sz <= nstates);

		double this_count = (double)(it->second);
		if (sz == 1)
			{
			// pattern had only one state, so add pattern count to appropriate constant site bin
			binv[last_state] += this_count;
			}
		else if (minbins)
			{
			// pattern had sz states, so add pattern count to appropriate variable site bin
			binv[nstates + sz - 2] += this_count;
			}
        else
            {
            // sz > 1 and not minbins, so need to use composition of state_set to determine bin index
            if (sz == 2)
                {
                if (seen[0] && seen[1])
                    binv[4] += this_count;
                else if (seen[0] && seen[2])
                    binv[5] += this_count;
                else if (seen[0] && seen[3])
                    binv[6] += this_count;
                else if (seen[1] && seen[2])
                    binv[7] += this_count;
                else if (seen[1] && seen[3])
                    binv[8] += this_count;
                else if (seen[2] && seen[3])
                    binv[9] += this_count;
                else
                    PHYCAS_ASSERT(0);
                }
            else if (sz == 3)
                {
                if (seen[0] && seen[1] && seen[2])
                    binv[10] += this_count;
                else if (seen[0] && seen[1] && seen[3])
                    binv[11] += this_count;
                else if (seen[0] && seen[2] && seen[3])
                    binv[12] += this_count;
                else if (seen[1] && seen[2] && seen[3])
                    binv[13] += this_count;
                else
                    PHYCAS_ASSERT(0);
                }
            else
                {
                PHYCAS_ASSERT(sz == 4);
                binv[14] += this_count;
                }

            }

		}
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the binned t function used in the Gelfand-Ghosh measure for the patterns currently stored in this
|	object. The binned t function looks like this for multinomial data, 4 DNA states and 4 taxa:
|>
|	Count   Bin description                     t calculation
|   ----------------------------------------------------------------------------------------------
|	n_0     all patterns containing only A      t_0 = (n_0 + eps)/(n + 1) log[(n_0 + eps)/(n + 1)]
|	n_1     all patterns containing only C      t_1 = (n_1 + eps)/(n + 1) log[(n_1 + eps)/(n + 1)]
|	n_2     all patterns containing only G      t_2 = (n_2 + eps)/(n + 1) log[(n_2 + eps)/(n + 1)]
|	n_3     all patterns containing only T      t_3 = (n_3 + eps)/(n + 1) log[(n_3 + eps)/(n + 1)]
|	n_4     patterns containing any 2 states    t_4 = (n_4 + eps)/(n + 1) log[(n_4 + eps)/(n + 1)]
|	n_5     patterns containing any 3 states    t_5 = (n_5 + eps)/(n + 1) log[(n_5 + eps)/(n + 1)]
|	n_6     patterns containing any 4 states    t_6 = (n_6 + eps)/(n + 1) log[(n_6 + eps)/(n + 1)]
|   ----------------------------------------------------------------------------------------------
|   n = n_0 + n_1 + ... + n_6                   t = t_0 + t_1 + ... + t_6
|>
|	where n_i is the count of the number of characters placed into bin i, n is the total number of characters, and
|	eps equals 1/7 (the inverse of the total number of bins). The total number of bins is (2*nstates - 1), regardless
|	of the number of taxa, which allows this approach to Gelfand-Ghosh to scale nicely to large problems. The drawback,
|	of course, is that information about frequencies of individual rare patterns is not used.
|
|	Warning: this function currently assumes no missing data! A missing state is treated as if it were one of the
|	other states in the pattern. For example, the count for the pattern 'AAACA?GAA' would be stuffed into the 3-state
|	bin, with the implicit assumption being that the ? equals either A, C or G.
*/
double SimData::calctBinned(
  unsigned nstates, /**< is the number of states */
  bool minbins)     /**< if true, minimum number of bins will be used; if false, a bin will be created for every permutation of every possible subset of states */
	{
	PHYCAS_ASSERT(nstates == 4);    // Currently assumes DNA or RNA sequence data

	// m is the number of distinct patterns
	//double m					= (double)sim_pattern_map.size();

	// s is the number of character states
	//double s					= (double)nstates;
	//double log_s				= std::log(s);

	// n is the number of characters (i.e. sum of counts of all distinct patterns)
	double n					= (double)total_count;
	double n_plus_one			= n + 1.0;
	double log_n_plus_one		= std::log(n_plus_one);

    unsigned nbins = 0;
    if (minbins)
        nbins = 7;
    else
        nbins = 15;

	// epsilon is the number of bins (nstates constant bins plus nstates - 1 additional bins for patterns with 2, 3, ..., nstates states)
	double epsilon				= 1.0/(double)nbins;

	// classify patterns and build up bin, the vector of bin counts
    buildBinVector(nstates, minbins);

	// now accumulate t by summing over bins
	double t = 0.0;
	for (std::vector<double>::iterator vit = binv.begin(); vit != binv.end(); ++vit)
		{
		double count				= (*vit);
		double count_plus_epsilon	= count + epsilon;
		double first				= count_plus_epsilon;
		double second				= std::log(count_plus_epsilon);
		double this_term			= first*second;
		t += this_term;
		}
    t /= n_plus_one;
    t -= log_n_plus_one;

	return t;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the t function used in the Gelfand-Ghosh measure for the patterns currently stored in this object. The
|	t function looks like this for multinomial data and 4 taxa:
|>
|	sum_{i=1}^{256} (x_i + eps)/(n + 1) log[(x_i + eps)/(n + 1)]
|>
|	where 256 is the total number of possible data patterns (=nstates^ntaxa), x_i is the count of the number of
|	characters having pattern i, n is the total number of characters, and eps equals 1/256 (the inverse of the total
|	number of possible patterns).
*/
double SimData::calct(
  unsigned nstates) /**< is the number of states */
	{
	// m is the number of distinct patterns
	double m					= (double)sim_pattern_map.size();
	//unused:  double log_m				= std::log(m);

	// s is the number of character states
	double s					= (double)nstates;
	double log_s				= std::log(s);

	// n is the number of characters (i.e. sum of counts of all distinct patterns)
	double n					= (double)total_count;

	double n_plus_one			= n + 1.0;
	double log_n_plus_one		= std::log(n_plus_one);

	// ntaxa is the number of taxa in the tree
	double ntaxa				= (double)pattern_length;
	double ntaxa_times_log_s	= ntaxa*log_s;

	// epsilon is the inverse of the total number of possible patterns
	double epsilon				= std::exp(-ntaxa_times_log_s);

	// added_on is the sum of all terms in t in which the pattern count equals zero
	double log_term				= -ntaxa_times_log_s - log_n_plus_one;
	double added_on				= (1.0 - m*epsilon)*log_term/n_plus_one;

	// now compute the sum of all terms in t in which the pattern count is greater than zero
	double sum					= 0.0;
	for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
		{
		double count				= (double)it->second;
		double count_plus_epsilon	= count + epsilon;
		double first				= count_plus_epsilon/n_plus_one;
		double second				= std::log(count_plus_epsilon) - log_n_plus_one;
		double this_term			= first*second;
		sum += this_term;
		}

	double t = sum + added_on;

	return t;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns simulated data stored in `sim_pattern_map' as a string in the form of a two-column table. The first column
|	is labeled "Count" and the second is labeled "Pattern". The Count column shows the number of times its associated
|	pattern was inserted using the insertPattern function. The Pattern column shows a representation of the pattern
|	itself, using symbols for states provided in the `state_symbols' argument. The `state_symbols' argument should be
|	a vector of single-character strings supplying a symbol to represent each state that might show up in any pattern.
|	Assumes that no state in any pattern stored in `sim_pattern_map' is greater than or equal to the length of the
|	`state_symbols' vector (because states are used as indices into `state_symbols').
*/
std::string SimData::patternTable(
  const StringVect & state_symbols) /**< is a vector of strings representing states (e.g. {"A", "C", "G", "T"}). Note that each state symbol should be a string of length 1 (i.e. a single character) */
	{
	PHYCAS_ASSERT(state_symbols.size() > 0);

	outstr.clear();

	if (sim_pattern_map.empty())
		{
		outstr = "Sorry, no patterns are stored";
		}
	else
		{
		outstr = "     Count  Pattern";

		for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
			{
			// Output the count first
			outstr << str(boost::format("\n%10.1f") % it->second) << "  ";

			// Now output the pattern
			std::transform(it->first.begin(), it->first.end(), std::back_inserter(outstr), LookupStateSymbol(state_symbols));
			}
		}
	return outstr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns row `i' of `sim_pattern_vect' as a vector. Returns vector by value.
*/
std::vector<unsigned> SimData::getPatternVectRow(
  unsigned i) const
    {
    unsigned nchar = (unsigned)sim_pattern_vect.size();
    std::vector<unsigned> v(nchar, 5);
    for (unsigned k = 0; k < nchar; ++k)
        v[k] = (unsigned)sim_pattern_vect[k][i];
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Saves simulated data stored in `sim_pattern_map' to a file named `filename'. If the file already exists, it will be
|	overwritten without warning. If the file does not yet exist, it will be created. The file written will be a valid
|	NEXUS data file suitable for executing in phylogenetic analysis software that reads the NEXUS file format. The
|	supplied `taxon_names' will be used in the matrix command of the NEXUS file to specify the names of the taxa.
|	Assumes that the number of elements in `taxon_names' equals `pattern_length'. The `datatype' argument should be the
|	correct NEXUS datatype (e.g. "dna", "standard") for the data simulated. The symbols used for the states are
|	supplied in the `state_symbols' vector. Each element of this vector should be a single-character string. Assumes
|	that no state in any pattern stored in `sim_pattern_map' is greater than or equal to the length of the
|	`state_symbols' vector (because states are used as indices into `state_symbols').
*/
void SimData::saveToNexusFile(
  const std::string filename,			/**< is the name of the file to create containing the simulated data stored in this object */
  const StringVect & taxon_names,		/**< is a vector containing the taxon names to use in the saved file */
  const std::string datatype,			/**< is a string to be used as the NEXUS datatype (e.g. "dna" or "standard") */
  const StringVect & state_symbols)		/**< is a vector of strings representing states (e.g. {"A", "C", "G", "T"}). Note that each state symbol should be a string of length 1 (i.e. a single character) */
	{
	//std::cerr << "taxon_names size = " << taxon_names.size() << std::endl;
	//std::cerr << "pattern_length   = " << pattern_length << std::endl;
	//std::cerr << "taxon_names: |";
	//std::copy(taxon_names.begin(), taxon_names.end(), std::ostream_iterator<std::string>(std::cerr, "|"));

	PHYCAS_ASSERT(state_symbols.size() > 0);
	PHYCAS_ASSERT(taxon_names.size() == pattern_length);

	// Find length of longest string in taxon_names vector; this is used later for formatting purposes
	// The 2 is included in case apostrophes are needed when taxon names are output in the matrix command
	unsigned length_of_longest_name = 2 + (unsigned)std::max_element(taxon_names.begin(), taxon_names.end(), StringLengthLess())->length();

	std::ofstream outf(filename.c_str());

	outf << "#nexus" << "\n\n";
	outf << "begin data;" << "\n";
	outf << str(boost::format("  dimensions ntax=%d nchar=%d;") % pattern_length % (unsigned)total_count) << "\n";
	outf << "  format datatype=" << datatype << ";" << "\n";
	outf << "  matrix" << "\n";

	// Create a format string to use with boost::format that left-justifies (the "-" flag) the
	// taxon names in a field of width length_of_longest_name
	std::string fmtstr = str(boost::format("    %%-%ds") % length_of_longest_name);

	for (unsigned i = 0; i < pattern_length; ++i)
		{
		if (taxon_names[i].find(' ') != std::string::npos)
			{
			std::string s = "'";
			s += taxon_names[i];
			s += "'";
			outf << str(boost::format(fmtstr) % s) << "  ";
			}
		else
			{
			outf << str(boost::format(fmtstr) % taxon_names[i]) << "  ";
			}

#if 1
        // Spit out characters in the order in which they were simulated. While this is a nice feature,
        // it currently requires storing the data twice (sim_pattern_vect and sim_pattern_map)
        unsigned nchar = (unsigned)sim_pattern_vect.size();
        for (unsigned k = 0; k < nchar; ++k)
            {
            unsigned j = (unsigned)sim_pattern_vect[k][i];
            char s = state_symbols[j][0];
            outf << s;
            }
#else
		for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
			{
			// The first member of it is the state, which must be converted from its coded form
			// to the standard symbol for this data type (e.g. A, C, G or T for DNA characters)
			int8_t j = (*it).first[i];

			PHYCAS_ASSERT(j < (int8_t)state_symbols.size());
			char s = state_symbols[j][0];	// use the first (and hopefully only) character in the string at position j

			// The second member of it is the pattern count
			unsigned n = (unsigned)it->second;	//@POL assuming counts not fractional
			for (unsigned k = 0; k < n; ++k)
				{
				outf << s;
				}
			}
#endif
		outf << "\n";
		}

	outf << "  ;" << "\n";
	outf << "end;" << "\n";

	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `tmp_pattern' to `sim_pattern_map'. Increments total_count to reflect the new total number of counts over all
|   patterns. After calling, you should pass `missing_state' to the wipePattern() function to fill `tmp_pattern' with
|   invalid values. Unlike insertPattern, insertPatternOnly does not add site indices to `pattern_to_sites_map' and does
|   not add `tmp_pattern' to `sim_pattern_vect', so only use this function if you do not care to keep track of the
|   order in which sites were simulated.
*/
void SimData::insertPatternOnly(
  pattern_count_t count)	/**< is the count to be associated with the pattern now stored in `tmp_pattern' */
	{
	// insert the pattern stored in tmp_pattern into sim_pattern_map
	// Add tmp_pattern to sim_pattern_map if it has not yet been seen, otherwise increment the count
	// for this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
	pattern_map_t::iterator lowb = sim_pattern_map.lower_bound(tmp_pattern);
	if (lowb != sim_pattern_map.end() && !(sim_pattern_map.key_comp()(tmp_pattern, lowb->first)))
		{
		// Pattern is already in sim_pattern_map, so just modify its count
		lowb->second += count;
		}
	else
		{
		// tmp_pattern has not yet been stored in sim_pattern_map
		sim_pattern_map.insert(lowb, pattern_map_t::value_type(tmp_pattern, count));
		}

    total_count += count;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `tmp_pattern' to `sim_pattern_map'. Increments total_count to reflect the new total number of counts over all
|   patterns. After calling, you should pass `missing_state' to the wipePattern() function to fill `tmp_pattern' with
|   invalid values.
*/
void SimData::insertPattern(
  const uint_vect_t & sitelist, /**< is a list of site indices corresponding to the current pattern */  //UINT_LIST
  pattern_count_t count)	/**< is the count to be associated with the pattern now stored in `tmp_pattern' */
	{
	// Insert the pattern stored in tmp_pattern into sim_pattern_map
	// Add tmp_pattern to sim_pattern_map if it has not yet been seen, otherwise increment the count
	// for this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
	pattern_map_t::iterator lowb = sim_pattern_map.lower_bound(tmp_pattern);
	if (lowb != sim_pattern_map.end() && !(sim_pattern_map.key_comp()(tmp_pattern, lowb->first)))
		{
		// Pattern is already in sim_pattern_map, so just modify its count
		lowb->second += count;
		}
	else
		{
		// tmp_pattern has not yet been stored in sim_pattern_map
		sim_pattern_map.insert(lowb, pattern_map_t::value_type(tmp_pattern, count));
		}

	// Insert the supplied sitelist into pattern_to_sites_map
	pattern_to_sites_map_t::iterator lowbb = pattern_to_sites_map.lower_bound(tmp_pattern);
	if (lowbb != pattern_to_sites_map.end() && !(pattern_to_sites_map.key_comp()(tmp_pattern, lowbb->first)))
		{
		// Pattern is already in pattern_to_sites_map, so just modify the existing list
		uint_vect_t & existing_sitelist = lowbb->second;//UINT_LIST
        unsigned existing_nsites = (unsigned)existing_sitelist.size();
        unsigned additional_sites = (unsigned)sitelist.size();
        existing_sitelist.resize(existing_nsites + additional_sites);
		std::copy(sitelist.begin(), sitelist.end(), std::back_inserter(existing_sitelist));
        }
	else
		{
		// tmp_pattern has not yet been stored in pattern_to_sites_map
		pattern_to_sites_map.insert(lowbb, pattern_to_sites_map_t::value_type(tmp_pattern, sitelist));
		}

    // Insert tmp_pattern into sim_pattern_vect at every position listed in sitelist vector (usually of length one)
    for (uint_vect_t::const_iterator i = sitelist.begin(); i != sitelist.end(); ++i)
        {
        unsigned pos = (unsigned)(*i);
        PHYCAS_ASSERT(sim_pattern_vect.size() > pos);
        sim_pattern_vect[pos] = tmp_pattern;
        total_count += count;
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `sim_pattern_map' to the patterns already in `other'. The value of `mult' is used to
|	modify the counts before they are added to `other'; that is, the count of each pattern added to `other' is the
|	original count multiplied by `mult'. Normally, `mult' would be specified to be 1.0, but in some cases it is
|	necessary to build up SimData objects that represent averages of other SimData objects, and it is these situations
|	where `mult' comes in handy. Assumes that `mult' is positive and non-zero. Assumes that `pattern_length' for this
|   SimData object is identical to the `pattern_length' of `other'.
*/
void SimData::addDataTo(
  SimData & other, 			/**< is the SimData object that will receive the data currently contained in this SimData object */
  pattern_count_t mult)		/**< is the factor multiplied by each pattern's count before pattern is stored in `other' */
	{
	PHYCAS_ASSERT(mult > 0.0);

	// If this object has no patterns, return immediately
	if (total_count == 0.0)
		return;

	// If other is empty, then it most likely needs to be initialized
	if (other.getTotalCount() == 0)
		{
		other.resetPatternLength(pattern_length);
		}
	PHYCAS_ASSERT(pattern_length == other.getPatternLength());

    pattern_map_t::iterator pit = sim_pattern_map.begin();
	for (; pit != sim_pattern_map.end(); ++pit)
		{
		pattern_count_t count = pit->second;

		// getCurrPattern returns a workspace for building up a pattern to be added to other
		int8_vect_t & other_pattern = other.getCurrPattern();

		// Copy pattern represented by *it to other's workspace
		std::copy(pit->first.begin(), pit->first.end(), other_pattern.begin());

		// Add the pattern in other's workspace to other's pattern map
		pattern_count_t mult_count = mult*count;

		other.insertPatternOnly(mult_count);
		}
    //#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `total_count' and `pattern_length' both to zero.
*/
SimData::SimData()
  : pattern_length(0), total_count(0.0)
    , _nchar(0)
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Resets this object to its just-constructed state.
*/
void SimData::clear()
    {
    tmp_pattern.clear();
    sim_pattern_map.clear();
    total_count = 0.0;
    pattern_length = 0;
    sim_pattern_vect.clear();
    pattern_to_sites_map.clear();
    sim_pattern_vect.resize(_nchar);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Resizes sim_pattern_vect to the supplied length.
*/
void SimData::resizePatternVect(
  unsigned sz)  /**< new size of sim_pattern_vect vector */
    {
    sim_pattern_vect.resize(sz);
    _nchar = sz;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the data member `sim_pattern_map'.
*/
const pattern_map_t & SimData::getSimPatternMap() const
    {
    return sim_pattern_map;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the data member `pattern_to_sites_map'.
*/
pattern_to_sites_map_t & SimData::getPatternToSitesMap()
    {
    return pattern_to_sites_map;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	If `ntaxa' equals `pattern_length', this function returns immediately without taking any action. If `ntaxa' differs
|	from `pattern_length', sets `pattern_length' to `ntaxa' and clears `sim_pattern_map' because it is invalidated when
|	`pattern_length' changes. Also refreshes the `tmp_pattern' data member to conform to the new value of
|	`pattern_length'. Assumes `ntaxa' is non-zero. Use clear() to set `ntaxa' to 0 and return the object to its
|	just-constructed state.
*/
void SimData::resetPatternLength(
  unsigned ntaxa)	/**< is the number of taxa (same as the number of elements in a pattern vector) */
    {
    PHYCAS_ASSERT(ntaxa > 0);
    if (ntaxa == pattern_length)
        return;

    clear();
    pattern_length = ntaxa;

    // Create a int8_vect_t vector with ntaxa elements all of which are -1
    // and swap into tmp_pattern so that tmp_pattern now has the correct size
    pattern_t v(ntaxa, missing_state);	// GCC 3.2.3 (Red Hat Linux 3.2.3-20) requires this split
    tmp_pattern.swap(v);					// into two lines
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets every element in `tmp_pattern' to the value of `missing_state'.
*/
void SimData::wipePattern()
    {
    tmp_pattern.assign((pattern_t::size_type)pattern_length, SimData::missing_state);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets all counts in `sim_pattern_map' to zero.
*/
void SimData::zeroCounts()
    {
    for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
        {
        it->second = 0.0;
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Saves all counts as a tab-delimited row appended to the file whose name is 'fn'. Note: no header is output
|	identifying which pattern goes with each count, so this will not be useful unless only the unidentified counts are
|	of interest.
*/
void SimData::appendCountsToFile(
  std::string fn,   /**< is the name of the file to which counts will be appended */
  bool binary)      /**< if true, any non-zero counts will be output as 1 */
    {
    std::ofstream outf(fn.c_str(), std::ios::out | std::ios::app);
    pattern_map_t::iterator it = sim_pattern_map.begin();
    unsigned count = (unsigned)(it->second);
    if (binary && count > 0)
        count = 1;
    outf << count;
    ++it;
    for (; it != sim_pattern_map.end(); ++it)
    {
        count = (unsigned)(it->second);
        if (binary && count > 0)
            count = 1;
        outf << '\t' << count;
    }
    outf << std::endl;
    outf.close();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Saves all counts as a tab-delimited row appended to the file whose name is 'fn'. Note: no header is output
|	identifying which pattern goes with each count, so this will not be useful unless only the unidentified counts are
|	of interest.
*/
void SimData::debugAppendCountsToFile(
  std::string row_name,     /**< is the name of the row, which is output before any counts */
  std::string fn)           /**< is the name of the file to which counts will be appended */
    {
    std::ofstream outf(fn.c_str(), std::ios::out | std::ios::app);
    outf << row_name;

#if 0

    for (unsigned i = 0; i < 4; ++i)
        {
        for (unsigned j = 0; j < 4; ++j)
            {
            for (unsigned k = 0; k < 4; ++k)
                {
                for (unsigned m = 0; m < 4; ++m)
                    {
                    int8_vect_t v;
                    v.push_back(i);
                    v.push_back(j);
                    v.push_back(k);
                    v.push_back(m);
                    pattern_map_t::iterator it = sim_pattern_map.find(v);
                    PatternCountType count = 0.0;
                    if (it != sim_pattern_map.end())
                        {
                        count = (PatternCountType)(it->second);
                        }
                    outf << '\t' << count;
                    }
                }
            }
        }

#else

    pattern_map_t::iterator it = sim_pattern_map.begin();
    pattern_count_t count = (pattern_count_t)(it->second);
    outf << count;
    ++it;
    for (; it != sim_pattern_map.end(); ++it)
        {
        count = (pattern_count_t)(it->second);
        outf << '\t' << count;
        }

#endif

    outf << std::endl;
    outf.close();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of patterns, each represented as a string. The patterns are in the same order as the counts
|	returned by the appendCountsToFile function, so this function can be used to identify the patterns belonging to the
|	pattern counts returned by that function. This is a slow function, so don't use it in situations where speed is
|	critical.
*/
std::vector<std::string> SimData::getPatterns(
  std::vector<std::string> symbols) /**< is a vector of state symbols */
    {
    std::vector<std::string> v;
    for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
        {
        // Create a string out of the pattern
        std::string s;
        for (pattern_t::const_iterator i = it->first.begin(); i != it->first.end(); ++i)
            {
            s += symbols.at(*i);
            }
        v.push_back(s);
        }
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Saves all counts as a string that can be used within Maple to create a 3D surface plot. To make a surface plot of
|	the following 2x3 matrix in Maple 9,
|>
|	1  2  3
|	4  6  8
|>
|	place the following commands in a file named maple_commands:
|>
|	with(linalg);
|	with(plots);
|	points := [[[0,0,1],[0,1,2],[0,2,3]],[[1,0,4],[1,1,6],[1,2,8]]];
|	surfdata(points, axes=framed, labels=["rows", "cols", "height"]);
|>
|	then read this file by typing 'read maple_commands;' in the Maple interpreter. In our case, rows will be separate
|	posterior predictive data sets and columns will be patterns (256 columns for 4 taxa and DNA data). The height
|	will be the number of counts. This function returns a string representing one row of the points matrix. For
|	example, if this data set were going to be the 2nd row (index 1) in the above example, the string returned would be
|>
|	[[1,0,4],[1,1,5],[1,2,6]]
|>
|
*/
std::string SimData::createMapleTuples(
  unsigned row,         /**< is the row */
  unsigned cutoff)      /**< is the cutoff */
    {
    bool use_cutoff = true;
    if (cutoff == 0)
        use_cutoff = false;

    std::string s = "[";
    unsigned col = 0;
    pattern_map_t::iterator it = sim_pattern_map.begin();
    unsigned count = (unsigned)(it->second);
    if (use_cutoff && count > cutoff)
        count = 0;
    double log_count = count > 0 ? std::log((double)count) : 0.0;
    s += str(boost::format("[%d,%d,%f]") % row % col % log_count);
    ++it;
    ++col;
    for (; it != sim_pattern_map.end(); ++it, ++col)
        {
        count = (unsigned)(it->second);
        if (use_cutoff && count > cutoff)
            count = 0;
        log_count = count > 0 ? std::log((double)count) : 0.0;
        s += str(boost::format(",[%d,%d,%f]") % row % col % log_count);
        }
    s += "]";
    return s;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `tmp_pattern'[`pos'] to the supplied `state'.
*/
void SimData::setState(
  unsigned pos, 	/**< is the position in `tmp_pattern' to set */
  int8_t state)		/**< is the state to assign to the element in `tmp_pattern' at position pos */
    {
    PHYCAS_ASSERT(pos < pattern_length);
    tmp_pattern[pos] = state;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the `tmp_pattern' data member.
*/
pattern_t & SimData::getCurrPattern()
    {
    return tmp_pattern;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the `pattern_length' data member.
*/
unsigned SimData::getPatternLength()
    {
    return pattern_length;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the `total_count' data member, which represents the total number of patterns added (not
|	`sim_pattern_map'.size()).
*/
pattern_count_t SimData::getTotalCount()
    {
    return total_count;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the `total_count' data member, which represents the sum of all pattern counts (not
|	`sim_pattern_map'.size()).
*/
void SimData::setTotalCount(
  pattern_count_t total)	/**< is the new total count */
    {
    total_count = total;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value `sim_pattern_map'.size(), the number of unique patterns currently stored.
*/
unsigned SimData::getNUniquePatterns()
    {
    return (unsigned)sim_pattern_map.size();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Divides the count associated with every pattern in `sim_pattern_map' by `factor'. Also divides `total_count' by
|	`factor'. Assumes `factor' is greater than zero.
*/
void SimData::divideBy(
  pattern_count_t factor)   /**< is the factor by which counts will be divided */
    {
    PHYCAS_ASSERT(factor > 0.0);

    total_count /= factor;

    for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
        {
        it->second /= factor;
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Multiplies the count associated with every pattern in `sim_pattern_map' by `factor'. Also multiplies `total_count'
|	`factor'. Assumes `factor' is greater than zero.
*/
void SimData::multBy(
  pattern_count_t factor)  /**< is the factor by which counts will be multiplied */
    {
    if (sim_pattern_map.empty())
        return;

    total_count *= factor;

    for (pattern_map_t::iterator it = sim_pattern_map.begin(); it != sim_pattern_map.end(); ++it)
        {
        it->second *= factor;
        }
    }

}	// namespace phycas
