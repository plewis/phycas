/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2012 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	define NO_IMPORT_ARRAY
#endif

#include <numeric>
#include <functional>
#include <cmath>
#include <boost/lambda/lambda.hpp>

#include "tree_length_distribution.hpp"
#include "dirichlet_distribution.hpp"

#include <fstream>  //temporary!
#include <iterator>  //temporary!

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	By default, let the tree length distribution be exponential with mean 10 and the edge length distribution
|   conditional on tree length be flat.
*/
TreeLengthDistribution::TreeLengthDistribution()
  : _alphaT(1.0), _betaT(0.1), _alpha(1.0), _c(1.0), _lot(&_myLot), _num_internal_edges(0), _num_external_edges(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a TreeLengthDistribution based on the parameter values supplied.
*/
TreeLengthDistribution::TreeLengthDistribution(
  double alphaT,     /**< is the shape parameter of the gamma distribution of tree length (mean = shape/scale) */
  double betaT,      /**< is the scale parameter of the gamma distribution of tree length (mean = shape/scale) */
  double alpha,      /**< is the Dirichlet parameter governing external edge length */
  double c)          /**< is the Dirichlet parameter governing internal edge length */
  :  _alphaT(alphaT), _betaT(betaT), _alpha(alpha), _c(c), _lot(&_myLot), _num_internal_edges(0), _num_external_edges(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameters to the values contained in their counterparts in `other'.
*/
TreeLengthDistribution::TreeLengthDistribution(
  const TreeLengthDistribution & other)	/* the tree length distribution to clone */
  : _alphaT(other._alphaT), _betaT(other._betaT), _alpha(other._alpha), _c(other._c), _lot(other._lot), _num_internal_edges(0), _num_external_edges(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the string "TreeLengthDistribution".
*/
std::string TreeLengthDistribution::GetDistributionName() const
	{
	return "TreeLengthDistribution";
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a string similar to "TreeLengthDistribution(1.0,0.1,1.0,1.0)".
*/
std::string TreeLengthDistribution::GetDistributionDescription() const
	{
	return boost::str(boost::format("TreeLengthDistribution(%#.5f,%#.5f,%#.5f,%#.5f)") % _alphaT % _betaT % _alpha % _c);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates two probability distributions that can be used to sample from this tree length distribution. The first is
|   a Gamma distribution stored in the data member `_tldist' and the second is a Dirichlet distribution stored in the
|   data member `_eldist'. The edge length distribution depends on the number of internal and external edges, so this
|   function must be called whenever either of these numbers changes. The number of internal and external edges last
|   used is stored in `_num_internal_edges' and `_num_external_edges', respectively, to avoid unnecessary
|   construction/destruction.
*/
void TreeLengthDistribution::SetupSamplingDistributions(
  unsigned num_external, /**< is the number of external edges */
  unsigned num_internal) /**< is the number of internal edges */
	{
    unsigned num_total = num_external + num_internal;
    PHYCAS_ASSERT(num_total > 0);
    if (num_external == _num_external_edges && num_internal == _num_internal_edges && _tldist && _eldist)
        {
        return;
        }

    //std::cerr << "Creating _tldist and _eldist for " << num_internal << " internal and " << num_external << " external edges" << std::endl;

    double_vect_t dirichlet_params(num_total, _alpha);

    //std::cerr << "\n\n***** dirichlet_params before *****" << std::endl;
    //std::copy(dirichlet_params.begin(), dirichlet_params.end(), std::ostream_iterator<double>(std::cerr, "|"));
    //std::cerr << "\n\n" << std::endl;

    std::transform(dirichlet_params.begin()+num_external, dirichlet_params.end(), dirichlet_params.begin()+num_external, boost::lambda::_1*_c);

    //std::cerr << "\n\n***** dirichlet_params after *****" << std::endl;
    //std::copy(dirichlet_params.begin(), dirichlet_params.end(), std::ostream_iterator<double>(std::cerr, "|"));
    //std::cerr << "\n\n" << std::endl;

    _tldist.reset(new GammaDistribution(_alphaT, 1.0/_betaT));   // Rannala-Zhu-Yang paper defines Gamma mean = shape/scale, which contrasts with other parts of phycas where Gamma mean = shape*scale
    _tldist->SetLot(_lot);

    _eldist .reset(new DirichletDistribution(dirichlet_params));
    _eldist->SetLot(_lot);

    _num_internal_edges = num_internal;
    _num_external_edges = num_external;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Samples and returns a point from this distribution. In the returned point, the `num_external' external edge lengths
|   will be first, followed by the `num_internal' internal edge lengths.
*/
double_vect_t TreeLengthDistribution::Sample(
  unsigned num_external, /**< is the number of external edges */
  unsigned num_internal) /**< is the number of internal edges */
	{
    SetupSamplingDistributions(num_external, num_internal);

    double tree_length = _tldist->Sample();
    double_vect_t edge_lengths = _eldist->Sample();
    std::transform(edge_lengths.begin(), edge_lengths.end(), edge_lengths.begin(), boost::lambda::_1*tree_length);

	return edge_lengths;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the probability density function evaluated for the edge lengths in the supplied
|   tree `t'. Assumes supplied tree `t' is unrooted. See Rannala, Zhu and Yang (2012) for details.
*/
double TreeLengthDistribution::GetLnPDF(
  TreeShPtr t) const 	/**< is the tree for which the density is to be evaulated */
	{
    PHYCAS_ASSERT(!t->IsRooted());

    double tree_length = 0.0;
    double sum_internal = 0.0;
    double sum_external = 0.0;
    unsigned num_total = 0;
    unsigned num_internal = 0;
    unsigned num_external = 0;
    for (preorder_iterator it = t->begin(); it != t->end(); ++it)
        {
        if (!it->IsAnyRoot())
            {
            double edge_length = it->GetEdgeLen();
            tree_length += edge_length;
            num_total++;
            if (it->IsTip() || it->IsSubroot())
                {
                num_external++;
                sum_external += (_alpha - 1.0)*log(edge_length);
                }
            else
                {
                num_internal++;
                sum_internal += (_alpha*_c - 1.0)*log(edge_length);
                }
            }
        }

    PHYCAS_ASSERT(num_total == num_internal + num_external);

    // logterm1 is first line of equation 36 in Rannala, Zhu, and Yang (2012)
    double logterm1 = _alphaT*log(_betaT) - _cdf.LnGamma(_alphaT) - _betaT*tree_length + (_alphaT - 1.0)*log(tree_length);

    // logterm2 is second line of equation 36 in Rannala, Zhu, and Yang (2012)
    double beta_term = _cdf.LnGamma(_alpha)*num_external + _cdf.LnGamma(_alpha*_c)*num_internal - _cdf.LnGamma(_alpha*num_external + _alpha*_c*num_internal);
    double logterm2 = sum_internal + sum_external - beta_term;

    // logterm3 is third line of equation 36 in Rannala, Zhu, and Yang (2012)
    double logterm3 = (1.0 - _alpha*num_external - _alpha*_c*num_internal)*log(tree_length);

    double lnPDF = logterm1 + logterm2 + logterm3;
    return lnPDF;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the probability density function evaluated for the edge lengths in the supplied
|   tree `t'. Assumes supplied tree `t' is unrooted. See Rannala, Zhu and Yang (2012) for details.
*/
double TreeLengthDistribution::GetRelativeLnPDF(
  TreeShPtr t) const  	/**< is the tree for which the density is to be evaulated */
	{
	return GetLnPDF(t);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object, calls the new object's SetLot member function (passing the
|	supplied Lot object `other'), and returns a pointer to it. The caller is expected to manage the new object.
*/
TreeLengthDistribution * TreeLengthDistribution::cloneAndSetLot(Lot * other) const
	{
    TreeLengthDistribution * clone = new TreeLengthDistribution(*this);
	clone->SetLot(other);
	return clone;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object.
*/
TreeLengthDistribution * TreeLengthDistribution::Clone() const
	{
    return new TreeLengthDistribution(*this);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the TreeLengthDistribution::Sample member function. The original
|	random number generator	(data member `myLot') can be replaced by calling the ProbabilityDistribution::ResetLot
|	function. Note that this object does not take ownership of the Lot object whose pointer is specified as `other'.
|	It is assumed that `other' is non-NULL.
*/
void TreeLengthDistribution::SetLot(
	Lot * other) /**< is a pointer to the random number generator object to be used subsequently by Sample */
	{
	if (other == NULL)
		throw XProbDist("attempt made to install a non-existent pseudorandom number generator");
	_lot = other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number seed for the `myLot' data member. Note that if TreeLengthDistribution::SetLot has been
|	called, calling TreeLengthDistribution::SetSeed is pointless because you will not be setting the seed for the
|	correct random number generator!
*/
void TreeLengthDistribution::SetSeed(
  unsigned rnseed)	/**< is the new seed value */
	{
	_myLot.SetSeed(rnseed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes the data member `_lot' (which is used as the random number generator by the member function
|	TreeLengthDistribution::Sample) point to the local data member _myLot. This function only needs to be called if
|	TreeLengthDistribution::SetLot has been called previously to replace the random number generator used by
|	TreeLengthDistribution::Sample.
*/
void TreeLengthDistribution::ResetLot()
	{
	_lot = &_myLot;
	}

} // namespace phycas


