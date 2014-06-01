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
#include "prior_distribution.hpp"
#include "topo_prior_calculator.hpp"
#include "tree_length_distribution.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   Constructor.
*/
PriorDistribution::PriorDistribution()
  : _edgelen_type(edgelen_none), _univar_value(0.0), _current_density(0.0)
    {
    _univar_dist.reset();
    _multivar_dist.reset();
    _topology_dist.reset();
    _treelen_dist.reset();
    _tree.reset();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Destructor.
*/
PriorDistribution::~PriorDistribution()
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Checks to see whether this is the first time this object has been initialized. Returns true if this is the first
|   first initialization, false otherwise.
*/
bool PriorDistribution::_neverBeforeInitialized()
    {
    if (!_univar_dist && !_multivar_dist && !_topology_dist && !_treelen_dist)
        return true;
    else
        return false;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Initializes object with supplied `name', univariate probability distribution `dist' and value `value'. Immediately
|   computes `_current_density'.
*/
void PriorDistribution::initUnivariate(
  std::string name, ProbDistShPtr dist, double value)
    {
    PHYCAS_ASSERT(_neverBeforeInitialized());
    _name = name;
    _univar_dist = dist;
    _univar_value = value;
    _recalculateLogDensity();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Initializes object with supplied `name', multivariate probability distribution `dist' and vector `value'.
|   Immediately computes `_current_density'.
*/
void PriorDistribution::initMultivariate(
  std::string name,
  MultivarProbDistShPtr dist,
  const double_vect_t & value)
    {
    PHYCAS_ASSERT(_neverBeforeInitialized());
    _name = name;
    _multivar_dist = dist;
    _multivar_value.assign(value.begin(), value.end());
    _recalculateLogDensity();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Initializes object with supplied `name', univariate probability distribution `dist' and tree `tree'. Object will be
|   responsible for managing the joint prior over all external edge lengths. Immediately computes `_current_density'.
*/
void PriorDistribution::initExternalEdgeLenDistribution(
  std::string name,
  ProbDistShPtr dist,
  TreeShPtr tree)
    {
    PHYCAS_ASSERT(_neverBeforeInitialized());
    _univar_dist = dist;
    _name = name;
    _tree = tree;
    _edgelen_type = edgelen_external;
    _recalculateLogDensity();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Initializes object with supplied `name', univariate probability distribution `dist' and tree `tree'. Object will be
|   responsible for managing the joint prior over all internal edge lengths. Immediately computes `_current_density'.
*/
void PriorDistribution::initInternalEdgeLenDistribution(
  std::string name,
  ProbDistShPtr dist,
  TreeShPtr tree)
    {
    PHYCAS_ASSERT(_neverBeforeInitialized());
    _univar_dist = dist;
    _name = name;
    _tree = tree;
    _edgelen_type = edgelen_internal;
    _recalculateLogDensity();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Initializes object with supplied `name', tree length distribution `dist' and tree `tree'. Object will be
|   responsible for managing the joint prior over all edge lengths. Immediately computes `_current_density'.
*/
void PriorDistribution::initTreeLenDistribution(
  std::string name,
  TreeLengthDistributionShPtr dist,
  TreeShPtr tree)
    {
    PHYCAS_ASSERT(_neverBeforeInitialized());
    _treelen_dist = dist;
    _name = name;
    _tree = tree;
    _recalculateLogDensity();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Initializes object with supplied `name', topology probability calculator `topo_prior_calculator' and tree `tree'.
|   Immediately computes `_current_density'.
*/
void PriorDistribution::initTopologyDistribution(
  std::string name,
  TopoProbCalculatorShPtr topo_prior_calculator,
  TreeShPtr tree)
    {
    PHYCAS_ASSERT(_neverBeforeInitialized());
    _topology_dist = topo_prior_calculator;
    _name = name;
    _tree = tree;
    _recalculateLogDensity();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Recalculates and returns `_current_density' based on the current stored value.
*/
double PriorDistribution::_recalculateLogDensity()
    {
    PHYCAS_ASSERT(!_neverBeforeInitialized());
    if (_univar_dist)
        {
        if (_tree)
            {
            // This object manages the joint prior on edge lengths (either all of them or just internal or just external)
            PHYCAS_ASSERT(_edgelen_type != edgelen_none);
            if (_edgelen_type == edgelen_all)
                {
                _current_density = 0.0;
                _univar_value = 0.0;
                preorder_iterator it = _tree->begin();  // always skip root node because it has no edge
                for (++it; it != _tree->end(); ++it)
                    {
                    double v = it->GetEdgeLen();
                    _univar_value += v;
                    _current_density += _univar_dist->GetLnPDF(v);
                    }
                }
            else if (_edgelen_type == edgelen_internal)
                {
                _current_density = 0.0;
                _univar_value = 0.0;
                preorder_iterator it = _tree->begin();  // always skip root node because it has no edge
                for (++it; it != _tree->end(); ++it)
                    {
                    if (it->IsInternal() && !it->IsSubroot())
                        {
                        double v = it->GetEdgeLen();
                        _univar_value += v;
                        _current_density += _univar_dist->GetLnPDF(v);
                        }
                    }
                }
            else
                {
                _current_density = 0.0;
                _univar_value = 0.0;
                preorder_iterator it = _tree->begin();  // always skip root node because it has no edge
                for (++it; it != _tree->end(); ++it)
                    {
                    if (it->IsTip() || it->IsSubroot())
                        {
                        double v = it->GetEdgeLen();
                        _univar_value += v;
                        _current_density += _univar_dist->GetLnPDF(v);
                        }
                    }
                }
            }
        else if (_multivar_value.size() > 0)
            {
            // This object manages the joint prior on some collection of univariate parameters (currently unused)
            _current_density = 0.0;
            for (double_vect_t::const_iterator it = _multivar_value.begin(); it != _multivar_value.end(); ++it)
                _current_density += _univar_dist->GetLnPDF(*it);
            }
        else
            {
            // This object manages the prior on a univariate parameter (e.g. kappa, gamma shape, pinvar, etc.)
            _current_density = _univar_dist->GetLnPDF(_univar_value);
            }
        }
    else if (_multivar_dist)
        {
        // This object manages the prior on a multivariate parameter (e.g. state frequencies, GTR exchangeabilities, subset relative rates, etc.)
        _current_density = _multivar_dist->GetLnPDF(_multivar_value);
        }
    else if (_treelen_dist)
        {
        // This object manages the joint prior governing both tree length and individual edges
        PHYCAS_ASSERT(_tree);
        _current_density = _treelen_dist->GetLnPDF(_tree);
        }
    else
        {
        // This object manages the prior on the tree topology
        PHYCAS_ASSERT(_tree);
        _current_density = _topology_dist->GetLnTopoProb(_tree);
        }
    return _current_density;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Shows the current value of the parameter on the supplied output stream (e.g. std::cerr).
*/
void PriorDistribution::debugShowCurrValue(
  std::ostream & out) const  /**< is the output stream on which to display the value */
    {
    const char * spacer = " ";
    if (_univar_dist)
        {
        if (_tree)
            {
            // This object manages the joint prior on edge lengths (either all of them or just internal or just external)
            PHYCAS_ASSERT(_edgelen_type != edgelen_none);
            if (_edgelen_type == edgelen_all)
                {
                preorder_iterator it = _tree->begin();  // always skip root node because it has no edge
                for (++it; it != _tree->end(); ++it)
                    {
                    double v = it->GetEdgeLen();
                    out << v << spacer;
                    }
                }
            else if (_edgelen_type == edgelen_internal)
                {
                preorder_iterator it = _tree->begin();  // always skip root node because it has no edge
                for (++it; it != _tree->end(); ++it)
                    {
                    if (it->IsInternal() && !it->IsSubroot())
                        {
                        double v = it->GetEdgeLen();
                        out << v << spacer;
                        }
                    }
                }
            else
                {
                preorder_iterator it = _tree->begin();  // always skip root node because it has no edge
                for (++it; it != _tree->end(); ++it)
                    {
                    if (it->IsTip() || it->IsSubroot())
                        {
                        double v = it->GetEdgeLen();
                        out << v << spacer;
                        }
                    }
                }
            }
        else if (_multivar_value.size() > 0)
            {
            // This object manages the joint prior on some collection of univariate parameters (currently unused)
            for (double_vect_t::const_iterator it = _multivar_value.begin(); it != _multivar_value.end(); ++it)
                {
                double v = *it;
                out << v << spacer;
                }
            }
        else
            {
            // This object manages the prior on a univariate parameter (e.g. kappa, gamma shape, pinvar, etc.)
            out << _univar_value;
            }
        }
    else if (_multivar_dist)
        {
        // This object manages the prior on a multivariate parameter (e.g. state frequencies, GTR exchangeabilities, subset relative rates, etc.)
        std::copy(_multivar_value.begin(), _multivar_value.end(), std::ostream_iterator<double>(out, spacer));
        }
    else if (_treelen_dist)
        {
        // This object manages the joint prior governing both tree length and individual edges
        PHYCAS_ASSERT(_tree);
        out << _tree->MakeNewick();
        }
    else
        {
        // This object manages the prior on the tree topology
        PHYCAS_ASSERT(_tree);
        out << _tree->MakeNewick();
        }
    }

}	// namespace phycas
