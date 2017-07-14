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
#include <sstream>
#include <iterator>
#include "basic_tree.hpp"
#include "topo_prior_calculator.hpp"
#include "joint_prior_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   Constructor.
*/
JointPriorManager::JointPriorManager()
  : _log_joint_prior(0.0), _skip_next_debug_check(false)
    {
    _internal_hyper_key     = "internal_hyper";
    _external_hyper_key     = "external_hyper";
    _universal_hyper_key    = "edgelen_hyper";
    _internal_edgelens_key  = "internal_edgelen";
    _external_edgelens_key  = "external_edgelen";
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Destructor.
*/
JointPriorManager::~JointPriorManager()
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns natural logarithm of the joint prior density.
*/
double JointPriorManager::getLogJointPrior() const
    {
    return _log_joint_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns natural logarithm of the topology prior probability.
*/
double JointPriorManager::getLogTopologyPrior() const
    {
    //prior_distr_map_t = std::map< std::string, PriorDistrShPtr >
    prior_distr_map_t::const_iterator it = _distr_map.find("tree_topology");
    PriorDistrShPtr p = it->second;
    return p->_current_density;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Recalculates `_log_joint_prior'. This function should be called periodically (e.g. after each cycle) in order to
|   prevent creep due to accumulated round-off errors.
*/
void JointPriorManager::recalcLogJointPrior()
    {
    unsigned n = (unsigned)_distr_vect.size();
    if (n > 0)
        {
        _log_joint_prior = 0.0;
        for(prior_distr_vect_t::const_iterator it = _distr_vect.begin(); it != _distr_vect.end(); ++it)
            {
            PriorDistrShPtr p = *it;
            double my_density = p->_recalculateLogDensity();
            _log_joint_prior += my_density;
            }
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Sets `_skip_next_debug_check' to true, which means next call to debugCheckLogJointPrior will be a no-op. This is
|   useful if two different kinds of priors are modified simultaneously, and updating the first would cause the debug
|   check to fail because the second prior would not yet have been updated.
*/
void JointPriorManager::skipNextDebugCheck()
    {
    _skip_next_debug_check = true;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns true if complete recalculation of `_log_joint_prior' results in the current value of `_log_joint_prior' to
|   within a tolerance of 1.e-10; otherwise, returns false. Should be used only in debug mode.
*/
bool JointPriorManager::debugCheckLogJointPrior(
  std::string called_from) const
    {
    double tol = 1.e-8;

    if (_skip_next_debug_check)
        {
        _skip_next_debug_check = false;
        return true;
        }

    _debug_tmp.clear();
    //std::cerr << "|~o~o~o~o~o~o~> debugCheckLogJointPrior <-- " << called_from << std::endl;

    double _debug_log_joint_prior = 0.0;
    unsigned n = (unsigned)_distr_vect.size();
    if (n == 0)
        {
        return true;
        }
    else
        {
        for(prior_distr_vect_t::const_iterator it = _distr_vect.begin(); it != _distr_vect.end(); ++it)
            {
            PriorDistrShPtr p = *it;
            //std::cerr << "--> " << p->_name;
            double curr_density = p->_current_density;
            double my_density = p->_recalculateLogDensity();
            if (fabs(my_density - curr_density) > tol)
                _debug_tmp.push_back(p->_name);
            //std::cerr << ": " << my_density << " (curr_density = " << curr_density << ", ";
            //p->debugShowCurrValue(std::cerr);
            //std::cerr << ")" << std::endl;
            _debug_log_joint_prior += my_density;
            }
        }
    double absolute_difference = fabs(_debug_log_joint_prior - _log_joint_prior);
    bool good = absolute_difference < tol ? true : false;
    if (!good)
        {
        std::cerr << ">>> bad <<<" << std::endl;
        std::cerr << "_log_joint_prior    = " << _log_joint_prior       << std::endl;
        std::cerr << "recalculated        = " << _debug_log_joint_prior << std::endl;
        std::cerr << "absolute difference = " << absolute_difference    << std::endl;
        std::cerr << "offending parameters: ";
        std::copy(_debug_tmp.begin(), _debug_tmp.end(), std::ostream_iterator<std::string>(std::cerr, " "));
        std::cerr << std::endl;
        std::cerr << "Breakdown:" << std::endl;
        std::cerr << debugPriorBreakdown(NULL, NULL, NULL) << std::endl;
        }
    //else
    //    {
    //    std::cerr << ">>> good! <<<" << std::endl;
    //    std::cerr << "_log_joint_prior    = " << _log_joint_prior       << std::endl;
    //    std::cerr << "recalculated        = " << _debug_log_joint_prior << std::endl;
    //    }

    //std::cerr << boost::str(boost::format("--> %.20f <--") % absolute_difference) << std::endl;

    return good;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a string containing a summary of the joint prior, showing the contributions of every parameter.
*/
std::string JointPriorManager::debugPriorBreakdown(
  const char * msg,
  const char * before_text,
  const char * after_text) const
    {
    std::string s = "";
    if (before_text)
        s += before_text;
    s += "Joint prior breakdown";
    if (msg)
        s += " (" + std::string(msg) + ")";
    s += ":";
    unsigned n = (unsigned)_distr_vect.size();
    if (n == 0)
        {
        s += "\n  No parameters are currently being managed";
        }
    else
        {
        double recalculated = 0.0;
        s += boost::str(boost::format("\n%20s %15s  %s") % "param" % "density" % "value");
        s += boost::str(boost::format("\n%20s %15s  %15s") % "--------------------" % "---------------" % "---------------");
        for(prior_distr_vect_t::const_iterator it = _distr_vect.begin(); it != _distr_vect.end(); ++it)
            {
            PriorDistrShPtr p = *it;
            recalculated += p->_current_density;
            if (p->_univar_dist)
                {
                if (p->_multivar_value.size() > 0)
                    {
                    double_vect_t const & v = p->_multivar_value;
                    std::stringstream ss;
                    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(ss, " "));
                    s += boost::str(boost::format("\n%20s %15.5f  %s") % p->_name % p->_current_density % ss.str());
                    }
                else
                    {
                    s += boost::str(boost::format("\n%20s %15.5f  %.5f") % p->_name % p->_current_density % p->_univar_value);
                    }
                }
            else if (p->_multivar_dist)
                {
                double_vect_t const & v = p->_multivar_value;
                std::stringstream ss;
                std::copy(v.begin(), v.end(), std::ostream_iterator<double>(ss, " "));
                s += boost::str(boost::format("\n%20s %15.5f  %s") % p->_name % p->_current_density % ss.str());
                }
            else
                {
                PHYCAS_ASSERT(p->_topology_dist || p->_treelen_dist);
                PHYCAS_ASSERT(p->_tree);
                s += boost::str(boost::format("\n%20s %15.5f  %s") % p->_name % p->_current_density % "---");    //p->_tree->MakeNewick()
                }
            }
        s += boost::str(boost::format("\n%20s %15s  %15s") % "--------------------" % "---------------" % "---------------");
        s += boost::str(boost::format("\n%20s %15.5f") % "Total" % recalculated);
        s += boost::str(boost::format("\n%20s %15.5f") % "_log_joint_prior" % _log_joint_prior);
        }
    if (after_text)
        s += after_text;
    return s;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a univariate probability distribution to `_distr_map' and `_distr_vect'. The key in the map will be the
|   `keystr' provided (which is assumed to be unique), while the value associated with the key is a PriorDistribution
|   object containing the supplied `dist'.
*/
void JointPriorManager::addUnivariateDistribution(
  std::string keystr,     /**< is the string used to identify this distribution */
  ProbDistShPtr dist,     /**< is the univariate distribution to manage */
  double initial_value)   /**< is a pointer to the initial value, pass NULL to draw from `dist' */
    {
    //std::cerr << "JointPriorManager::addUnivariateDistribution: " << keystr << std::endl;
    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());
    p->initUnivariate(keystr, dist, initial_value);
    _distr_map[keystr] = p;
    _distr_vect.push_back(p);
    _log_joint_prior += p->_recalculateLogDensity();
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("addUnivariateDistribution (%s)") % keystr)));
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a multivariate probability distribution to `_distr_map' and `_distr_vect'. The key in the map will be the
|   `keystr' provided (which is assumed to be unique), while the value associated with the key is a PriorDistribution
|   object containing the supplied `dist'.
*/
void JointPriorManager::addMultivariateDistribution(
  std::string keystr,                      /**< is the string used to identify this distribution */
  MultivarProbDistShPtr dist,              /**< is the multivariate distribution to manage */
  double_vect_t const & initial_value)     /**< is a pointer to the initial value vector, pass NULL to draw from `dist' */
    {
    //std::cerr << "JointPriorManager::addMultivariateDistribution: " << keystr << std::endl;
    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());
    p->initMultivariate(keystr, dist, initial_value);
    _distr_map[keystr] = p;
    _distr_vect.push_back(p);
    _log_joint_prior += p->_recalculateLogDensity();
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("addMultivariateDistribution (%s)") % keystr)));
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a TreeLengthDistribution to `_distr_map' and `_distr_vect'. The key in the map will be the `keystr' provided
|   (which is assumed to be unique), while the value associated with the key is a PriorDistribution object containing the
|   supplied `dist'.
*/
void JointPriorManager::addTreeLengthDistribution(
  std::string keystr,                      /**< is the string used to identify this distribution */
  TreeLengthDistributionShPtr dist,        /**< is the TreeLengthDistribution to manage */
  TreeShPtr tree)                          /**< is the tree to use in computing initial probability */
    {
    //std::cerr << "JointPriorManager::addTreeLengthDistribution: " << keystr << std::endl;
    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());
    p->initTreeLenDistribution(keystr, dist, tree);
    _distr_map[keystr] = p;
    _distr_vect.push_back(p);
    _log_joint_prior += p->_recalculateLogDensity();
    PHYCAS_ASSERT(debugCheckLogJointPrior("addTreeLengthDistribution"));
    }

#define IGNORING_TOPOLOGY_PRIOR false

/*----------------------------------------------------------------------------------------------------------------------
|   Returns `_topology_dist' data member for the topology prior distribution object. The return value is a shared
|   pointer to either a TopoProbCalculator or a PolytomyTopoPriorCalculator.
*/
TopoProbCalculatorShPtr JointPriorManager::getTopoProbCalculator()
    {
    PriorDistrShPtr p = _distr_map["tree_topology"];
    PHYCAS_ASSERT(p);
    return p->_topology_dist;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a PriorDistribution node to both `_distr_map' and `_distr_vect' that has the ability to compute the probability
|   of a tree topology. If the supplied TopoProbCalculatorShPtr `topo_prior_calculator' is NULL, uses the default
|   discrete uniform distribution over fully-resolved tree topologies to compute the probability of a tree topology.
*/
void JointPriorManager::addTopologyDistribution(
  std::string keystr,                               /**< is the string used to identify this distribution */
  TopoProbCalculatorShPtr topo_prior_calculator,    /**< is the topology prior calculator (if NULL, uses discrete uniform distribution over binary trees) */
  TreeShPtr tree)                                   /**< is the tree to use in computing initial probability */
    {
    if (IGNORING_TOPOLOGY_PRIOR)
        {
        std::cerr << "IGNORING JointPriorManager::addTopologyDistribution: " << keystr << std::endl;
        return;
        }

    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());

    // following line calls PriorDistribution::_recalculateLogDensity
    p->initTopologyDistribution(keystr, topo_prior_calculator, tree);

    _distr_map[keystr] = p;
    _distr_vect.push_back(p);

    // following line calls PriorDistribution::_recalculateLogDensity
    _log_joint_prior += p->_recalculateLogDensity();

    // following line calls PriorDistribution::_recalculateLogDensity in debug mode
    PHYCAS_ASSERT(debugCheckLogJointPrior("addTopologyDistribution"));
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a univariate probability distribution to `_distr_map' and `_distr_vect' that applies to each of the internal
|   edge lengths. The `tree' serves as the source of starting values.
*/
void JointPriorManager::addInternalEdgelenDistribution(
  std::string keystr,     /**< is the string used to identify this distribution */
  ProbDistShPtr dist,     /**< is the univariate distribution that applies to each internal edge length */
  TreeShPtr tree)         /**< is the tree from which to obtain the starting value */
    {
    PHYCAS_ASSERT(_internal_edgelens_key == keystr);
    //std::cerr << "JointPriorManager::addInternalEdgelenDistribution: " << keystr << std::endl;
    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());
    p->initInternalEdgeLenDistribution(keystr, dist, tree);
    _distr_map[keystr] = p;
    _distr_vect.push_back(p);
    _log_joint_prior += p->_recalculateLogDensity();
    PHYCAS_ASSERT(debugCheckLogJointPrior("addInternalEdgelenDistribution"));
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a univariate probability distribution to `_distr_map' and `_distr_vect' that applies to each of the external
|   edge lengths. The `tree' serves as the source of starting values.
*/
void JointPriorManager::addExternalEdgelenDistribution(
  std::string keystr,     /**< is the string used to identify this distribution */
  ProbDistShPtr dist,     /**< is the univariate distribution that applies to each external edge length */
  TreeShPtr tree)         /**< is the tree from which to obtain the starting value */
    {
    PHYCAS_ASSERT(_external_edgelens_key == keystr);
    //std::cerr << "JointPriorManager::addExternalEdgelenDistribution: " << keystr << std::endl;
    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());
    p->initExternalEdgeLenDistribution(keystr, dist, tree);
    _distr_map[keystr] = p;
    _distr_vect.push_back(p);
    _log_joint_prior += p->_recalculateLogDensity();
    PHYCAS_ASSERT(debugCheckLogJointPrior("addExternalEdgelenDistribution"));
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Creates and returns a string showing the keys currently stored in `_distr_map', used for debugging purposes.
*/
std::string JointPriorManager::_debugShowDistrMapKeys()
    {
    std::string s = "|";
    for (prior_distr_map_t::const_iterator it = _distr_map.begin(); it != _distr_map.end(); ++it)
        {
        prior_distr_map_t::value_type p = *it;
        s += p.first;
        s += "|";
        }
    return s;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns true if supplied `name' is a key in `_distr_map', false otherwise.
*/
bool JointPriorManager::_checkDistrMapKey(
  std::string name) const
    {
    prior_distr_map_t::const_iterator lowb = _distr_map.lower_bound(name);
    if (lowb != _distr_map.end() && !(_distr_map.key_comp()(name, lowb->first)))
        {
        // name is a key
        return true;
        }
    else
        {
        // name is not a key
        //std::cerr << "\n_distr_map keys: " << _debugShowDistrMapKeys() << std::endl;
        return false;
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current univariate value with supplied value `v'. Recalculates `_current_density' for the parameter and
|   updates `_log_joint_prior'.
*/
double JointPriorManager::univariateModified(
  std::string name,
  double v)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    p->_univar_value = v;
    double prev_logd = p->_current_density;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("univariateModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Identical to JointPriorManager::univariateModified except that no call to debugCheckLogJointPrior is issued. This is
|   necessary when many univariate priors must be updated before a debug check is made.
*/
double JointPriorManager::univariateModifiedNoDebugCheck(
  std::string name,
  double v)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    p->_univar_value = v;
    double prev_logd = p->_current_density;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current multivariate value with supplied vector `v'. Recalculates `_current_density' for the parameter and
|   updates `_log_joint_prior'.
*/
double JointPriorManager::multivariateModified(
  std::string name,
  const double_vect_t & v)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    PHYCAS_ASSERT(p);
    p->_multivar_value.assign(v.begin(), v.end());
    double prev_logd = p->_current_density;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("multivariateModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Determines whether edge lengths are represented individually (fixed-tree-topology analyses) or as a group (variable
|   tree topology analyses, or fixed-tree-topology analyses using the compound dirichlet tree length prior), and calls
|   the appropriate function to recalculate `_current_density' and update `_log_joint_prior'.
*/
double JointPriorManager::allEdgeLensModified(
  TreeShPtr t)
    {
    double new_logd = 0.0;
    if (_checkDistrMapKey(_external_edgelens_key))
        {
        // Variable tree topology
        PHYCAS_ASSERT(_checkDistrMapKey(_internal_edgelens_key));   // if "external_edgelen" is a key, then "internal_edgelen" should be one too
        new_logd = externalAndInternalEdgeLensModified(_external_edgelens_key, _internal_edgelens_key, t);
        }
    else if (_checkDistrMapKey("tree_length"))
        {
        // Tree length prior being used
        new_logd = treeLengthModified("tree_length", t);
        }
    else
        {
        // Fixed tree topology: edge lengths under the control of individual edge length parameters
        preorder_iterator it = t->begin();  // always skip root node because it has no edge
        for (++it; it != t->end(); ++it)
            {
            // construct correct key
            unsigned node_number = it->GetNodeNumber();

            std::string k;
            if (it->IsInternal() && !it->IsSubroot())
                k = boost::str(boost::format("intedge_%d") % node_number);
            else
                k = boost::str(boost::format("extedge_%d") % node_number);

            double v = it->GetEdgeLen();

            PHYCAS_ASSERT(_checkDistrMapKey(k));
            new_logd = univariateModified(k, v);
            }
        }
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Notifies manager that at least one external edge length has been modified. Recalculates `_current_density' and
|   updates `_log_joint_prior'.
*/
double JointPriorManager::externalEdgeLensModified(
  std::string name,
  TreeShPtr t)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    double prev_logd = p->_current_density;
    p->_tree = t;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("externalEdgeLensModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Notifies manager that at least one internal edge length has been modified. Recalculates `_current_density' and
|   updates `_log_joint_prior'.
*/
double JointPriorManager::internalEdgeLensModified(
  std::string name,
  TreeShPtr t)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    double prev_logd = p->_current_density;
    p->_tree = t;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("internalEdgeLensModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Notifies manager that at least one external and at least one internal edge length has been modified. Recalculates
|   `_current_density' and updates `_log_joint_prior'. This function is essentially just externalEdgeLensModified
|   followed immediately by internalEdgeLensModified, but was created because debugCheckLogJointPrior
|   will fail if both external and internal edges have changed but debugCheckLogJointPrior is called inside
|   externalEdgeLensModified before internalEdgeLensModified has a chance to add its contribution to `_log_joint_prior'.
*/
double JointPriorManager::externalAndInternalEdgeLensModified(
  std::string name1,
  std::string name2,
  TreeShPtr t)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name1));
    PHYCAS_ASSERT(_checkDistrMapKey(name2));
    PriorDistrShPtr p1 = _distr_map[name1];
    PriorDistrShPtr p2 = _distr_map[name2];
    double prev_logd = p1->_current_density;
    prev_logd += p2->_current_density;
    p1->_tree = t;
    p2->_tree = t;
    double new_logd = p1->_recalculateLogDensity();
    new_logd += p2->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior("externalAndInternalEdgeLensModified"));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Notifies manager that at least one edge length has been modified. Should only be called when the tree length prior
|   is in use; otherwise use internalEdgeLensModified and externalEdgeLensModified instead. Recalculates
|   `_current_density' and updates `_log_joint_prior'.
*/
double JointPriorManager::treeLengthModified(
  std::string name,
  TreeShPtr t)
    {
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    double prev_logd = p->_current_density;
    p->_tree = t;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("treeLengthModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current internal edge length prior mean with supplied value `new_mean'. Calls internalEdgeLensModified
|   after replacing prior mean to force recalculation of the internal edge length prior.
*/
void JointPriorManager::_replaceInternalEdgelenPriorMean(
  double new_mean,  /**< is the new internal edge length prior mean */
  TreeShPtr tree)   /**< is the tree */
    {
    // Replace the mean of the prior on internal edge lengths
    if (_checkDistrMapKey(_internal_edgelens_key))
        {
        PriorDistrShPtr p = _distr_map[_internal_edgelens_key];
        p->_univar_dist->SetMeanAndVariance(new_mean, new_mean*new_mean);

        // Force recalculation of the internal edge length prior
        internalEdgeLensModified(_internal_edgelens_key, tree);
        }
    else
        {
        // Fixed tree topology, so there is a separate prior in the joint prior manager for each edge length parameter
        // Find all the entries corresponding to internal edge lengths and change their prior means
        for(prior_distr_vect_t::const_iterator it = _distr_vect.begin(); it != _distr_vect.end(); ++it)
            {
            PriorDistrShPtr p = *it;
            if (p->_name.find("intedge_") != std::string::npos)
                {
                p->_univar_dist->SetMeanAndVariance(new_mean, new_mean*new_mean);
                univariateModifiedNoDebugCheck(p->_name, p->_univar_value);
                }
            }

        // Now perform the debug prior check suppressed above
        PHYCAS_ASSERT(debugCheckLogJointPrior("_replaceInternalEdgelenPriorMean"));
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current external edge length prior mean with supplied value `new_mean'. Calls externalEdgeLensModified
|   after replacing prior mean to force recalculation of the external edge length prior.
*/
void JointPriorManager::_replaceExternalEdgelenPriorMean(
  double new_mean,  /**< is the new external edge length prior mean */
  TreeShPtr tree)   /**< is the tree */
    {
    // Replace the mean of the prior on external edge lengths
    if (_checkDistrMapKey(_external_edgelens_key))
        {
        PriorDistrShPtr p = _distr_map[_external_edgelens_key];
        p->_univar_dist->SetMeanAndVariance(new_mean, new_mean*new_mean);

        // Force recalculation of the external edge length prior
        externalEdgeLensModified(_external_edgelens_key, tree);
        }
    else
        {
        // Fixed tree topology, so there is a separate prior in the joint prior manager for each edge length parameter
        // Find all the entries corresponding to external edge lengths and change their prior means
        for(prior_distr_vect_t::const_iterator it = _distr_vect.begin(); it != _distr_vect.end(); ++it)
            {
            PriorDistrShPtr p = *it;
            if (p->_name.find("extedge_") != std::string::npos)
                {
                p->_univar_dist->SetMeanAndVariance(new_mean, new_mean*new_mean);
                univariateModifiedNoDebugCheck(p->_name, p->_univar_value);
                }
            }

        // Now perform the debug prior check suppressed above
        PHYCAS_ASSERT(debugCheckLogJointPrior("_replaceExternalEdgelenPriorMean"));
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current external edge length prior mean with supplied value `new_external_mean' and replaces current
|   internal edge length prior mean with supplied value `new_internal_mean'. Calls externalEdgeLensModified and
|   internalEdgeLensModified after replacing prior means to force recalculation of the external and internal edge length
|   prior.
*/
void JointPriorManager::_replaceExternalAndInternalEdgelenPriorMean(
  double new_external_mean,  /**< is the new external edge length prior mean */
  double new_internal_mean,  /**< is the new internal edge length prior mean */
  TreeShPtr tree)   /**< is the tree */
    {
    // Replace the mean of the prior on external edge lengths
    if (_checkDistrMapKey(_external_edgelens_key) && _checkDistrMapKey(_internal_edgelens_key))
        {
        PriorDistrShPtr p1 = _distr_map[_external_edgelens_key];
        p1->_univar_dist->SetMeanAndVariance(new_external_mean, new_external_mean*new_external_mean);

        // Replace the mean of the prior on internal edge lengths
        PriorDistrShPtr p2 = _distr_map[_internal_edgelens_key];
        p2->_univar_dist->SetMeanAndVariance(new_internal_mean, new_internal_mean*new_internal_mean);

        // Force recalculation of the internal edge length prior
        externalAndInternalEdgeLensModified(_external_edgelens_key, _internal_edgelens_key, tree);
        }
    else
        {
        // Fixed tree topology, so there is a separate prior in the joint prior manager for each edge length parameter
        // Find all the entries corresponding to external edge lengths and change their prior means
        for(prior_distr_vect_t::const_iterator it = _distr_vect.begin(); it != _distr_vect.end(); ++it)
            {
            PriorDistrShPtr p = *it;
            if (p->_name.find("extedge_") != std::string::npos)
                {
                p->_univar_dist->SetMeanAndVariance(new_external_mean, new_external_mean*new_external_mean);
                univariateModifiedNoDebugCheck(p->_name, p->_univar_value);
                }
            else if (p->_name.find("intedge_") != std::string::npos)
                {
                p->_univar_dist->SetMeanAndVariance(new_internal_mean, new_internal_mean*new_internal_mean);
                univariateModifiedNoDebugCheck(p->_name, p->_univar_value);
                }
            }

        // Now perform the debug prior check suppressed above
        PHYCAS_ASSERT(debugCheckLogJointPrior("_replaceExternalAndInternalEdgelenPriorMean"));
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current edge length hyperparameter value with supplied value `v'. Recalculates `_current_density'
|   for the edge length hyperparameter and and forces recalculation of the internal and/or external edge length
|   priors. Updates `_log_joint_prior'.
*/
double JointPriorManager::edgeLenHyperparamModified(
  std::string name,
  TreeShPtr t,
  double v)
    {
    PHYCAS_ASSERT(name == _internal_hyper_key || name == _external_hyper_key || name == _universal_hyper_key);
    PHYCAS_ASSERT(_checkDistrMapKey(name));
    if (name == _internal_hyper_key)
        _replaceInternalEdgelenPriorMean(v, t);
    else if (name == _external_hyper_key)
        _replaceExternalEdgelenPriorMean(v, t);
    else if (name == _universal_hyper_key)
        {
        _replaceExternalAndInternalEdgelenPriorMean(v, v, t);
        }
    PriorDistrShPtr p = _distr_map[name];
    p->_univar_value = v;
    double prev_logd = p->_current_density;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("edgeLenHyperparamModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds a univariate probability distribution to `_distr_map' and `_distr_vect' that applies to an edge length
|   hyperparameter.
*/
void JointPriorManager::addEdgelenHyperprior(
  std::string keystr,     /**< is the string used to identify this distribution */
  ProbDistShPtr dist,     /**< is the univariate distribution that applies to the universal edge length hyperparameter */
  TreeShPtr t,            /**< is the tree (needed for updating edge length prior(s) after setting initial value of hyperparameter) */
  double initial_value)   /**< is the starting value */
    {
    //std::cerr << "JointPriorManager::addEdgelenHyperprior: " << keystr << std::endl;
    PriorDistrShPtr p = PriorDistrShPtr(new PriorDistribution());
    p->initUnivariate(keystr, dist, initial_value);
    _distr_map[keystr] = p;
    _distr_vect.push_back(p);
    _log_joint_prior += p->_recalculateLogDensity();

    edgeLenHyperparamModified(keystr, t, initial_value);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Replaces current univariate value with supplied vector `v'. Recalculates `_current_density' for the parameter and
|   updates `_log_joint_prior'.
*/
double JointPriorManager::topologyModified(
  std::string name,
  TreeShPtr t)
    {
    if (IGNORING_TOPOLOGY_PRIOR)
        {
        return 0.0;
        }

    PHYCAS_ASSERT(_checkDistrMapKey(name));
    PriorDistrShPtr p = _distr_map[name];
    double prev_logd = p->_current_density;
    p->_tree = t;
    double new_logd = p->_recalculateLogDensity();
    _log_joint_prior += (new_logd - prev_logd);
    PHYCAS_ASSERT(debugCheckLogJointPrior(boost::str(boost::format("topologyModified (%s)") % name)));
    return new_logd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns true if a tree length prior is being used (specifically, returns true if "tree_length" is a key in
|   `_distr_map').
*/
bool JointPriorManager::isTreeLengthPrior() const
    {
    return _checkDistrMapKey("tree_length");
    }

}	// namespace phycas
