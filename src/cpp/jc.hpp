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

#if ! defined(JC_HPP)
#define JC_HPP

#include "states_patterns.hpp"
#include "model.hpp"
#include <boost/shared_ptr.hpp>

namespace phycas{

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the Jukes and Cantor (1969) model.
*/
class JC: public Model
{
public:
    JC();
    ~JC()
    {
        //std::cerr << "JC dying..." << std::endl;
    }

    virtual std::string             getModelName() const;
    double                          calcUniformizationLambda() const;
    void                            calcPMat(double * * pMat, double edgeLength) const;
    double                          calcLMat(double * * lMat) const;
    double                          calcUMat(double * * uMat) const;
    virtual std::string             paramHeader() const;
    virtual std::string             paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;

    virtual void                    appendFreeParamNames(std::vector<std::string> & names, std::string prefix = "") const;
    virtual void                    appendParamNames(std::vector<std::string> & names, std::string prefix = "") const;
    virtual void                    appendUntransformedParamValues(std::vector<double> & values) const;
    virtual void                    appendTransformedParamValues(std::vector<double> & values) const;
    virtual bool                    setParamValueFromTransformed(std::string parameter_name, double transformed_value, TreeShPtr tree);
    virtual double                  calcLogDetJacobian() const;
};

typedef boost::shared_ptr<JC> JCShPtr;

}

#endif
