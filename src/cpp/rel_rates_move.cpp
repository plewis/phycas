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
#include <boost/cast.hpp>       //temporary!
#include <boost/math/special_functions/fpclassify.hpp>
#include "probability_distribution.hpp"
#include "relative_rate_distribution.hpp"
#include "model.hpp"
#include "basic_tree_node.hpp"
#include "tree_likelihood.hpp"
#include "xlikelihood.hpp"
#include "mcmc_chain_manager.hpp"
#include "dirichlet_move.hpp"
#include "rel_rates_move.hpp"
#include "basic_tree.hpp"
#include "tree_manip.hpp"
#include "gtr.hpp"
#include "hky.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
RelRatesMove::RelRatesMove() : DirichletMove()
	{
	dim = 6;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the relative rates of the associated GTR model to those in the supplied vector `v'.
*/
void RelRatesMove::sendCurrValuesToModel(const double_vect_t & v)
	{
	PHYCAS_ASSERT(dim == v.size());
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	PHYCAS_ASSERT(gtr_model);
	gtr_model->setRelRates(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the supplied vector `v'.
*/
void RelRatesMove::getCurrValuesFromModel(double_vect_t & v) const
	{
	PHYCAS_ASSERT(dim > 0);
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	if (gtr_model)
		{
		const std::vector<double> & rrates = gtr_model->getRelRates();
		PHYCAS_ASSERT(dim == rrates.size());
		v.resize(rrates.size());
		std::copy(rrates.begin(), rrates.end(), v.begin());
    	}
    else
    	{
    	v.assign(dim, 1.0/(double)dim);
    	}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, returning them as an anonymous vector.
*/
double_vect_t RelRatesMove::listCurrValuesFromModel()
	{
	PHYCAS_ASSERT(dim > 0);
	double_vect_t v(dim);
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	if (gtr_model)
		{
		const double_vect_t & rrates = gtr_model->getRelRates();
		PHYCAS_ASSERT(dim == rrates.size());
		std::copy(rrates.begin(), rrates.end(), v.begin());
    	}
    else
    	{
    	v.assign(dim, 1.0/(double)dim);
    	}
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the data member `orig_params'.
*/
void RelRatesMove::getParams()
	{
	getCurrValuesFromModel(orig_params);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the relative rates in the model with those supplied in the vector `v'.
*/
void RelRatesMove::setParams(
  const std::vector<double> & v)    /*< is the vector of parameter values to send to the model */
	{
    GTR * gtr_model = dynamic_cast<GTR *>(model.get());
    gtr_model->setRelRates(v);
	}

