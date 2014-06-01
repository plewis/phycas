/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "square_matrix.hpp"
#include "rectangular_matrix.hpp"

#include <boost/format.hpp>
#include <iostream> //temporary!
#include <iterator> //temporary!

#include "ncl/nxsallocatematrix.h"

namespace phycas
{

unsigned RectangularMatrix::_k = 0;	//temporary!

/*----------------------------------------------------------------------------------------------------------------------
|	RectangularMatrix default constructor. Sets `_nrows', '_ncols' and `_m' data members to 0.
*/
RectangularMatrix::RectangularMatrix()
  :	 _id(++_k),  _m(0), _nrows(0), _ncols(0), _varcov(0), _dirty_mean(true), _dirty_varcov(true)
	{
	//std::cerr << "====> constructing default RectangularMatrix " << _id << " <====" << (this) << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	RectangularMatrix constructor. Creates rectangular matrix `_m' with row dimension `_nrows', column dimension
|   '_ncols', and sets all elements of `_m' to `value'.
*/
RectangularMatrix::RectangularMatrix(
  unsigned nrows,   /**< the row dimension of the matrix */
  unsigned ncols,   /**< the column dimension of the matrix */
  double value)		/**< the value to which to set all elements */
  : _id(++_k), _m(0), _nrows(nrows), _ncols(ncols), _varcov(0), _dirty_mean(true), _dirty_varcov(true)
	{
	//std::cerr << "====> constructing RectangularMatrix " << _id << " of size " << _nrows << " rows, " << _ncols << " columns, with all values set to " << value << " <====" << (this) << std::endl;
	CreateMatrix(_nrows, _ncols, value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	RectangularMatrix copy constructor. This is necessary because the standard vector resize function in some
|   implementations of the standard template library creates only one element with the default constructor, then creates
|   the remaining elements using the copy constructor.
*/
RectangularMatrix::RectangularMatrix(
  const RectangularMatrix & other)	/**< is the RectangularMatrix to copy */
  : _id(++_k), _m(0), _nrows(other._nrows), _ncols(other._ncols), _varcov(0), _dirty_mean(true), _dirty_varcov(true)
  	{
	//std::cerr << "====> copy constructing RectangularMatrix " << _id << " from " << other._id << " <====" << (this) << std::endl;
	if (_nrows > 0 && _ncols > 0)
		{
		CreateMatrix(_nrows, _ncols, 0.0);
		double * pother = &other._m[0][0];
		double * p = &_m[0][0];
		unsigned last = _nrows*_ncols;
		for (unsigned i = 0; i < last; ++i)
			*p++ = *pother++;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	RectangularMatrix destructor. If `_m' exists, deletes it.
*/
RectangularMatrix::~RectangularMatrix()
	{
	//std::cerr << "====> destroying RectangularMatrix " << _id << " <====" << (this) << std::endl;
	if (_m)
		DeleteTwoDArray<double>(_m);
	if (_varcov)
		delete _varcov;
	_m = 0;
	_nrows = 0;
	_ncols = 0;
	_varcov = 0;
    _dirty_mean = true;
    _dirty_varcov = true;
	_id = UINT_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for `_m' data member.
*/
void RectangularMatrix::CreateMatrix(
  unsigned nrows,   /**< is the row dimension of the rectangular matrix to be created */
  unsigned ncols,   /**< is the column dimension of the rectangular matrix to be created */
  double value)		/**< is the value to which each element will be set */
	{
	//std::cerr << "----> RectangularMatrix::CreateMatrix " << _id << ", _nrows = " << nrows << ", _ncols = " << ncols << ", value = " << value << " <----" << std::endl;
    _dirty_mean = true;
    _dirty_varcov = true;
    _nrows = nrows;
    _ncols = ncols;
	if (_m)
		DeleteTwoDArray<double>(_m);
	_m = 0;
	if (_varcov)
		delete _varcov;
	_varcov = 0;
	if (_nrows > 0 && _ncols > 0)
		{
		_m = NewTwoDArray<double>(_nrows, _ncols);
		Fill(value);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the dimensions of the data member `_m' as a std::pair<unsigned,unsigned> (number of rows, number of columns).
*/
std::pair<unsigned, unsigned> RectangularMatrix::GetDimensions() const
	{
	return std::pair<unsigned, unsigned>(_nrows, _ncols);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of rows of the data member `_m'.
*/
unsigned RectangularMatrix::GetNRows() const
	{
	return _nrows;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of columns of the data member `_m'.
*/
unsigned RectangularMatrix::GetNCols() const
	{
	return _ncols;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that simply returns `_m'. This allows `_m' to be passed to functions that expect a two-dimensional
|	array of doubles, but keep in mind that this is somewhat unsafe.
*/
double * * RectangularMatrix::GetMatrixAsRawPointer() const
	{
	return _m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Operator that allows access to this RectangularMatrix object as if it were a two-dimensional array of doubles.
*/
double * RectangularMatrix::operator[](
  unsigned row) const
	{
	PHYCAS_ASSERT(row < _nrows);
	return _m[row];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns element at row `i', column `'j of the matrix.
*/
double RectangularMatrix::GetElement(unsigned i, unsigned j) const
	{
	PHYCAS_ASSERT(i < _nrows);
	PHYCAS_ASSERT(j < _ncols);
    return _m[i][j];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets element at row `i', column `j' of the matrix to the value `v'.
*/
void RectangularMatrix::SetElement(unsigned i, unsigned j, double v)
	{
	PHYCAS_ASSERT(i < _nrows);
	PHYCAS_ASSERT(j < _ncols);
    _m[i][j] = v;
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets element at row `i', column `j' of the matrix to the current value plus `v'.
*/
void RectangularMatrix::AddToElement(unsigned i, unsigned j, double v)
	{
	PHYCAS_ASSERT(i < _nrows);
	PHYCAS_ASSERT(j < _ncols);
    _m[i][j] += v;
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns row `i' of the matrix as a vector.
*/
std::vector<double> RectangularMatrix::GetRow(unsigned i) const
	{
	PHYCAS_ASSERT(i < _nrows);
    double * first = &_m[i][0];
    double * last = first + _ncols;
    return std::vector<double>(first, last);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets row `i' of the matrix to the values in the vector `v'.
*/
void RectangularMatrix::SetRow(unsigned i, const std::vector<double> & v)
	{
	PHYCAS_ASSERT(i < _nrows);
	PHYCAS_ASSERT(v.size() == _ncols);
    unsigned j = 0;
    for (std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it)
        {
        double value = *it;
        _m[i][j++] = value;
        }
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entire matrix as a vector by row.
*/
std::vector<double> RectangularMatrix::GetMatrix() const
	{
    std::vector<double> v;
    double * ptr = (double *)&_m[0];
    for (unsigned i = 0; i < _nrows*_ncols; ++i)
        v.push_back(ptr[i]);
    return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies supplied vector `v' into this matrix. Vector `v' is expected to have length `nrows'*`ncols', and matrix is
|   read from `v' by row. The data member `_nrows' is set equal to `nrows', and the data member '_ncols' is set equal to
|   'ncols'.
*/
void RectangularMatrix::SetMatrix(unsigned nrows, unsigned ncols, std::vector<double> v)
	{
	PHYCAS_ASSERT(nrows*ncols == (unsigned)v.size());
    CreateMatrix(nrows, ncols, 0.0);
    double * ptr = (double *)&_m[0][0];
    for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
        {
        double val = *it;
        *ptr++ = val;
        }
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fills all cells of matrix `_m' with `value'.
*/
void RectangularMatrix::Fill(
  double value) /**< is the value to which every element in matrix is to be set */
	{
    std::pair<unsigned,unsigned> dim = GetDimensions();
	unsigned last = dim.first*dim.second;
	double * p = &_m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ = value;
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes '_nrows' > 0, '_ncols' > 0, and `_m' exists. Multiplies each element of `_m' by the supplied `scalar'.
*/
void RectangularMatrix::ScalarMultiply(
  double scalar)	/**< is the scalar value multiplied by each element */
	{
	//std::cerr << "----> RectangularMatrix::ScalarMultiply " << _id << ", scalar = " << scalar << " <----" << std::endl;
	PHYCAS_ASSERT(_nrows > 0);
	PHYCAS_ASSERT(_ncols > 0);
	unsigned last = _nrows*_ncols;
	double * p = &_m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ *= scalar;
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assumes '_nrows' > 0, '_ncols' > 0, and `_m' exists. Sets each element of `_m' to the supplied `scalar'.
*/
void RectangularMatrix::SetToScalar(
  double scalar)	/**< is the scalar value to which each element is set */
	{
	//std::cerr << "----> RectangularMatrix::SetToScalar " << _id << ", scalar = " << scalar << " <----" << std::endl;
	PHYCAS_ASSERT(_nrows > 0);
	PHYCAS_ASSERT(_ncols > 0);
	unsigned last = _nrows*_ncols;
	double * p = &_m[0][0];
	for (unsigned i = 0; i < last; ++i)
		*p++ = scalar;
    _dirty_mean = true;
    _dirty_varcov = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves a representation of the matrix to the supplied string for use in debugging. The caller should supply a format
|	string suitable for representing a single value of the matrix: e.g. "%12.5f\t".
*/
void RectangularMatrix::MatrixToString(
  std::string & s,		 /**< the string to which a representation will be saved */
  std::string fmt) const /**< the format string */
	{
    std::pair<unsigned,unsigned> dim = GetDimensions();
	//s += boost::str(boost::format("This is RectangularMatrix %d\n") % _id);
	for (unsigned i = 0; i < dim.first; ++i)
		{
		for (unsigned j = 0; j < dim.second; ++j)
			{
			s += str(boost::format(fmt) % _m[i][j]);
			}
		s += '\n';
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns string representation of the matrix.
*/
std::string RectangularMatrix::GetStringRepresentation() const
	{
    std::string s;
    std::string fmt = "\t%g";
    MatrixToString(s, fmt);
    return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns string representation of the matrix.
*/
std::string RectangularMatrix::GetFormattedStringRepresentation(std::string fmt = "%g") const
	{
    std::string s;
    MatrixToString(s, fmt);
    return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the data member `_mean'. Calls ComputeMean() first if `_mean' is not up-to-date.
*/
std::vector<double> RectangularMatrix::GetMean()
	{
    if (_dirty_mean || _mean.size() != _ncols)
        ComputeMean();
    return _mean;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a copy of the data member `_varcov'. Calls ComputeVarCovMatrix() first if  `_varcov' is not up-to-date.
*/
SquareMatrix * RectangularMatrix::GetVarCovMatrix()
	{
    if (_dirty_varcov || !_varcov || _varcov->GetDimension() != _ncols)
        ComputeVarCovMatrix();
    PHYCAS_ASSERT(_varcov);
    return new SquareMatrix(*_varcov);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes the mean vector `_mean'.
*/
void RectangularMatrix::ComputeMean()
	{
    if (_nrows == 0)
        {
        _mean.clear();
        _dirty_mean = false;
        return;
        }

    // compute sample mean vector
    double n = (double)(_nrows);
    _mean.assign(_ncols, 0.0);
    for (unsigned i = 0; i < _nrows; ++i)
        {
        for (unsigned j = 0; j < _ncols; ++j)
            {
            double v = _m[i][j];
            _mean[j] += v/n;
            }
        }

    _dirty_mean = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes the variance-covariance matrix `_varcov'.
*/
void RectangularMatrix::ComputeVarCovMatrix()
	{
    if (_nrows == 0)
        {
        _varcov->Clear();
        }
    else if (_nrows == 1)
        {
        _varcov->SetToScalar(0.0);
        }
    else
        {
        // compute sample mean vector
        double n = (double)(_nrows);
        //_mean.assign(_ncols, 0.0);
        //std::vector<double> sample_sum(_ncols, 0.0);
        //std::vector<double> sample_ss(_ncols, 0.0);
        //for (unsigned i = 0; i < _nrows; ++i)
        //    {
        //    for (unsigned j = 0; j < _ncols; ++j)
        //        {
        //        double v = _m[i][j];
        //        _mean[j] += v/n;
        //        sample_sum[j] += v;
        //        sample_ss[j] += v*v;
        //        }
        //    }

        if (_dirty_mean)
            ComputeMean();

        // compute sample variance-covariance matrix
        if (_varcov && _varcov->GetDimension() == _ncols)
            _varcov->SetToScalar(0.0);
        else
            {
            delete _varcov;
            _varcov = new SquareMatrix(_ncols, 0.0);
            }

        for (unsigned i = 0; i < _ncols; ++i)
            {
            for (unsigned j = i; j < _ncols; ++j)
                {
                for (unsigned k = 0; k < _nrows; ++k)
                    {
                    double tmp = (_m[k][i] - _mean[i])*(_m[k][j] - _mean[j])/(n - 1.0);
                    _varcov->AddToElement(i, j, tmp);
                    }
                if (j != i)
                    {
                    double tmp = _varcov->GetElement(i, j);
                    _varcov->SetElement(j, i, tmp);
                    }
                }
            }
        }
    _dirty_varcov = false;
	}

} // namespace phycas
