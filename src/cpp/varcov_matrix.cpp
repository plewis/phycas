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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include "ncl/nxsallocatematrix.h"
#include "varcov_matrix.hpp"
#include "xlikelihood.hpp"
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|   Constructs an empty VarCovMatrix object.
*/
VarCovMatrix::VarCovMatrix()
  : dimension(0), vcmat(0), w(0), z(0), fv(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Destructor.
*/
VarCovMatrix::~VarCovMatrix()
	{
	clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns object to its just-constructed state.
*/
void VarCovMatrix::clear()
	{
	dimension = 0;

	if (vcmat)
		{
		DeleteTwoDArray<double>(vcmat);
		vcmat = 0;
		}

	if (w)
		{
		delete [] w;
		w = 0;
		}

	if (fv)
		{
		delete [] fv;
		fv = 0;
		}

	if (z)
		{
		DeleteTwoDArray<double>(z);
		z = 0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Generates a fixed 3x3 variance covariance matrix and computes and shows its eigenvalues and eigenvectors.
*/
void VarCovMatrix::example()
	{
    // Variance-covariance matrix for testing:
    //  1  2  3
    //  2  1  4
    //  3  4  1
    //
    // Eigenvalues:     7.07467    -3.18788    -0.886791
    //
    // Eigenvectors:   -0.505785   -0.255232    0.824038
    //                 -0.584374   -0.601302   -0.544925
    //                 -0.634577    0.757161   -0.154979
    redimension(3);
    vcmat[0][0] = 1.0;
    vcmat[0][1] = 2.0;
    vcmat[0][2] = 3.0;
    vcmat[1][0] = 2.0;
    vcmat[1][1] = 1.0;
    vcmat[1][2] = 4.0;
    vcmat[2][0] = 3.0;
    vcmat[2][1] = 4.0;
    vcmat[2][2] = 1.0;

    calcEigensystem();

    std::cerr << "\nVariance-Covariance matrix:" << std::endl;
    std::cerr << showVarCovMatrix() << std::endl;

    std::cerr << "\nEigenvalues:" << std::endl;
    std::vector<double> eval = getEigenValues();
    std::copy(eval.begin(), eval.end(), std::ostream_iterator<double>(std::cerr, " "));

    std::cerr << "\nEigenvectors:" << std::endl;
    //std::vector<double> evec = getEigenVectors();
    //std::copy(evec.begin(), evec.end(), std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << showEigenVectorsAsMatrix() << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears object and then sets `vcmat' to the supplied `flatvcmat'. Assumes that the supplied 1-dimensional matrix is
|   a 2-dimensional variance-covariance matrix stored by rows. Enforces symmetry by averaging upper and lower triangle
|   entries.
*/
void VarCovMatrix::setVarCovMatrix(std::vector<double> flatvcmat)
    {
    unsigned sz = (unsigned)flatvcmat.size();
    unsigned d = (unsigned)sqrt(sz);
    assert(d*d == sz);
    redimension(d);
    unsigned row = 0;
    unsigned col = 0;
    BOOST_FOREACH(double x, flatvcmat)
        {
        vcmat[row][col] = x;
        col++;
        if (col == d)
            {
            col = 0;
            row++;
            }
        }

    calcEigensystem();

    std::cerr << "\nVariance-Covariance matrix:" << std::endl;
    std::cerr << showVarCovMatrix() << std::endl;

    std::cerr << "\nEigenvalues:" << std::endl;
    std::vector<double> eval = getEigenValues();
    std::copy(eval.begin(), eval.end(), std::ostream_iterator<double>(std::cerr, " "));

    std::cerr << "\nEigenvectors:" << std::endl;
    std::cerr << showEigenVectorsAsMatrix() << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `dimension', which is the number of rows (and
|	columns) in the variance-covariance matrix.
*/
unsigned VarCovMatrix::getDimension()
	{
	return dimension;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reallocates all arrays and vectors (except `rr', `pi' and `sqrtPi') according to the value of `new_dim' and sets
|	`dimension' equal to `new_dim' before returning. Assumes `new_dim' is greater than zero. Note that `vcmat', `w', `z'
|	and `fv' are allocated but not filled with any values when this function returns.
*/
void VarCovMatrix::redimension(unsigned new_dim)
	{
	PHYCAS_ASSERT(new_dim > 0);
    clear();
	dimension = new_dim;

	// Create and fill vcmat with default (Mk model) values
	vcmat = NewTwoDArray<double>(dimension, dimension);
	PHYCAS_ASSERT(vcmat);

	// Create w (array of eigenvalues)
	w = new double[dimension];

	// Create z (two-dimensional array of eigenvectors)
	z = NewTwoDArray<double>(dimension, dimension);

	// Create fv (workspace used by EigenRealSymmetric)
	fv = new double[dimension];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes eigenvalues and eigenvectors corresponding to the variance-covariance matrix.
*/
void VarCovMatrix::calcEigensystem()
	{
	// Calculate eigenvalues (w) and eigenvectors (z)
	int err_code = EigenRealSymmetric(dimension, vcmat, w, z, fv);

	if (err_code != 0)
		{
		clear();
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues of the variance-covariance matrix");
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvalues in the data member `w' as a vector.
*/
std::vector<double> VarCovMatrix::getEigenValues()
	{
	std::vector<double> v(dimension, 0.0);
	double * w_begin = &w[0];
	double * w_end = w_begin + dimension;
	std::copy(w_begin, w_end, v.begin());
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvectors in the data member `z' as a NumArray. Calls recalcQMatrix first to ensure that `z' is up
|	to date.
*/
std::vector<double> VarCovMatrix::getEigenVectors()
	{
	std::vector<double> p;
	flattenTwoDMatrix(p, z, dimension);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats `vcmat' nicely for output.
*/
std::string VarCovMatrix::showVarCovMatrix()
	{
	unsigned flat_length = dimension*dimension;
	double * vcmat_begin = &vcmat[0][0];
	double * vcmat_end   = vcmat_begin + flat_length;
	return showMatrixImpl(vcmat_begin);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats the matrix of eigenvectors nicely for output.
*/
std::string VarCovMatrix::showEigenVectorsAsMatrix()
	{
	unsigned flat_length = dimension*dimension;
	double * z_begin = &z[0][0];
	double * z_end   = z_begin + flat_length;
	return showMatrixImpl(z_begin);
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
std::string VarCovMatrix::showMatrixImpl(const double * q) const
	{
	unsigned i, j;

	// Output one column for row labels and a label for every column of vcmat
	std::string s = str(boost::format("%12s ") % " ");
	for (i = 0; i < dimension; ++i)
		{
		s += str(boost::format("%12d ") % i);
		}
	s += "\n";

	// Output rows of q
	double * pq = (double *)q;	//PELIGROSO
	for (i = 0; i < dimension; ++i)
		{
		s += str(boost::format("%12d ") % i);
		for (j = 0; j < dimension; ++j)
			{
			s += str(boost::format("%12.5f ") % *pq++);
			}
		s += "\n";
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores a flattened version of the supplied 2-dimensional array `twoDarr', storing the result in the supplied std::vector<double>
|	reference variable `p'. The supplied `twoDarr' should be laid out so that rows occupy contiguous memory.
*/
void VarCovMatrix::flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const
	{
	unsigned flat_length = dim*dim;
	p.resize(flat_length);
	double * twoD_begin = &twoDarr[0][0];
	double * twoD_end   = twoD_begin + flat_length;
	std::copy(twoD_begin, twoD_end, p.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entries in the data member `vcmat' as a single dimensional vector.
*/
std::vector<double> VarCovMatrix::getVarCovMatrix()
	{
	std::vector<double> p;
	flattenTwoDMatrix(p, vcmat, dimension);
	return p;
	}


