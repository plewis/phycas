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

#if ! defined(RECTANGULAR_MATRIX_HPP)
#define RECTANGULAR_MATRIX_HPP

//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include <boost/shared_ptr.hpp>
#include <vector>
#include "square_matrix.hpp"

namespace phycas
{

class RectangularMatrix
	{
	public:
										RectangularMatrix();
										RectangularMatrix(unsigned nrows, unsigned ncols, double value);
										RectangularMatrix(const RectangularMatrix & other);
										~RectangularMatrix();
		void							CreateMatrix(unsigned nrows, unsigned ncols, double value);
		void							ScalarMultiply(double scalar);
		void							SetToScalar(double scalar);
		std::pair<unsigned,unsigned>	GetDimensions() const;
		unsigned                        GetNRows() const;
		unsigned                        GetNCols() const;
		double * *						GetMatrixAsRawPointer() const;
		double *						operator[](unsigned i) const;
		void							MatrixToString(std::string & s, std::string fmt) const;
		void							Fill(double value);
        std::string                     GetStringRepresentation() const;
        std::string                     GetFormattedStringRepresentation(std::string fmt) const;
        double                          GetElement(unsigned i, unsigned j) const;
        void                            AddToElement(unsigned i, unsigned j, double v);
        void                            SetElement(unsigned i, unsigned j, double v);
        std::vector<double>             GetRow(unsigned i) const;
        void                            SetRow(unsigned i, const std::vector<double> & v);
        std::vector<double>             GetMatrix() const;
        void                            SetMatrix(unsigned nrows, unsigned ncols, std::vector<double>);
        std::vector<double>             GetMean();
        SquareMatrix *                  GetVarCovMatrix();

    protected:

        void                            ComputeMean();
        void                            ComputeVarCovMatrix();

	protected:
		static unsigned					_k;		//temporary
		unsigned						_id;    //temporary
		double * *						_m;		/**< the two-dimensional matrix of doubles */
		unsigned						_nrows;	/**< the row dimension of the matrix */
		unsigned						_ncols;	/**< the column dimension of the matrix */
        std::vector<double>             _mean;      /**< the mean vector */
        SquareMatrix *                  _varcov;    /**< the variance-covariance matrix */
        bool                            _dirty_mean; /**< true if _m has changed but _mean has not been updated */
        bool                            _dirty_varcov; /**< true if _m has changed but _varcov has not been updated */
	};

typedef boost::shared_ptr<RectangularMatrix> RectangularMatrixShPtr;


}	// namespace phycas
#endif
