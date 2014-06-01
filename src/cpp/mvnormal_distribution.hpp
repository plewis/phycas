/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2009 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#if !defined(MVNORMAL_DISTRIBUTION_HPP)
#define MVNORMAL_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include "ncl/nxsallocatematrix.h"
#include "multivariate_probability_distribution.hpp"

extern "C"
{
#include "linalg.h"
}

namespace phycas
{

class MVNormalDistribution : public MultivariateProbabilityDistribution
	{
	public:
											MVNormalDistribution();
											MVNormalDistribution(const std::vector<double> & meanvec, const std::vector<double> & varcovmat);
                        					MVNormalDistribution(const MVNormalDistribution & other);
											~MVNormalDistribution() {}
	
        MVNormalDistribution * 				cloneAndSetLot(Lot * other) const;
        MVNormalDistribution * 				Clone() const;
        
        bool								Fit(unsigned nrows, unsigned ncols, std::vector<double> data);
        
		void								SetLot(Lot * other);
		void								ResetLot();
		void								SetSeed(unsigned rnseed);

		bool								IsDiscrete() const;
		std::string 						GetDistributionName() const;
		std::string 						GetDistributionDescription() const;
		std::string 						GetDescriptionForPython() const;
		std::vector<double>					GetMean() const;
		std::vector<double> 				GetVar() const;
		std::vector<double> 				GetStdDev() const;
		std::vector<double>					Sample() const;
		double								ApproxCDF(const std::vector<double> &x, unsigned nsamples = 10000) const;
		double								GetLnPDF(const std::vector<double> &x) const;
		double								GetRelativeLnPDF(const std::vector<double> &x) const;
		void 								SetMeanAndVariance(const std::vector<double> & m, const std::vector<double> & v);
		std::vector<double>					GetVarCovarMatrix();

		unsigned							GetNParams() const;
		
		// non-virtual member functions
		bool 								InitMVNorm();
		std::string 						DebugMVNorm() const;
		void 								MakeVarCovSymmetric();

    private:
    
		// parameters
		std::vector<double>					mean;
		ScopedTwoDMatrix<double>			varcov;
		
		// used for drawing multivariate normal deviates
		ScopedTwoDMatrix<double>			L;
		ScopedTwoDMatrix<double>			invVarCov;
		NormalDistribution					stdnorm;
		double								logDetVarCov;
		double								nLogTwoPi;
	};

} // namespace phycas

#endif
