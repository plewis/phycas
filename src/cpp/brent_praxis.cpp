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

#if 0

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include <cmath>
#include "brent_praxis.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the Init() member function.
*/
BrentPraxis::BrentPraxis()
	{
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the Init() member function after setting `r' to the supplied `rnd' and `func' to the supplied `f'.
*/
BrentPraxis::BrentPraxis(
  LotShPtr rnd,             /**< is the random number generator object to use */
  FuncToMinimizeShPtr f)    /**< is the function to minimize */
  : _func(f), _r(rnd)
	{
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to do.
*/
BrentPraxis::~BrentPraxis()
	{
    if (_v)
        DeleteTwoDArray<double>(_v)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by constructors to initialize the object.
*/
void BrentPraxis::Init()
	{
    _n = 0
    _v = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the functor `func' to the supplied FuncToMinimize object representing the function to be minimized. Calls 
|   Init() to reset this object to its just-constructed state.
*/
void BrentPraxis::AttachFunc(
  FuncToMinimize f)
	{
	_func = f;
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number generator pointer `r' to the supplied random number generator pointer `rnd'. Calls Init() to
|	reset this object to its just-constructed state.
*/
void BrentPraxis::AttachRandomNumberGenerator(
  LotShPtr rnd)     /**< is the random number generator to attach */
	{
	_r = rnd;
	Init();
	}
    
/*----------------------------------------------------------------------------------------------------------------------
|	Sort eigenvalues _d in descending order, and keep corresponding eigenvectors (columns of _v) in sync with _d.
*/
void BrentPraxis::Sort()
    {
    int k, i, j;
    double s;

    for (i = 0; i < _n-1; ++i) 
        {
        // find largest of the downstream values of _d that are greater than _d[i]
        k = i; 
        s = _d[i];
        for (j = i+1; j < _n; ++j) 
            {
            if (_d[j] > s) 
                {
                k = j;
                s = _d[j];
                }
            }
        // if _d[k] is largest of the downstream values larger than _d[i], 
        // swap _d[i] and _d[k] and swap ith. and kth. columns of _v
        if (k > i) 
            {
            _d[k] = _d[i];
            _d[i] = s;
            for (j = 0; j < _n; ++j) 
                {
                s = _v[j][i];
                _v[j][i] = _v[j][k];
                _v[j][k] = s;
                }
            }
        }
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
double BrentPraxis::flin(double l, int j)
    {
    int i;
    std::vector<double> tflin(_n, 0.0);

    if (j != -1) /* linear search */
        {		
        for (i = 0; i < _n; ++i)
            tflin[i] = _minimum_point[i] + l*v[i][j];
        }   
    else    /* search along parabolic space curve */
        {			
        qa = l*(l-qd1)/(qd0*(qd0+qd1));
        qb = (l+qd0)*(qd1-l)/(qd0*qd1);
        qc = l*(l+qd0)/(qd1*(qd0+qd1));
        for (i = 0; i < _n; ++i)
            tflin[i] = qa*q0[i]+qb*_minimum_point[i]+qc*q1[i];
        }
    nf++;
    return (*_func)(tflin);
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
double BrentPraxis::min(int j, int nits, double *d2, double *x1, double f1, int fk)
    {
    int k, i, dz;
    double x2, xm, f0, f2, fm, d1, t2, s, sf1, sx1;

    sf1 = f1; 
    sx1 = *x1;
    k = 0; 
    xm = 0.0; 
    fm = f0 = fx; 
    dz = *d2 < macheps;
    
    // find step size
    // let s be the length of the _minimum_point vector
    s = 0.0;
    for (i = 0; i < _n; ++i) 
        s += _minimum_point[i]*_minimum_point[i];
    s = sqrt(s);
    
    if (dz)
        t2 = m4*sqrt(fabs(fx)/dmin + s*ldt) + m2*ldt;
    else
        t2 = m4*sqrt(fabs(fx)/(*d2) + s*ldt) + m2*ldt;
    s = s*m4 + t;
    if (dz && t2 > s) 
        t2 = s;
    if (t2 < small) 
        t2 = small;
    if (t2 > 0.01*h) 
        t2 = 0.01 * h;
    if (fk && f1 <= fm) 
        {
        xm = *x1;
        fm = f1;
        }
    if (!fk || fabs(*x1) < t2) 
        {
        *x1 = (*x1 > 0 ? t2 : -t2);
        f1 = flin(*x1, j);
        }
    if (f1 <= fm) 
        {
        xm = *x1;
        fm = f1;
        }
        
    next:
    if (dz) 
        {
        x2 = (f0 < f1 ? -(*x1) : 2*(*x1));
        f2 = flin(x2, j);
        if (f2 <= fm) 
            {
            xm = x2;
            fm = f2;
            }
        *d2 = (x2*(f1-f0) - (*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
        }
    d1 = (f1-f0)/(*x1) - *x1**d2; dz = 1;
    if (*d2 <= small) 
        {
        x2 = (d1 < 0 ? h : -h);
        }
    else 
        {
        x2 = - 0.5*d1/(*d2);
        }
    if (fabs(x2) > h)
        x2 = (x2 > 0 ? h : -h);
        
    tryme:
    f2 = flin(x2, j);
    if ((k < nits) && (f2 > f0)) 
        {
        k++;
        if ((f0 < f1) && (*x1*x2 > 0.0))
            goto next;
        x2 *= 0.5;
        goto tryme;
        }
    nl++;
    if (f2 > fm) 
        x2 = xm; 
    else 
        fm = f2;
    if (fabs(x2*(x2-*x1)) > small) 
        {
        *d2 = (x2*(f1-f0) - *x1*(fm-f0))/(*x1*x2*(*x1-x2));
        }
    else 
        {
        if (k > 0) *d2 = 0;
        }
    if (*d2 <= small) 
        *d2 = small;
    *x1 = x2; fx = fm;
    if (sf1 < fx) 
        {
        fx = sf1;
        *x1 = sx1;
        }
    if (j != -1)
        for (i = 0; i < _n; ++i)
            _minimum_point[i] += (*x1)*v[i][j];
}
    

/*----------------------------------------------------------------------------------------------------------------------
|	This is where all the work takes place. Starts with a point defined by the supplied vector `starting_point' and 
|   returns the function value at the achieved minimum. The final point is stored in the data member `minimum_point'.
*/
double BrentPraxis::Minimize(const std::vector<double> starting_point)
    {
    // init global extern variables and parameters
    double macheps = 1.0e-8;    //originally EPSILON
    double h = step; 
    double t = 1.0e-16;         //originally tol=SQREPSILON
    
    // copy starting_point to minimum_point
    _n = (unsigned)starting_point.size();
    _minimum_point.resize(_n);
    std::copy(starting_point.begin(), starting_point.end(), _minimum_point.begin());
    
    double      small   = macheps*macheps; 
    double      vsmall  = small*small;
    double      large   = 1.0/small; 
    double      vlarge  = 1.0/vsmall;
    double      m2      = sqrt(macheps); 
    double      m4      = sqrt(m2);
    double      ldfac   = (_illc ? 0.1 : 0.01);
    int         nl      = 0;
    int         kt      = 0; 
    int         nf      = 1; 
    double      fx      = (*_func)(_minimum_point); 
    double      qf1     = fx;
    double      t2      = small + fabs(t); 
    t = t2; 
    double      dmin    = small;
    
    _v = new NewTwoDArray<double>(_n, _n); //PELIGROSO
    std::vector<double> q0(_n, 0.0);
    std::vector<double> q1(_n, 0.0);
    _d.asign(_n, 0.0);
    std::vector<double> y(_n, 0.0);
    std::vector<double> z(_n, 0.0);

    if (h < 100.0*t)
        h = 100.0*t;
    ldt = h;
    for (i = 0; i < _n; ++i)
        for (j = 0; j < _n; j++)
            _v[i][j] = (i == j ? 1.0 : 0.0);
    _d[0] = 0.0; 
    double qd0 = 0.0;
    for (i = 0; i < _n; ++i) 
        q1[i] = _minimum_point[i];
    if (_prin > 1) 
        {
        std::cerr << "\n------------- enter function praxis -----------\n";
        std::cerr << "... current parameter settings ...\n";
        std::cerr << boost::str(boost::format("... scaling ... %20.10e\n") % _scbd);
        std::cerr << boost::str(boost::format("...   tol   ... %20.10e\n") % t);
        std::cerr << boost::str(boost::format("... maxstep ... %20.10e\n") % h);
        std::cerr << boost::str(boost::format("...   illc  ... %20u\n") % _illc);
        std::cerr << boost::str(boost::format("...   ktm   ... %20u\n") % _ktm);
        std::cerr << boost::str(boost::format("... maxfun  ... %20u\n") % _maxfun);
        }
    if (_prin > 0)
        std::cerr << std::endl;

    bool done = false;
    while (!done)
        {
        sf = _d[0];
        s = _d[0] = 0.0;

        // minimize along first direction
        min(0, 2, &_d[0], &s, fx, 0);
        if (s <= 0.0)
            for (i=0; i < _n; i++)
                _v[i][0] = -_v[i][0];
        if ((sf <= (0.9 * _d[0])) || ((0.9 * sf) >= _d[0]))
            for (i=1; i<_n; i++)
                _d[i] = 0.0;
        for (k=1; k<_n; k++) 
            {
            for (i=0; i<_n; i++)
                y[i] = _minimum_point[i];
            sf = fx;
            illc = illc || (kt > 0);
            next:
               kl = k;
               df = 0.0;
               if (illc) 
                    {        
                    // random step to get off resolution valley
                    for (i=0; i<_n; i++) 
                        {
                        z[i] = (0.1 * ldt + t2 * pow(10.0,(double)kt)) * (random() - 0.5);
                        s = z[i];
                        for (j=0; j < _n; j++)
                            _minimum_point[j] += s * _v[j][i];
                        }
                    fx = (*_func)(_minimum_point, _n);
                    nf++;
                    }
                    
                // minimize along non-conjugate directions
                for (k2=k; k2<_n; k2++) 
                    {  
                    sl = fx;
                    s = 0.0;
                    min(k2, 2, &_d[k2], &s, fx, 0);
                    if (illc)
                        {
                        double szk = s + z[k2];
                        s = _d[k2] * szk*szk;
                        }
                    else 
                        s = sl - fx;
                    if (df < s)
                        {
                        df = s;
                        kl = k2;
                        }
                    }
                if (!illc && (df < fabs(100.0 * macheps * fx))) 
                    {
                    illc = 1;
                    goto next;
                    }
                if ((k == 1) && (prin > 1))
                    vecprint("\n... New Direction ...",_d,_n);
                    
                // minimize along conjugate directions
                for (k2=0; k2<=k-1; k2++) 
                    {
                    s = 0.0;
                    min(k2, 2, &_d[k2], &s, fx, 0);
                    }
                f1 = fx;
                fx = sf;
                lds = 0.0;
                for (i=0; i<_n; i++) 
                    {
                    sl = _minimum_point[i];
                    _minimum_point[i] = y[i];
                    y[i] = sl - y[i];
                    sl = y[i];
                    lds = lds + sl*sl;
                    }
                lds = sqrt(lds);
                if (lds > small) 
                    {
                    for (i=kl-1; i>=k; i--) 
                        {
                        for (j=0; j < _n; j++)
                            _v[j][i+1] = _v[j][i];
                        _d[i+1] = _d[i];
                        }
                    _d[k] = 0.0;
                    for (i=0; i < _n; i++)
                        _v[i][k] = y[i] / lds;
                    min(k, 4, &_d[k], &lds, f1, 1);
                    if (lds <= 0.0) 
                        {
                        lds = -lds;
                        for (i=0; i<_n; i++)
                            _v[i][k] = -_v[i][k];
                        }
                    }
                ldt = ldfac * ldt;
                if (ldt < lds)
                    ldt = lds;
                if (prin > 1)
                    print();
                t2 = 0.0;
                for (i=0; i<_n; i++)
                    t2 += _minimum_point[i]*_minimum_point[i];
                t2 = m2 * sqrt(t2) + t;
                if (ldt > (0.5 * t2))
                    kt = 0;
                else 
                    kt++;
                if (kt > ktm)
                    goto fret; 
                }
                
        //  try quadratic extrapolation in case
        //  we are stuck in a curved valley 
        quad();
        dn = 0.0;
        for (i=0; i<_n; i++) 
            {
            _d[i] = 1.0 / sqrt(_d[i]);
            if (dn < _d[i])
                dn = _d[i];
            }
       if (prin > 2)
          matprint("\n... New Matrix of Directions ...",_v,_n);
            for (j=0; j<_n; j++) 
                {
                s = _d[j] / dn;
                for (i=0; i < _n; i++)
                    _v[i][j] *= s;
                }
        if (scbd > 1.0) 
            {       
            /* scale axis to reduce condition number */
            s = vlarge;
            for (i=0; i<_n; i++) 
                {
                sl = 0.0;
                for (j=0; j < _n; j++)
                sl += _v[i][j]*_v[i][j];
                z[i] = sqrt(sl);
                if (z[i] < m4)
                z[i] = m4;
                if (s > z[i])
                s = z[i];
                }
            for (i=0; i<_n; i++) 
                {
                sl = s / z[i];
                z[i] = 1.0 / sl;
                if (z[i] > scbd) 
                    {
                    sl = 1.0 / scbd;
                    z[i] = scbd;
                    }
                }
            }
        for (i=1; i<_n; i++)
            for (j=0; j<=i-1; j++) 
                {
                s = _v[i][j];
                _v[i][j] = _v[j][i];
                _v[j][i] = s;
                }
        minfit(_n, macheps, vsmall, _v, _d);
        if (scbd > 1.0) 
           {
            for (i=0; i<_n; i++) 
                {
                s = z[i];
                for (j=0; j<_n; j++)
                    _v[i][j] *= s;
                }
            for (i=0; i<_n; i++) 
                {
                s = 0.0;
                for (j=0; j<_n; j++)
                    s += _v[j][i]*_v[j][i];
                s = sqrt(s);
                _d[i] *= s;
                s = 1.0 / s;
                for (j=0; j<_n; j++)
                    _v[j][i] *= s;
                }
           }
       for (i=0; i<_n; i++) 
           {
           if ((dn * _d[i]) > large)
              _d[i] = vsmall;
           else if ((dn * _d[i]) < small)
              _d[i] = vlarge;
           else 
              _d[i] = pow(dn * _d[i],-2.0);
           }
       Sort();               /* sort the new eigenvalues and eigenvectors */
       dmin = _d[_n-1];
       if (dmin < small)
          dmin = small;
       illc = (m2 * _d[0]) > dmin;
       if ((prin > 2) && (scbd > 1.0))
          vecprint("\n... Scale Factors ...",z,_n);
       if (prin > 2)
          vecprint("\n... Eigenvalues of A ...",_d,_n);
       if (prin > 2)
          matprint("\n... Eigenvectors of A ...",_v,_n);

       if ((maxfun > 0) && (nl > maxfun)) 
           {
            if (prin)
                printf("\n... maximum number of function calls reached ...\n");
            goto fret;
           }
       }    //main loop

fret:
   if (prin > 0) {
         vecprint("\n... Final solution is ...", _minimum_point, _n);
         printf("\n... ChiSq reduced to %20.10e ...\n", fx);
	 printf("... after %20u function calls.\n", nf);
   }
   
   return(fx);
    }
    
#endif
