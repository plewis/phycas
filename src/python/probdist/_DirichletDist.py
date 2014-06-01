from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Dirichlet(DirichletDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the multivariate Dirichlet probability distribution.
    *** Not finished documenting this class yet ***

    """
    def __init__(self, c):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the parameters of the Dirichlet object as a tuple. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,1,1,1))
        >>> print d.getMean()
        (0.25, 0.25, 0.25, 0.25)

        """
        DirichletDistBase.__init__(self, c)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Dirichlet distribution.
        
        """
        return DirichletDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the Dirichlet distribution is continuous.

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,1,1,1))
        >>> print d.isDiscrete()
        False
        
        """
        return DirichletDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Dirichlet'

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,1,1,1))
        >>> print d.getDistName()
        Dirichlet
        
        """
        return DirichletDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Dirichlet
        object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.__str__()
        Dirichlet((1.00000, 2.00000, 3.00000, 4.00000))
        
        """
        return DirichletDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Dirichlet
        object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.__repr__()
        Dirichlet((1.00000, 2.00000, 3.00000, 4.00000))
        >>> print d
        Dirichlet((1.00000, 2.00000, 3.00000, 4.00000))
        
        """
        return DirichletDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d = Dirichlet((1.0, 1.0, 1.0, 1.0))
        >>> d.setLot(g)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.006552381150
        0.421429842993
        0.270456715211
        0.301561060645
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.602547509254
        0.024508953948
        0.328473170470
        0.044470366328
        >>> g.setSeed(1357)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.006552381150
        0.421429842993
        0.270456715211
        0.301561060645
        
        In this example, only one random number generator (g) was involved;
        however, one could pass g to several different probability distribu-
        tions, thus ensuring that the entire sequence of random numbers could
        be recreated by keeping track of only one seed value. If setLot is not
        used, each distribution object maintains its own random number
        generator that is initialized using the system clock at the time the
        object is created, making it difficult to replicate results.

        """
        return DirichletDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect.

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> d.setSeed(135)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.000104630137
        0.270690528796
        0.037251633232
        0.691953207834
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.234069762243
        0.170795104732
        0.191374394925
        0.403760738099
        >>> d.setSeed(135)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.000104630137
        0.270690528796
        0.037251633232
        0.691953207834

        """
        return DirichletDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return DirichletDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution). Because this is a multivariate distribution,
        the object returned is a tuple.

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.getMean()
        (0.10000000000000001, 0.20000000000000001, 0.29999999999999999, 0.40000000000000002)

        """
        return DirichletDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution). Because this is a multivariate
        distribution, the object returned is a tuple.

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.getVar()
        (0.0081818181818181825, 0.014545454545454545, 0.019090909090909092, 0.02181818181818182)

        """
        return DirichletDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution). Because this is
        a multivariate distribution, the object returned is a tuple.

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.getStdDev()
        (0.090453403373329092, 0.12060453783110545, 0.13816985594155148, 0.14770978917519928)

        """
        return DirichletDistBase.getStdDev(self)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the Dirichlet distribution specified
        by this Dirichlet object. Python list comprehensions can be used
        to store many simulated samples for use in subsequent calculations.

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> d.setSeed(97531)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.120887631014
        0.013728524332
        0.512400278022
        0.352983566632
        
        """
        return DirichletDistBase.sample(self)
        
    def approxCDF(self, x, n = 10000):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Approximates the cumulative distribution function evaluated at the
        supplied point x. The precision of the approximation is controlled
        by nsamples. The approximation is done using a brute force approach:
        nsamples samples are drawn from this Dirichlet distribution, and the
        proportion of those samples that are inside the region defined by x
        is returned as the approximated CDF. The supplied point x should be
        a tuple of length k, where k is one fewer than the number of
        parameters of the Dirichlet distribution. If x has length greater
        than k, the extra elements will be ignored. In the following example,
        the result returned from approxCDF for a Dirichlet((1,1))
        distribution is compared to the exact result for the equivalent
        univariate Beta(1,1) distribution (the setSeed call is needed
        to ensure that the approximated CDF will be the same every time this
        example is run):

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,1))
        >>> d.setSeed(1357)
        >>> print d.approxCDF(d.getMean())
        0.5059
        >>> b = Beta(1,1)
        >>> print b.getCDF(b.getMean())
        0.5
         
        """
        assert n > 0, 'n must be greater than zero in Dirichlet.approxCDF function'
        nparams = DirichletDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to approxCDF has length %d but length expected was %d' % (len(x), nparams)
        return DirichletDistBase.approxCDF(self, x, n)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.getLnPDF((0.4,0.3,0.2,0.1))
        -1.01368307788

        For example, the 4-parameter Dirichlet density is

        Gamma(a+b+c+d) p^(a-1) q^(b-1) r^(c-1) (1-p-q-r)^(d-1)
        ------------------------------------------------------
                Gamma(a) Gamma(b) Gamma(c) Gamma(d)
        
        where a = 1, b = 2, c = 3, d = 4, p = 0.4, q = 0.3, r = 0.2 and
        Gamma is the Gamma function, not the Gamma probability distribution.
        Note that the argument x is a tuple, which in this example would be

        x = (p, q, r, 1-p-q-r) = (0.4, 0.3, 0.2, 0.1)
        
        The natural logarithm of the relative density is thus

        (a-1) log(p) + (b-1) log(q) + (c-1) log(r) + (d-1) log(1-p-q-r)
          - log(Gamma(a)) - log(Gamma(b)) - log(Gamma(c)) - log(Gamma(d))
          + log(Gamma(a+b+c+d)) 
        
        For the example given, this equals

        (1-1) log(0.4) + (2-1) log(0.3) + (3-1) log(0.2) + (4-1) log(0.1)
          - log(Gamma(1)) - log(Gamma(2)) - log(Gamma(3)) - log(Gamma(4))
          + log(Gamma(1+2+3+4))
        = 0.0 + log(0.3) + 2 log(0.2) + 3 log(0.1)
          - 0.0 - 0.0 - log(2!) - log(3!)
          + log(9!)
        = -1.01368307788
        
        """
        nparams = DirichletDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to getLnPDF has length %d but length expected was %d' % (len(x), nparams)
        return DirichletDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the relative probability density function at the supplied
        value x. Returns the natural logarithm of the relative density at x.
        Use this function if speed is important but normalization is not, say
        in MCMC calculations. Use getLnPDF instead if you need to have a
        correctly normalized density value (i.e. from a density function that
        integrates to 1.0)

        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,2,3,4))
        >>> print d.getRelativeLnPDF((0.4,0.3,0.2,0.1))
        -11.3306039082

        For example, the 4-parameter Dirichlet density is

        Gamma(a+b+c+d) p^(a-1) q^(b-1) r^(c-1) (1-p-q-r)^(d-1)
        ------------------------------------------------------
                Gamma(a) Gamma(b) Gamma(c) Gamma(d)
        
        where a = 1, b = 2, c = 3, d = 4, p = 0.4, q = 0.3, r = 0.2 and
        Gamma is the Gamma function, not the Gamma probability distribution.
        The relative density requires only the four terms containing p, q,
        and r in the numerator, so the natural logarithm of the relative
        density is

        (a-1) log(p) + (b-1) log(q) + (c-1) log(r) + (d-1) log(1-p-q-r)
        
        For the example given, this equals
        (1-1) log(0.4) + (2-1) log(0.3) + (3-1) log(0.2) + (4-1) log(0.1)
        = log(0.3) + 2 log(0.2) + 3 log(0.1)
        = -11.3306039082
        
        """
        nparams = DirichletDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to getRelativeLnPDF has length %d but length expected was %d' % (len(x), nparams)
        return DirichletDistBase.getRelativeLnPDF(self, x)
        
    # Uncomment this version if numarray is re-introduced
    #def setMeanAndVariance(self, mean, variance):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Sets the (multivariate) mean and variance of this distribution.
    #    Note: mean and variance are numarray array objects rather than simple
    #    tuples. This was an exercise to see if numarray could be used for this
    #    purpose, but tuples would be better. Note that the variances are
    #    given as a one-dimensional array, leaving out the covariances (which
    #    are not needed to fully specify the Dirichlet distribution).
    #    For example:
    #
    #    >>> from phycas.ProbDist import *
    #    >>> from numarray import array
    #    >>> d = Dirichlet((1,1,1))
    #    >>> m = array([1./9., 3./9., 5./9.])
    #    >>> print m
    #    [ 0.11111111  0.33333333  0.55555556]
    #    >>> v = array([8./810.,18./810.,20./810.])
    #    >>> print v
    #    [ 0.00987654  0.02222222  0.02469136]
    #    >>> d.setMeanAndVariance(m,v)
    #    >>> d.getMean()
    #    (0.1111111111111111, 0.33333333333333331, 0.55555555555555558)
    #    >>> d.getVar()
    #    (0.009876543209876543, 0.022222222222222223, 0.024691358024691357)
    #    
    #    """
    #    nparams = DirichletDistBase.getNParams(self)
    #    assert len(mean.shape) == 1, 'Mean vector supplied to Dirichlet.setMeanAndVariance should be a single-dimensional array'
    #    assert len(variance.shape) == 1, 'Variance vector supplied to Dirichlet.setMeanAndVariance should be a single-dimensional array'
    #    assert nparams == len(mean), 'Mean vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(mean))
    #    assert nparams == len(variance), 'Variance vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(variance))
    #    return DirichletDistBase.setMeanAndVariance(self, mean, variance)
        
    # Comment out this version if numarray is re-introduced
    def setMeanAndVariance(self, mean, variance):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the (multivariate) mean and variance of this distribution.
        Note that the variances are given as a one-dimensional tuple, leaving
        out the covariances (which are not needed to fully specify the
        Dirichlet distribution). For example:
        
        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,1,1))
        >>> m = (1./9., 3./9., 5./9.)
        >>> print m
        (0.1111111111111111, 0.33333333333333331, 0.55555555555555558)
        >>> v = (8./810.,18./810.,20./810.)
        >>> print v
        (0.009876543209876543, 0.022222222222222223, 0.024691358024691357)
        >>> d.setMeanAndVariance(m,v)
        >>> d.getMean()
        (0.1111111111111111, 0.33333333333333331, 0.55555555555555558)
        >>> d.getVar()
        (0.009876543209876543, 0.022222222222222223, 0.024691358024691357)
        
        """
        nparams = DirichletDistBase.getNParams(self)
        assert nparams == len(mean), 'Mean vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(mean))
        assert nparams == len(variance), 'Variance vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(variance))
        return DirichletDistBase.setMeanAndVariance(self, mean, variance)
        
    # Uncomment this version if numarray is re-introduced
    #def getVarCovarMatrix(self):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Returns the variance-covariance matrix of this distribution. The
    #    matrix is returned in the form of a numarray object. Letting c be the
    #    sum of the n Dirichlet parameters, c_i be the ith. parameter, and
    #    denom be c*c*(c + 1), then
    #
    #    Var_i = c_i*(c - c_i)/denom
    #    Cov_ij = -c_i*c_j/denom
    #
    #    In this example,
    #
    #    n      = 3
    #    c_1    = 1 
    #    c_2    = 1 
    #    c_3    = 1
    #    c      = 3
    #    denom  = 3*3*4 = 36
    #    Var_i  = 1*2/36 = 0.05555556
    #    Cov_ij = -1*1/36 = -0.02777778
    #
    #    >>> from phycas.ProbDist import *
    #    >>> from numpy.numarray import array
    #    >>> d = Dirichlet((1,1,1))
    #    >>> print d.getVarCovarMatrix()
    #    [[ 0.05555556 -0.02777778 -0.02777778]
    #     [-0.02777778  0.05555556 -0.02777778]
    #     [-0.02777778 -0.02777778  0.05555556]]
    #    
    #    """
    #    return DirichletDistBase.getVarCovarMatrix(self)

    # Comment out this version if numarray is re-introduced
    def printSquareMatrix(self, matrix):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Utility function used to interpret a vector as a square matrix and
        print it out.
        
        """
        total_len = len(matrix)
        import math
        row_len = int(math.sqrt(total_len))
        assert row_len*row_len == total_len, 'Attempting to print a matrix that is not square'
        k = 0
        for i in range(row_len):
            for j in range(row_len):
                print '% .8f ' % matrix[k],
                k += 1
            print
        
    # Comment out this version if numarray is re-introduced
    def getVarCovarMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance-covariance matrix of this distribution. The
        matrix is returned in the form of a numarray object. Letting c be the
        sum of the n Dirichlet parameters, c_i be the ith. parameter, and
        denom be c*c*(c + 1), then

        Var_i = c_i*(c - c_i)/denom
        Cov_ij = -c_i*c_j/denom
        
        In this example,

        n      = 3
        c_1    = 1 
        c_2    = 1 
        c_3    = 1
        c      = 3
        denom  = 3*3*4 = 36
        Var_i  = 1*2/36 = 0.05555556
        Cov_ij = -1*1/36 = -0.02777778
        
        >>> from phycas.ProbDist import *
        >>> d = Dirichlet((1,1,1))
        >>> d.printSquareMatrix(d.getVarCovarMatrix())
         0.05555556  -0.02777778  -0.02777778 
        -0.02777778   0.05555556  -0.02777778 
        -0.02777778  -0.02777778   0.05555556 
        
        """
        return DirichletDistBase.getVarCovarMatrix(self)
        
    def getNParams(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of parameters in the Dirichlet distribution. 
        For example:
        
        >>> from phycas.ProbDist import *
        >>> d1 = Dirichlet((1,2,3,4))
        >>> print d1.getNParams()
        4
        >>> d2 = Dirichlet((1,1))
        >>> print d2.getNParams()
        2
        
        """
        return DirichletDistBase.getNParams(self)
    
