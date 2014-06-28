from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class RelativeRateDistribution(RelRateDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the multivariate relative rate probability distribution.
    *** Not finished documenting this class yet ***

    """
    def __init__(self, params, coeffs):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the parameters and coefficients (probabilities of rate 
        categories) of the RelativeRateDistribution object as a tuple. e.g.,

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((10,20,30,40),(0.25,0.25,0.25,0.25))
        >>> print '(%.1f, %.1f, %.1f, %.1f)' % tuple([x for x in d.getMean()])
        (0.4, 0.8, 1.2, 1.6)
        """
        RelRateDistBase.__init__(self, params, coeffs)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this RelativeRateDistribution.
        
        """
        return RelRateDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the relative rate distribution is 
        continuous.

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,1,1,1))
        >>> print d.isDiscrete()
        False
        
        """
        return RelRateDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'RelativeRateDistribution'

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,1,1,1))
        >>> print d.getDistName()
        RelativeRateDistribution
        
        """
        return RelRateDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another 
        RelativeRateDistribution object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> print d.__str__()
        RelativeRateDistribution((1.00000,2.00000,3.00000,4.00000))
        
        """
        return RelRateDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another 
        RelativeRateDistribution object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> print d.__repr__()
        RelativeRateDistribution((1.00000,2.00000,3.00000,4.00000))
        >>> print d
        RelativeRateDistribution((1.00000,2.00000,3.00000,4.00000))
        
        """
        return RelRateDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d = RelativeRateDistribution((1.0, 1.0, 1.0, 1.0))
        >>> d.setLot(g)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.026209524600
        1.685719371974
        1.081826860846
        1.206244242581
        
        >>> for x in d.sample():
        ...     print "%.12f" % x
        2.410190037018
        0.098035815791
        1.313892681879
        0.177881465312
        
        >>> g.setSeed(1357)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.026209524600
        1.685719371974
        1.081826860846
        1.206244242581
            
        In this example, only one random number generator (g) was involved;
        however, one could pass g to several different probability distribu-
        tions, thus ensuring that the entire sequence of random numbers could
        be recreated by keeping track of only one seed value. If setLot is not
        used, each distribution object maintains its own random number
        generator that is initialized using the system clock at the time the
        object is created, making it difficult to replicate results.

        """
        return RelRateDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect.

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> d.setSeed(135)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.000418520550
        1.082762115186
        0.149006532930
        2.767812831334
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.936279048974
        0.683180418928
        0.765497579702
        1.615042952397
        >>> d.setSeed(135)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.000418520550
        1.082762115186
        0.149006532930
        2.767812831334
    
        """
        return RelRateDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return RelRateDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution). Because this is a multivariate distribution,
        the object returned is a tuple.

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> print '(%.1f, %.1f, %.1f, %.1f)' % tuple([x for x in d.getMean()])
        (0.4, 0.8, 1.2, 1.6)

        """
        return RelRateDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution). Because this is a multivariate
        distribution, the object returned is a tuple.

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> print '(%.5f, %.5f, %.5f, %.5f)' % tuple([x for x in d.getVar()])
        (0.13091, 0.23273, 0.30545, 0.34909)

        """
        return RelRateDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution). Because this is
        a multivariate distribution, the object returned is a tuple.

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> print '(%.5f, %.5f, %.5f, %.5f)' % tuple([x for x in d.getStdDev()])
        (0.36181, 0.48242, 0.55268, 0.59084)

        """
        return RelRateDistBase.getStdDev(self)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the relative rate distribution 
        specified by this RelativeRateDistribution object. Python list 
        comprehensions can be used to store many simulated samples for use in 
        subsequent calculations.

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3,4))
        >>> d.setSeed(97531)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        0.483550524056
        0.054914097328
        2.049601112088
        1.411934266528
            
        """
        return RelRateDistBase.sample(self)
        
    def approxCDF(self, x, n = 10000):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Approximates the cumulative distribution function evaluated at the
        supplied point x. The precision of the approximation is controlled
        by nsamples. The approximation is done using a brute force approach:
        nsamples samples are drawn from this RelativeRateDistribution, and the
        proportion of those samples that are inside the region defined by x
        is returned as the approximated CDF. The supplied point x should be
        a tuple of length k, where k is one fewer than the number of
        parameters of the RelativeRateDistribution. If x has length greater
        than k, the extra elements will be ignored. The setSeed call is needed
        in the following example to ensure that the approximated CDF will be 
        the same every time this example is run:

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,1))
        >>> d.setSeed(1357)
        >>> print d.approxCDF(d.getMean())
        0.5059
         
        """
        assert n > 0, 'n must be greater than zero in RelativeRateDistribution.approxCDF function'
        nparams = RelRateDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to approxCDF has length %d but length expected was %d' % (len(x), nparams)
        return RelRateDistBase.approxCDF(self, x, n)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,2,3))
        >>> d.setCoefficients((1.0/6.0,1.0/3.0,1.0/2.0))
        >>> print '%.12f' % d.getLnPDF((0.5, 1.0))
        -0.972632485808

        The 3-parameter RelativeRateDistribution(a,b,c) density defined in the
        exampel above is

              (x1*p1)^(a-1) (x2*p2)^(b-1) (x3*p3)^(c-1)
        p1 p2 -----------------------------------------
              Gamma(a) Gamma(b) Gamma(c) / Gamma(a+b+c)
        
        where a = 1, b = 2, c = 3, p1 = 1/6, p2 = 1/3, p3 = 1/2 and
        Gamma is the Gamma function, not the Gamma probability distribution.
        Note that the argument x is a tuple, which in this example could be
        supplied as either

        x = (0.5, 1.0)
        
        or
        
        x = (0.5, 1.0, 7/6)
        
        (The final value x3 is determined by x1, x2, p1, p2 and p3 and thus
        does not need to be supplied.)
                
        """
        nparams = RelRateDistBase.getNParams(self)
        assert len(x) >= nparams - 1, 'Vector supplied to getLnPDF has length %d but should be greater than or equal to %d' % (len(x), nparams - 1)
        return RelRateDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function simply calls getLnPDF.
        
        """
        nparams = RelRateDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to getRelativeLnPDF has length %d but length expected was %d' % (len(x), nparams)
        return RelRateDistBase.getRelativeLnPDF(self, x)
        
    #def setCoefficients(self, p):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Allows user to provide coefficients for this distribution, overriding
    #    the default coefficients, which are all 1/dim, where dim is the 
    #    dimension (number of parameters).
    #    
    #    """
    #    nparams = RelRateDistBase.getNParams(self)
    #    assert nparams == len(p), 'Vector supplied to setCoefficients has length %d but length expected was %d' % (len(x), nparams)
    #    return RelRateDistBase.setCoefficients(self, p)
        
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
    #    >>> from phycas.probdist import *
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
    #    nparams = RelRateDistBase.getNParams(self)
    #    assert len(mean.shape) == 1, 'Mean vector supplied to Dirichlet.setMeanAndVariance should be a single-dimensional array'
    #    assert len(variance.shape) == 1, 'Variance vector supplied to Dirichlet.setMeanAndVariance should be a single-dimensional array'
    #    assert nparams == len(mean), 'Mean vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(mean))
    #    assert nparams == len(variance), 'Variance vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(variance))
    #    return RelRateDistBase.setMeanAndVariance(self, mean, variance)
        
    # Comment out this version if numarray is re-introduced
    def setMeanAndVariance(self, mean, variance):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the (multivariate) mean and variance of this distribution.
        Note that the variances are given as a one-dimensional tuple, leaving
        out the covariances (which are not needed to fully specify the
        RelativeRateDistribution distribution). For example:
        
        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1.0,1.0,1.0))
        >>> m = (1.2,1.0,0.8)
        >>> print '(%.1f, %.1f, %.1f)' % m
        (1.2, 1.0, 0.8)
        >>> v = (0.069677419,0.064516129,0.056774194)
        >>> print '(%.6f, %.6f, %.6f)' % v
        (0.069677, 0.064516, 0.056774)
        >>> d.setMeanAndVariance(m,v)
        >>> print '(%.1f, %.1f, %.1f)' % tuple(d.getMean())
        (1.2, 1.0, 0.8)
        >>> print '(%.6f, %.6f, %.6f)' % tuple(d.getVar())
        (0.069677, 0.064516, 0.056774)
        
        """
        nparams = RelRateDistBase.getNParams(self)
        assert nparams == len(mean), 'Mean vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(mean))
        assert nparams == len(variance), 'Variance vector supplied to setMeanAndVariance should have %d elements, but %d were found' % (nparams, len(variance))
        return RelRateDistBase.setMeanAndVariance(self, mean, variance)
        
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
    #    >>> from phycas.probdist import *
    #    >>> from numpy.numarray import array
    #    >>> d = Dirichlet((1,1,1))
    #    >>> print d.getVarCovarMatrix()
    #    [[ 0.05555556 -0.02777778 -0.02777778]
    #     [-0.02777778  0.05555556 -0.02777778]
    #     [-0.02777778 -0.02777778  0.05555556]]
    #    
    #    """
    #    return RelRateDistBase.getVarCovarMatrix(self)

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
        sum of the n RelativeRateDistribution parameters, c_i be the ith. 
        parameter, and denom be c*c*(c + 1), then

        Var_i = n*c_i*(c - c_i)/denom
        Cov_ij = -n*c_i*c_j/denom
        
        In this example,

        n      = 3
        c_1    = 1 
        c_2    = 1 
        c_3    = 1
        c      = 3
        denom  = 3*3*4 = 36
        Var_i  = 3*1*2/36 = 0.00091827
        Cov_ij = -3*1*1/36 = -0.083333
        
        >>> from phycas.probdist import *
        >>> d = RelativeRateDistribution((1,1,1))
        >>> d.printSquareMatrix(d.getVarCovarMatrix())
         0.50000000  -0.25000000  -0.25000000 
        -0.25000000   0.50000000  -0.25000000 
        -0.25000000  -0.25000000   0.50000000 
            
        """
        return RelRateDistBase.getVarCovarMatrix(self)
        
    def getNParams(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of parameters in the Dirichlet distribution. 
        For example:
        
        >>> from phycas.probdist import *
        >>> d1 = RelativeRateDistribution((1,2,3,4))
        >>> print d1.getNParams()
        4
        >>> d2 = RelativeRateDistribution((1,1))
        >>> print d2.getNParams()
        2
        
        """
        return RelRateDistBase.getNParams(self)
    
