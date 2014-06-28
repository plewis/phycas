from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *
import math

class MVNormal(MVNormalDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the multivariate normal probability distribution.
    *** Not finished documenting this class yet ***

    """
    def __init__(self, m = None, v = None):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        If called with no arguments, constructs a default object representing
        a bivariate standard normal distribution. e.g.,

        >>> from phycas.probdist import *
        >>> d = MVNormal()
        >>> print d.getDistName()
        Bivariate Normal
        
        Alternatively, the constructor can be provided with two arguments to
        create a multivariate normal distribution of arbitrary dimension.
        The first argument should be a list or tuple representing the mean
        vector. The second argument should be a two dimensional list or tuple
        representing the variance/covariance matrix. The variance/covariance
        matrix should be square, symmetric and have a dimension that matches
        the mean vector. e.g.,

        >>> from phycas.probdist import *
        >>> m = (1,2,3)
        >>> v = ((2,1,1),(1,2,1),(1,1,2))
        >>> d = MVNormal(m, v)
        >>> print d.getDistName()
        Multivariate Normal

        """
        MVNormalDistBase.__init__(self)
        dim = 0
        if m is None:
            m = (0.0, 0.0)
            dim = 2
        else:
            dim = len(m)
            if dim < 1:
                raise RuntimeError('cannot create an MVNormal distribution with a dimension less than 1')
        if v is None:
            v = self._identityMatrix(dim)
        else:
            self._checkVarCovMatrix(v, dim)
        self.setMeanAndVariance(m, self._flatten(v))
        
    def _flatten(self, m):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Flattens a two-dimensional matrix m rowwise into a one-dimensional
        tuple. 
        
        """        
        nrows = len(m)
        ncols = len(m[0])
        flat = []
        for i in range(nrows):
            for j in range(ncols):
                flat.append(m[i][j])
        return tuple(flat)
        
    def _checkVarCovMatrix(self, v, dim):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs sanity checks on a variance/covariance matrix supplied by
        the user to ensure that it is square, of the correct dimensions (as
        indicated by dim), and symmetric. If it violates any of these 
        requirements, a RuntimeError exception will be raised.
        
        """        
        if len(v) != dim:
            raise RuntimeError('Number of rows of the supplied variance/covariance matrix (%d) differs from expected (%d)' % (len(v), dim))
        for i in range(dim):
            if len(v[i]) != dim:
                raise RuntimeError('Number of columns of variance/covariance matrix (%d) differs from expected (%d) in row %d' % (len(v[i]), dim, i+1))
            for j in range(dim):
                if math.fabs(v[i][j] - v[j][i]) > 0.00001:
                    raise RuntimeError('Supplied variance/covariance matrix is not symmetrical for row %d, column %d' % (i+1, j+1))
        
    def _identityMatrix(self, dim):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates an identity matrix with dim rows and dim columns.
        
        """
        I = [tuple([i==j and 1.0 or 0.0 for j in range(dim)]) for i in range(dim)]
        return tuple(I)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this MVNormal distribution.
        
        """
        return MVNormalDistBase.clone(self)
        
    def fit(self, data):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fits this MVNormal distribution to the supplied data, which should be
        in the form of a two-dimensional tuple with rows equal to observations
        and columns equal to variables. Returns True if successful and False
        if there were problems (e.g. the variance-covariance matrix could not
        be inverted).
        
        """
        nrows = len(data)
        ncols = len(data[0])
        return MVNormalDistBase.fit(self, nrows, ncols, self._flatten(data))
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the MVNormal distribution is continuous.

        """
        return MVNormalDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Bivariate Normal' if the number of dimensions is
        2, otherwise returns 'Multivariate Normal'

        >>> from phycas.probdist import *
        >>> m = (0,0)
        >>> v = ((1,0),(0,1))
        >>> d = MVNormal(m, v)
        >>> print d.getDistName()
        Bivariate Normal
        >>> m = (0,0,0)
        >>> v = ((1,0,0),(0,1,0),(0,0,1))
        >>> d = MVNormal(m, v)
        >>> print d.getDistName()
        Multivariate Normal
        
        """
        return MVNormalDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another MVNorm
        object identical to this one. Same as __repr__() function.

        >>> from phycas.probdist import *
        >>> m = (1,2,3)
        >>> v = ((1,0,0),(0,1,0),(0,0,1))
        >>> d = MVNormal(m,v)
        >>> print d.__str__()
        mu = (
          1.00000,
          2.00000,
          3.00000
          )
        Sigma = (
          (1.00000, 0.00000, 0.00000),
          (0.00000, 1.00000, 0.00000),
          (0.00000, 0.00000, 1.00000)
          )
        MVNormal(mu, Sigma)
    
        """
        return MVNormalDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another MVNorm
        object identical to this one. Same as __str__() function.

        >>> from phycas.probdist import *
        >>> m = (1,2,3)
        >>> v = ((1,0,0),(0,1,0),(0,0,1))
        >>> d = MVNormal(m,v)
        >>> print d.__str__()
        mu = (
          1.00000,
          2.00000,
          3.00000
          )
        Sigma = (
          (1.00000, 0.00000, 0.00000),
          (0.00000, 1.00000, 0.00000),
          (0.00000, 0.00000, 1.00000)
          )
        MVNormal(mu, Sigma)
    
        """
        return MVNormalDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. In the following example, two samples are drawn from a
        MVNormal distribution before the seed is reset to its initial value.
        Note that the next sample drawn is identical to the first.

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d = MVNormal()
        >>> d.setLot(g)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        -2.303677281009
        -0.008075169733
        >>> for x in d.sample():
        ...     print "%.12f" % x
        -0.368038623909
        -0.283937813839
        >>> g.setSeed(1357)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        -2.303677281009
        -0.008075169733
        
        In this example, only one random number generator (g) was involved;
        however, one could pass g to several different probability distribu-
        tions, thus ensuring that the entire sequence of random numbers could
        be recreated by keeping track of only one seed value. If setLot is not
        used, each distribution object maintains its own random number
        generator that is initialized using the system clock at the time the
        object is created, making it difficult to replicate results.

        """
        return MVNormalDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect instead.

        >>> from phycas.probdist import *
        >>> d = MVNormal((0,0,0),((1,0,0),(0,1,0),(0,0,1)))
        >>> d.setSeed(135)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        -3.073855199517
        0.698607737112
        -2.471972303472
        >>> for x in d.sample():
        ...     print "%.12f" % x
        1.390212404902
        1.429551075843
        0.149676695051
        >>> d.setSeed(135)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        -3.073855199517
        0.698607737112
        -2.471972303472
    
        """
        return MVNormalDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return MVNormalDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. Because this is a multivariate 
        distribution, the object returned is a tuple.

        >>> from phycas.probdist import *
        >>> d = MVNormal(m=(1,2,3),v=((1,0,0),(0,1,0),(0,0,1)))
        >>> for x in d.getMean():
        ...     print '%.5f ' % x
        1.00000 
        2.00000 
        3.00000 
         
        """
        return MVNormalDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance/covariance matrix of the distribution as a two-
        dimensional tuple.

        >>> from phycas.probdist import *
        >>> d = MVNormal(m=(1,2,3),v=((1,2,3),(2,3,4),(3,4,5)))
        >>> print d.getVar()
        ((1.0, 2.0, 3.0), (2.0, 3.0, 4.0), (3.0, 4.0, 5.0))

        """
        varcov = MVNormalDistBase.getVar(self)
        dimsq = len(varcov)
        dim = int(math.sqrt(dimsq))
        k = 0
        v = []
        for i in range(dim):
            tmp = []
            for j in range(dim):
                tmp.append(varcov[k])
                k += 1
            v.append(tuple(tmp))
        return tuple(v) 
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the square root of the variance/covariance matrix, where the
        square root is applied to each element of the variance/covariance
        matrix. The result is returned as a two-dimensional tuple.

        >>> from phycas.probdist import *
        >>> d = MVNormal(m=(1,2,3),v=((1,2,3),(2,3,4),(3,4,5)))
        >>> sd = d.getStdDev()
        >>> for i in range(3):
        ...     for j in range(3):
        ...         print '%4.2f' % sd[i][j],
        ...     print
        1.00 1.41 1.73
        1.41 1.73 2.00
        1.73 2.00 2.24
        
        """
        sd = MVNormalDistBase.getStdDev(self)
        dimsq = len(sd)
        dim = int(math.sqrt(dimsq))
        k = 0
        v = []
        for i in range(dim):
            tmp = []
            for j in range(dim):
                tmp.append(sd[k])
                k += 1
            v.append(tuple(tmp))
        return tuple(v) 
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the multivariate normal distribution
        specified by this MVNormal object.

        >>> from phycas.probdist import *
        >>> d = MVNormal(m=(10,20,30,40),v=((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)))
        >>> d.setSeed(97531)
        >>> for x in d.sample():
        ...     print "%.12f" % x
        10.717002157630
        17.743363252531
        31.577667572968
        40.268951828105
    
        """
        return MVNormalDistBase.sample(self)
        
    def approxCDF(self, x, n = 10000):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Approximates the cumulative distribution function evaluated at the
        supplied point x. The precision of the approximation is controlled
        by n. The approximation is done using a brute force approach: n 
        samples are drawn from this MVNormal distribution, and the proportion
        of those samples that are less than x is returned as the approximated
        CDF. The supplied point x should be a tuple of length k, where k is
        equal to the number of parameters of the MVNormal distribution. If x
        has length different than k, a RuntimeException will be raised. 

        >>> from phycas.probdist import *
        >>> d = MVNormal()
        >>> d.setSeed(1357)
        >>> print d.approxCDF((0.0,0.0))
        0.2536

        """
        assert n > 0, 'n must be greater than zero in MVNormal.approxCDF function'
        nparams = MVNormalDistBase.getNParams(self)
        if nparams != len(x):
            raise RuntimeException('Vector supplied to approxCDF has length %d but length expected was %d' % (len(x), nparams))
        return MVNormalDistBase.approxCDF(self, x, n)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.probdist import *
        >>> d = MVNormal()
        >>> print d.getLnPDF((0.0, 0.0))
        -1.83787706641
        
        The density function is (latex representation):
        
        f(x) = (2 \pi)^{-n/2} \det(\Sigma)^{-0.5} 
               \exp^{ -0.5 (x - \mu)^T \Sigma^{-1} (x-\mu)}
        
        where n is the number of dimensions, \mu is the mean vector and 
        \Sigma is the variance-covariance matrix. The log of f(x) is thus:
        
        \log[f(x)] = -0.5 \left[ n \log(2 \pi) + \log \det(\Sigma) 
                   + (x - \mu)^T \Sigma^{-1} (x-\mu) \right]
        
        For the bivariate standard normal example above, n = 2, \mu = (0,0),
        \Sigma = ((1,0),(0,1)), \det(\Sigma) = 1, and x-\mu = (0,0), so

        \log(f(x)) = -0.5 [ 2 (1.837877) + 0 + \exp^{0} ] = -1.837877
        
        """
        nparams = MVNormalDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to getLnPDF has length %d but length expected was %d' % (len(x), nparams)
        return MVNormalDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Because there is very little extra cost to computing a fully 
        normalized density for this distribution, this function simply
        returns getLnPDF(x).
        
        """
        nparams = MVNormalDistBase.getNParams(self)
        assert nparams == len(x), 'Vector supplied to getRelativeLnPDF has length %d but length expected was %d' % (len(x), nparams)
        return MVNormalDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, m, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean vector to m and the variance covariance matrix to v.
        The vector m should be a one-dimensional tuple of length dim, and
        the matrix v should be supplied as a two-dimensional tuple in which
        the lengths of both rows and columns is also dim.
        
        """
        return MVNormalDistBase.setMeanAndVariance(self, m, v)
        
    def _printSquareMatrix(self, matrix):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Used for debugging to print, say, the variance-covariance matrix .
        
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
        
    def getVarCovarMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns variance-covariance matrix in the form of a two-dimensional
        tuple.
        
        >>> from phycas.probdist import *
        >>> d = MVNormal(m=(1,2,3), v=((1,2,3),(2,3,4),(3,4,5)))
        >>> v = d.getVarCovarMatrix()
        >>> for i in range(3):
        ...     for j in range(3):
        ...         print '%5.1f' % v[i][j],
        ...     print
          1.0   2.0   3.0
          2.0   3.0   4.0
          3.0   4.0   5.0
        
        """
        dim = self.getNParams()
        v = MVNormalDistBase.getVarCovarMatrix(self)
        varcov = []
        k = 0
        for i in range(dim):
            tmp = []
            for j in range(dim):
                tmp.append(v[k])
                k += 1
            varcov.append(tmp)
        return varcov
        
    def getNParams(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of parameters in the MVNormal distribution. 
        
        """
        return MVNormalDistBase.getNParams(self)
    
    def debugMVNorm(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns string showing information about this object. 
        
        """
        return MVNormalDistBase.debugMVNorm(self)
    
