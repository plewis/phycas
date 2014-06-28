from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Uniform(UniformDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate uniform probability distribution. This
    distribution has two parameters, a and b, representing the left and
    right bounds, respectively, on the real line. The mean is simply the
    center of the interval defined by a and b, and the variance is
    (b-a)^2/12.

    """
    def __init__(self, a, b):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the left and right bounds of the Uniform object. e.g.,

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getMean()
        2.5

        """
        UniformDistBase.__init__(self, a, b)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Uniform distribution.
        
        """
        return UniformDistBase.clone(self)

    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the uniform distribution is continuous.
        
        """
        return UniformDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Uniform'
        
        """
        return UniformDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Uniform
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.__str__()
        Uniform(2.00000, 3.00000)
        
        """
        return UniformDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Uniform
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.__repr__()
        Uniform(2.00000, 3.00000)
        
        """
        return UniformDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Uniform(2.0, 3.0)
        >>> d2 = Uniform(1.0, 3.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample()
        2.010620383045
        >>> print "%.12f" % d2.sample()
        1.993557016767
        
        In this example, both d1 and d2 use the same random number generator
        when their sample functions are called. Generating all random numbers
        from a single random number generator has the advantage of allowing
        an entire analysis to be recreated by simply specifying the original
        seed (here, 1357) to the master Lot object. If setLot is not used,
        each distribution object maintains its own random number generator,
        which is initialized using the system clock at the time the object is
        created. This behavior is convenient, but makes it very hard to
        replicate results.

        """
        return UniformDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect. e.g.,

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> d.setSeed(135)
        >>> print d.sample()
        2.00105655985
        >>> print d.sample()
        2.75760139644
        >>> d.setSeed(135)
        >>> print d.sample()
        2.00105655985

        """
        return UniformDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return UniformDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getMean()
        2.5

        """
        return UniformDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getVar()
        0.0833333333333

        """
        return UniformDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getStdDev()
        0.288675134595

        """
        return UniformDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from the
        left bound up to the value x.

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getCDF(2.5)
        0.5

        """
        return UniformDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the uniform distribution specified by
        this Uniform object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> d.setSeed(97531)
        >>> print [d.sample()] * 3
        [2.763313617332142, 2.763313617332142, 2.763313617332142]
        
        """
        return UniformDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.probdist import *
        >>> import math
        >>> d = Uniform(2, 3)
        >>> print math.fabs(d.getLnPDF(2.5))  # see note about fabs below
        0.0
        >>> d = Uniform(1, 3)
        >>> print d.getLnPDF(2.5)
        -0.69314718056

        The uniform density is 1/(b-a), where a is the left bound and b the
        right bound. In the second case above, b-a = 2, so the density is
        0.5 everywhere. The log of 0.5 is -0.69314718056, and this is what is
        returned. This function returns the most negative floating point
        number possible if x is out of bounds.
        
        You might wonder why math.fabs is needed in the first case illustrated
        in the above example. The answer is that on some platforms, -1.0*0.0
        equals -0.0, where as on other platforms, -1.0*0.0 equals 0.0!
        Python's doctest module is used to check all examples given in
        documentation like this, so if math.fabs were not used, doctest would
        say that the above example fails on platforms in which -1.0*0.0 equals
        -0.0. In this case, the log of the uniform(2,3) probability density is
        log(1/(3-2)) = log(0) - log(1) = -log(1), which is 0.0 on some
        platforms and -0.0 on others.

        """
        return UniformDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the natural logarithm of the relative probability density
        function at the supplied value x. This value is always 0.0 because
        the density function is constant. Use getLnPDF if you want instead
        the log of the actual normalized density function. Use this one if
        speed is critical but normalization is not (e.g. MCMC). This function
        returns the most negative floating point number possible if x is out
        of bounds.

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getRelativeLnPDF(1.5)
        -1.79769313486e+308
        >>> print d.getRelativeLnPDF(2.5)
        0.0

        """
        return UniformDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. The left and right
        bound parameters (a and b) can be calculated from mean and var as
        follows:

        c = b-a = sqrt(12*var)
        a = mean - c/2
        b = a + c

        >>> from phycas.probdist import *
        >>> d = Uniform(2, 3)
        >>> print d.getMean()
        2.5
        >>> print d.getVar()
        0.0833333333333
        >>> d.setMeanAndVariance(5, 2)
        >>> print d.getMean()
        5.0
        >>> print d.getVar()
        2.0
        
        """
        return UniformDistBase.setMeanAndVariance(self, mean, var)
    
