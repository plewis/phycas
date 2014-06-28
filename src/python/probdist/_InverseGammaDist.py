from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class InverseGamma(InverseGammaDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate inverse gamma probability distribution.

    Notes:
      - the inverse gamma distribution with shape a and scale b is the
        distribution of the inverse of a gamma random variable having
        shape a and scale b
      - the mean of an inverse gamma distribution having shape a and scale
        b is 1/[b*(a-1)]; the mean is undefined unless a > 1
      - the variance of an inverse gamma distribution having shape a and
        scale b is 1/[b^2*(a-1)^2*(a-2)]; the variance is undefined
        unless a > 2
        
    """
    def __init__(self, shape, scale):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the shape and scale parameters of the InverseGamma
        object. e.g.,

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.getMean()
        4.0

        Trying to create an InverseGamma object with a shape less than
        2 results in a VAlueError exception:

        >>> from phycas.probdist import *
        >>> d = InverseGamma(0.5, 3)
        Traceback (most recent call last):
        ...
        Exception: variance undefined for shape less than or equal to 2

        """
        InverseGammaDistBase.__init__(self, shape, scale)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Inverse Gamma distribution.
        
        """
        return InverseGammaDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the inverse gamma distribution is continuous.
        
        """
        return InverseGammaDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'InverseGamma'
        
        """
        return InverseGammaDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another InverseGamma
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.__str__()
        InverseGamma(2.25000, 0.20000)
        
        """
        return InverseGammaDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        InverseGamma object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.__repr__()
        InverseGamma(2.25000, 0.20000)
        
        """
        return InverseGammaDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = InverseGamma(2.25, 0.2)
        >>> d2 = InverseGamma(2.1, 10)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample()
        23.288879229004
        >>> print "%.12f" % d2.sample()
        0.056588582432
        
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
        return InverseGammaDistBase.setLot(self, lot)
        
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
        >>> d = InverseGamma(2.25, 0.2)
        >>> d.setSeed(135)
        >>> print d.sample()
        67.7996039852
        >>> print d.sample()
        1.63982940907
        >>> d.setSeed(135)
        >>> print d.sample()
        67.7996039852

        """
        return InverseGammaDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return InverseGammaDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.getMean()
        4.0

        """
        return InverseGammaDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.getVar()
        64.0

        """
        return InverseGammaDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.getStdDev()
        8.0

        """
        return InverseGammaDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 0.0
        up to the value x.

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.getCDF(1.5)
        0.198360740619

        """
        return InverseGammaDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the inverse gamma distribution specified by
        this InverseGamma object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> d.setSeed(97531)
        >>> print "%.12f" % d.sample()
        1.621624150641
        >>> print "%.12f" % d.sample()
        21.959303981216
        >>> print "%.12f" % d.sample()
        1.005590399350
        
        """
        return InverseGammaDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.probdist import *
        >>> d = InverseGamma(3, 2)
        >>> print d.getLnPDF(1.5)
        -4.72778248801

        The inverse gamma density is

        (1/x)^{a+1} * exp(-1/(x*b))
        ---------------------------
               b^a * Gamma(a)
        
        where a = 3, b = 2, x = 1.5, and Gamma is the Gamma function, not
        the probability distribution. The density is thus

        (1/1.5)^{3+1} * exp(-1/(1.5*2))   (2/3)^4 * exp(-1/3)
        ------------------------------- = ------------------- = 0.008846065
              (2)^(3) * Gamma(3)                  8*2

        and the log density is thus -4.727782488
        
        """
        return InverseGammaDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the natural logarithm of the relative probability density
        function at the supplied value x. Use this function if speed is
        important but normalization is not, say in MCMC calculations. Use
        getLnPDF instead if you need to have a correctly normalized density
        value (i.e. from a density function that integrates to 1.0)

        >>> from phycas.probdist import *
        >>> d = InverseGamma(3, 2)
        >>> print d.getRelativeLnPDF(1.5)
        -1.95519376577

        The inverse gamma density is

        (1/x)^{a+1} * exp(-1/(x*b))
        ---------------------------
               b^a * Gamma(a)
        
        where a = 3, b = 2, x = 1.5, and Gamma is the Gamma function, not
        the probability distribution. The log density is thus

        -(a+1)*log(x) - 1/(x*b) - a*log(b) - logGamma(a)
        
        The third and fourth terms are not needed if only the relative PDF
        is required because they do not contain x, so the return value 
        comprises just the first two terms: -4*log(1.5) - 1/3, which
        equals -1.95519376577
        
        """
        return InverseGammaDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. To create an inverse
        gamma distribution with specific mean and variance m and v,
        respectively, first set z = m^2/v. The shape parameter a = 2+z and
        the scale parameter b = 1/(m + m*z). Suppose you wish the mean to be
        5 and the variance 2.

        z = 5^2/2 = 25/2 = 12.5
        shape a = 2 + 12.5 = 14.5
        scale b = 1/(5 + 5*12.5) = 1/67.5 = 0.014814815

        >>> from phycas.probdist import *
        >>> d = InverseGamma(2.25, 0.2)
        >>> print d.getMean()
        4.0
        >>> print d.getVar()
        64.0
        >>> d.setMeanAndVariance(5, 2)
        >>> print d.getMean()
        5.0
        >>> print d.getVar()
        2.0
        >>> d = InverseGamma(14.5, 0.014814815)
        >>> print d.getMean()
        4.9999999375
        >>> print d.getVar()
        1.99999995

        """
        return InverseGammaDistBase.setMeanAndVariance(self, mean, var)
    
