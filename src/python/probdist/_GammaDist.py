from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Gamma(GammaDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate gamma probability distribution.

    Notes:
      - the gamma distribution has two parameters, the shape and scale
      - a gamma distribution with shape parameter 1 and scale parameter
        b is identical to an exponential distribution with mean b
      - the mean of a gamma distribution having shape a and scale b is
        a*b
      - the variance of a gamma distribution having shape a and scale b
        is a*b^2
      - the equation z = x/(x + y) can be used to generate z ~ Beta(a,b)
        given x ~ Gamma(shape=a, scale=1), y ~ Gamma(shape=b, scale=1).

    """
    def __init__(self, shape, scale):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the shape and scale parameters of the Gamma object. e.g.,

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getMean()
        6.0

        """
        GammaDistBase.__init__(self, shape, scale)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Gamma distribution.
        
        """
        return GammaDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the gamma distribution is continuous.
        
        """
        return GammaDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Gamma'
        
        """
        return GammaDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Gamma
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.__str__()
        Gamma(2.00000, 3.00000)
        
        """
        return GammaDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Gamma
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.__repr__()
        Gamma(2.00000, 3.00000)
        
        """
        return GammaDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Gamma(0.5, 2.0)
        >>> d2 = Gamma(5.0, 0.2)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample()
        0.000177184566
        >>> print "%.12f" % d2.sample()
        0.930716753542
        
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
        return GammaDistBase.setLot(self, lot)
        
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
        >>> d = Gamma(2, 3)
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.140064780761
        >>> print "%.12f" % d.sample()
        8.204670625640
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.140064780761

        """
        return GammaDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return GammaDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getMean()
        6.0

        """
        return GammaDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getVar()
        18.0

        """
        return GammaDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getStdDev()
        4.24264068712

        """
        return GammaDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 0.0
        up to the value x.

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getCDF(1.5)
        0.090204010431

        """
        return GammaDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the gamma distribution specified by
        this Gamma object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> d.setSeed(97531)
        >>> samples = [d.sample() for i in range(3)]
        >>> for x in samples:
        ...    print round(x,6)
        8.302218
        0.490945
        13.733702
        
        
        
        
        """
        return GammaDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> ln_pdf = d.getLnPDF(1.5)
        >>> print round(ln_pdf, 9)
        -2.291759469

        The gamma density is

        x^{a-1} * exp(-x/b)
        -------------------
           b^a * Gamma(a)
        
        where a = 2, b = 3, x = 1.5, and Gamma is the Gamma function, not
        the probability distribution. The density is thus

        (1.5)^{2-1} * exp(-1.5/3)   (1.5) * exp(-0.5)
        ------------------------- = ----------------- = 0.101088443
              3^2 * Gamma(2)             (9)(1)

        and the log density is thus -2.291759469
        Oops - looks like getLnPDF returns only relative PDF
        
        """
        return GammaDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the relative probability density function at the supplied
        value x. Returns the natural logarithm of the relative density at x.
        Use this function if speed is important but normalization is not, say
        in MCMC calculations. Use getLnPDF instead if you need to have a
        correctly normalized density value (i.e. from a density function that
        integrates to 1.0)

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getRelativeLnPDF(1.5)
        -0.0945348918918

        The gamma density is

        x^{a-1} * exp(-x/b)
        -------------------
           b^a * Gamma(a)
        
        where a = 2, b = 3, x = 1.5, and Gamma is the Gamma function, not
        the probability distribution. The log density is thus

        (a-1)*log(x) - x/b - a*log(b) - logGamma(a)
        
        The third and fourth terms are not needed if only the relative PDF
        is required because they do not contain x, so the return value 
        comprises just the first two terms: log(1.5) - 0.5, which
        equals -0.0945348918918
        
        """
        return GammaDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. 

        >>> from phycas.probdist import *
        >>> d = Gamma(2, 3)
        >>> print d.getMean()
        6.0
        >>> print d.getVar()
        18.0
        >>> d.setMeanAndVariance(5, 2)
        >>> print d.getMean()
        5.0
        >>> print d.getVar()
        2.0
        
        """
        return GammaDistBase.setMeanAndVariance(self, mean, var)
    
