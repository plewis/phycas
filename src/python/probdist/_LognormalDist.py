from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Lognormal(LognormalDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate lognormal probability distribution.

    """
    def __init__(self, logm, logstddev):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the mean and standard deviation of the log of a lognormal
        random variable when creating a new Lognormal object. e.g.,  
        X ~ Lognormal(0.0, 1.0) results in a distribution such that the mean
        of log(X) is 0.0 and the standard deviation of log(X) is 1.0.

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(0.0, 1.0)
        >>> print '%.5f' % d.getMean()
        1.64872
        >>> print '%.5f' % d.getStdDev()
        2.16120

        """
        LognormalDistBase.__init__(self, logm, logstddev)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Lognormal distribution.
        
        """
        return LognormalDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the lognormal distribution is continuous.
        
        """
        return LognormalDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Lognormal'
        
        """
        return LognormalDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Lognormal
        object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> print d.__str__()
        Lognormal(2.00000, 3.00000)
        
        """
        return LognormalDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Lognormal
        object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> print d.__repr__()
        Lognormal(2.00000, 3.00000)
        
        """
        return LognormalDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Lognormal(0.0, 1.0)
        >>> d2 = Lognormal(10.0, 1.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample()
        0.099890840821
        >>> print "%.12f" % d2.sample()
        21849.314570637751
        
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
        return LognormalDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.000730657606
        >>> print "%.12f" % d.sample()
        60.088784571183
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.000730657606

        """
        return LognormalDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return LognormalDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution). The mean of a Lognormal(logmean, logsd) 
        random variable is exp(logmean + logsd^2/2)

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> print '%.5f' % d.getMean()
        665.14163

        """
        return LognormalDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution). The variance of a Lognormal(logmean, 
        logsd) random variable is (exp(logsd^2) - 1)*exp(2*logmean + logsd^2).

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> print '%.5f' % d.getVar()
        3584470432.73958

        """
        return LognormalDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> print '%.5f' % d.getStdDev()
        59870.44707

        """
        return LognormalDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 
        negative infinity up to the value x.

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(0, 3)
        >>> print '%.5f' % d.getCDF(0.0)
        0.00000

        """
        return LognormalDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single value from the lognormal distribution specified
        by this Lognormal object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(0, 3)
        >>> d.setSeed(97531)
        >>> samples = [d.sample() for i in range(3)]
        >>> for x in samples:
        ...    print '%.5f' % x
        8.59350
        0.00115
        113.63627

        """
        return LognormalDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(0, 3)
        >>> ln_pdf = d.getLnPDF(1)
        >>> print '%.5f' % ln_pdf
        -2.01755

        The lognormal density f(x) is

                1              / -(log(x) - logmean)^2 \
        ------------------ exp| ---------------------- |
        x*logsd*sqrt(2*pi)     \      2*logsd^2        /
        
        If logmean = 0.5, logsd = 3.0 and x = 1.0, the density would thus be

            1        /  -1  \    1.013985788
        -------- exp| -----  | = ----------- = 0.134840601
        7.51988      \ 4*18 /    7.519884824

        The natural log (i.e., logarithm base e) density is

        (log(x) - logmean)^2
        -------------------- - log(x) - log(logsd) - 0.5*log(2*pi)
               2*logsd^2

        and thus the log density for the above example is

        -(0.5)^2/18.0 - log(1) - log(3.0) - 0.5*log(2*3.1415927)
        = -2.00366
        
        """
        return LognormalDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simply returns the log of the probability density function (i.e. 
        this function is identical to getLnPDF.
                
        """
        return LognormalDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. 

        >>> from phycas.ProbDist import *
        >>> d = Lognormal(2, 3)
        >>> print '%.5f' % d.getMean()
        665.14163
        >>> print'%.5f' %  d.getVar()
        3584470432.73958
        >>> d.setMeanAndVariance(5, 4)
        >>> print '%.5f' % d.getMean()
        5.00000
        >>> print '%.5f' % d.getVar()
        4.00000
        
        """
        return LognormalDistBase.setMeanAndVariance(self, mean, var)
    
