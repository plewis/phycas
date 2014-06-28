from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Binomial(BinomialDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate binomial probability distribution.

    """
    def __init__(self, n, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the number of trials n and the probability p of a success on
        any given trial when initializing a Binomial object. e.g.,

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getMean()
        1.0

        """
        BinomialDistBase.__init__(self, n, p)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Binomial distribution.
        
        """
        return BinomialDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the binomial distribution is discrete.
        
        """
        return BinomialDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Binomial'
        
        """
        return BinomialDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Binomial
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.__str__()
        Binomial(10, 0.10000)
        
        """
        return BinomialDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Binomial
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.__repr__()
        Binomial(10, 0.10000)
        
        """
        return BinomialDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Binomial(10, 0.5)
        >>> d2 = Binomial(10, 0.2)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print [d1.sample(), d2.sample()]
        [7.0, 0.0]
        
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
        return BinomialDistBase.setLot(self, lot)
        
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
        >>> d = Binomial(10, 0.1)
        >>> d.setSeed(1537)
        >>> d.sample()
        2.0
        >>> d.sample()
        0.0
        >>> d.setSeed(1537)
        >>> d.sample()
        2.0

        """
        return BinomialDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return BinomialDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getMean()
        1.0

        """
        return BinomialDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getVar()
        0.9

        """
        return BinomialDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getStdDev()
        0.948683298051

        """
        return BinomialDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the sum of the probabilities of all positive whole
        numbers from 0 up to the value x.

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getCDF(3)
        0.9872048016

        """
        return BinomialDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the binomial distribution specified by
        this Binomial object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> d.setSeed(97531)
        >>> print [d.sample()] * 3
        [1.0, 1.0, 1.0]
        
        """
        return BinomialDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability of the supplied value x. Returns the natural
        logarithm of the probability of x. e.g.,

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getLnPDF(3)
        -2.8577871458

        The binomial probability is

            n!
        --------- p^y (1-p)^(n-y)
        y! (n-y)!
        
        where n = 10, p = 0.1 and y = 3. The probability is thus

         10!                    
        ----- (0.1)^3 (0.9)^7 = 0.057395628
        3! 7!                    

        and the log density is thus -2.857787145804875166398780295405
        Oops - looks like getLnPDF returns only relative PDF
        
        """
        return BinomialDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the relative probability at the supplied value x. Returns
        the natural logarithm of the relative density at x. Use this function
        if speed is important but normalization is not, say in MCMC
        calculations. Use getLnPDF instead if you need to have a correctly
        normalized density value (i.e. from a density function that integrates
        to 1.0)

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getRelativeLnPDF(3)
        -7.64527888859
        
        The binomial probability is

            n!
        --------- p^y (1-p)^(n-y)
        y! (n-y)!
        
        where n = 10, p = 0.1 and y = 3. Leaving off terms that do not involve
        the parameter p provides the relative PDF:

        p^y (1-p)^(n-y)

        The log of this is y*log(p) + (n-y)*log(1-p), which for this case
        equals -7.6452788885869211606464812299283
        
        """
        return BinomialDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. The mean is n*p
        and thus p can be computed as mean/n. This distribution is determined
        entirely by the mean and sample size (which was set when the BinomDist
        object was first created), so the var argument of this function is
        simply ignored. The reason this function requires both mean and
        variance is for compatibility with functions of the same name in
        other distributions.

        >>> from phycas.probdist import *
        >>> d = Binomial(10, 0.1)
        >>> print d.getMean()
        1.0
        >>> print d.getVar()
        0.9
        >>> d.setMeanAndVariance(3, 0)
        >>> print d.getMean()
        3.0
        >>> print d.getVar()
        2.7
        
        """
        return BinomialDistBase.setMeanAndVariance(self, mean, var)
    
