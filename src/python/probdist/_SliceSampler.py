from _ProbDistExt import *

class SliceSampler(SliceSamplerBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Implements the univariate slice sampler described in Neal, Radford M.
    2003. Slice sampling. Annals of Statistics 31:705-741. To use, create
    a SliceSampler object by passing in a Lot object (for pseudorandom
    number generation) and a function object representing the probability
    density from which to sample. See documentation under the __init__
    function for examples of how to initialize a SliceSampler object.
    Once the SliceSampler object is created, drawing samples is as easy
    as calling SliceSampler.sample(). Call SliceSampler.debugSample() to
    get details of the sampling process.

    """
    
    def __init__(self, rng, func):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes a new SliceSampler object from a random number generator
        object (rng) and a function object representing a (not necessarily
        normalized) log probability density function. The parameter rng should
        be a ProbDist.Lot object, and func must be either one of the standard
        univariate continuous probability distributions defined in the
        ProbDist module or an object derived from ProbDist.AdHocDensityBase.
        Note that discrete distributions (e.g. Binomial) or multivariate
        distributions (e.g. Dirichlet) cannot be used to create a
        SliceSampler object.

        Below is an example where the function represents one of the standard
        univariate continuous distributions (the Beta distribution):

        >>> from phycas.ProbDist import *
        >>> rng = Lot(13579)
        >>> func = Beta(2.0, 1.0)
        >>> print round(2.0/3.0, 12) # expected mean
        0.666666666667
        >>> ss = SliceSampler(rng, func)
        >>> samples = [ss.sample() for x in range(1000)]
        >>> obs_mean = sum(samples)/1000.0
        >>> print round(obs_mean, 12) # observed mean
        0.662050296713

        In the above example (as well as in all of the examples that follow),
        the random number generator was instantiated with a seed: e.g.
        rng = Lot(13579). This is to ensure that the answer you get from 
        running these examples is identical to the result reported in the 
        example output provided. Ordinarily, you would simply use rng =
        Lot(), which sets the starting seed using the system clock.
        The desire for consistency also explains the use of round(x, n) to
        round floating point values to n digits of precision. Using round()
        ensures that the output generated will not differ from what is
        presented in the example, even across different CPU architectures.

        Here is a more complicated example involving a bimodal mixture
        distribution composed of equal parts Beta(2,19) and
        Beta(19,2). Because this mixture distribution is not represented
        in any of the distributions provided, we must create a new class
        (MyBimodalDist) derived from ProbDist.AdHocDensityBase. Any function
        used to create a SliceSampler must define a function (here getLnPDF)
        that returns the natural log of the probability density for any
        supplied value. If the supplied value has zero probability density
        (or is outside the range of the variate), this function should return
        a very large negative number. The standalone function in ProbDist
        called getEffectiveLnZero is used here to store the largest possible
        negative number in the data member lnzero. The __call__ function must
        be defined so that SliceSampler can treat an instance of MyBimodalDist
        as if it were a function.
        
        >>> from phycas.ProbDist import *
        >>> import math
        >>> class MyBimodalDist(AdHocDensityBase):
        ...    def __init__(self):
        ...        # Very important: must call base class constructor!
        ...        AdHocDensityBase.__init__(self) 
        ...        self.b1 = Beta(2,19)
        ...        self.b2 = Beta(19,2)
        ...        self.lnzero = getEffectiveLnZero()
        ...
        ...    def getLnPDF(self, x):
        ...        if x > 0.0 and x < 1.0:
        ...            p1 = math.exp(self.b1.getLnPDF(x))
        ...            p2 = math.exp(self.b2.getLnPDF(x))
        ...            p = (0.5*p1 + 0.5*p2)
        ...            lnp = math.log(p)
        ...        else:
        ...            lnp = self.lnzero
        ...        return lnp
        ...
        ...    def __call__(self, x):
        ...        return self.getLnPDF(x)
        >>> ss = SliceSampler(Lot(13579), MyBimodalDist())
        >>> samples = [ss.sample() for i in range(100)]
        >>> for x in samples:
        ...    print round(x, 12) # doctest: +ELLIPSIS
        0.852238673297
        0.878160025763
        ...
        0.067862903854
        0.039882597226

        """
        SliceSamplerBase.__init__(self, rng, func)

    def setXValue(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets current value of slice sampler to the supplied value x.
        
        """
        SliceSamplerBase.setXValue(self, x)        
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a sample from the target distribution using slice sampling. This
        version is meant to be relatively efficient. Use debugSample() if you
        want to get some details of the slice sampling process.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(91735), Exponential(10.0))
        >>> print round(ss.sample(), 12)
        0.044758783301
        
        """
        return SliceSamplerBase.sample(self)
    
        
    def debugSample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a sample from the target distribution using slice sampling.
        Same as sample(), but returns information about the slice sampling
        process in the form of a list of which the elements are (log means
        logarithm base e):
             0: sampled x value
             1: x-coord of vertical line
             2: log of y-coord of top of vertical line
             3: x-coord of left edge of slice
             4: x-coord of right edge of slice
             5: log of y-coord of slice
             6: original slice interval width
            7+: x-coord of failed sampling attempts
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(91735), Exponential(10.0))
        >>> info = ss.debugSample()
        >>> print round(info[0], 12) # sampled value
        0.044758783301
        >>> print round(info[1], 12) # x-coord of vertical line
        0.1
        >>> print round(info[2], 12) # log(y-coord of top of vertical line)
        1.30258509299
        >>> print round(info[3], 12) # left edge of slice
        -0.249143314404
        >>> print round(info[4], 12) # right edge of slice
        0.481004689004
        >>> print round(info[5], 12) # log of slice y-coord
        0.97123250238
        >>> print round(info[6], 12) # slice unit width
        1.0
        >>> for x in info[7:]:       # failed samples
        ...    print round(x, 12)
        -0.249143314404

        Notes:
        1. In the above case, the original slice of width 1.0 was cropped on
           the left end because first sampling attempt failed. This is why the
           left edge of the slice is identical to the x-coord of the failed
           sample, and why right edge minus left edge is less than the slice
           unit width.
        2. The y-coord of the top of the vertical line at x-coord = 0.1 is
           -1.0 because constant part of density function is ignored. The
           log of the exponential density function is thus -lambda*x rather
           than log(lambda) - lambda*x
        
        """
        return SliceSamplerBase.debugSample(self)

    def overrelaxedSample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a sample from the target distribution using overrelaxed slice
        sampling. This sampling method is in general more expensive than
        sample(), and using overrelaxedSample() with a distribution that is
        not unimodal may raise exceptions. Thus, use sample() unless there is
        a demonstrated need for overrelaxed sampling.

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(91735), Exponential(10.0))
        >>> print round(ss.overrelaxedSample(), 12)
        0.033135032654
        
        """
        return SliceSamplerBase.overrelaxedSample(self)

    def debugOverrelaxedSample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a sample from the target distribution using overrelaxed slice
        sampling. Same as overrelaxedSample(), but returns information about
        the overrelaxed slice sampling process in the form of a list of which
        the elements are (log means logarithm base e):
             0: sampled x value
             1: x-coord of vertical line
             2: log of y-coord at top of vertical line
             3: x-coord of left edge of slice
             4: x-coord of right edge of slice
             5: log of y-coord of slice
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(91735), Exponential(10.0))
        >>> info = ss.debugOverrelaxedSample()
        >>> print "%.12f" % info[0] # sampled value
        0.033135032654
        >>> print "%.1f"  % info[1] # x-coord of vertical line
        0.1
        >>> print "%.1f"  % info[2] # log(y-coord of top of vertical line)
        1.3
        >>> print "%.12f" % info[3] # left edge of slice
        -0.000000381470
        >>> print "%.12f" % info[4] # right edge of slice
        0.133135414124
        >>> print "%.12f" % info[5] # log of slice y-coord
        0.971232502380

        Notes:
        1. In overrelaxed slice sampling, the left and right edges of the
           slice are determined precisely using bisection. This is why the
           slice unit width is not reported, as it is not used.
        2. *** begin again here *** discuss how 0.033135032654 is left edge
           plus right edge - initial x-coord
        
        """
        return SliceSamplerBase.debugOverrelaxedSample(self)
    
    def attachFunc(self, func):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Swaps in a new probability density function for the one specified at
        object construction. The supplied func must either be one of the
        standard univariate continuous probability distributions defined in
        the ProbDist module or an object derived from the class
        ProbDist.AdHocDensityBase. Note that discrete distributions (e.g.
        Binomial) or multivariate distributions (e.g. Dirichlet)
        cannot be used to create a SliceSampler object.

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(91735), Exponential(10.0))
        >>> print round(ss.sample(), 12)
        0.044758783301
        >>> func = Beta(100.0, 10.0)
        >>> ss.attachFunc(func)
        >>> print round(ss.sample(), 12)
        0.846151407446
        
        """
        SliceSamplerBase.attachFunc(self, func)

    def attachRandomNumberGenerator(self, rng):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Swaps in a new pseudorandom number generator (an instance of the class
        Lot) for the one installed at object creation.
        
        >>> from phycas.ProbDist import *
        >>> # Create SliceSampler with r1
        >>> r1 = Lot(73915)
        >>> ss = SliceSampler(r1, Exponential(0.1))
        >>> print round(ss.sample(), 12)
        0.35589724917
        >>> # Switch to r2
        >>> r2 = Lot(51937)
        >>> ss.attachRandomNumberGenerator(r2)
        >>> print round(ss.sample(), 12)
        8.54091425806
        >>> # Switch back to r1 and reset r1's seed to original value
        >>> ss.attachRandomNumberGenerator(r1)
        >>> r1.setSeed(73915)
        >>> print round(ss.sample(), 12)
        0.35589724917

        """
        SliceSamplerBase.attachRandomNumberGenerator(self, rng)

    def adaptSimple(self, multiplier):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets slice unit width to a value equal to the average cropped slice
        interval width times the multiplier provided. At least one sample must
        be drawn before calling this method (because otherwise there is no
        basis for computing the average cropped slice interval width).
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> ss.adaptSimple(1.0)
        Traceback (most recent call last):
            ...
        AssertionError: must draw at least one sample before calling adaptSimple()
        >>> samples = [ss.sample() for i in range(10)]
        >>> print ss.getSliceUnitWidth()
        1.0
        >>> ss.adaptSimple(1.0)
        >>> print ss.getSliceUnitWidth()
        24.3843066307
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling adaptSimple()'
        SliceSamplerBase.adaptSimple(self, multiplier)

    def adaptNeal(self, multiplier):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets slice unit width to the average distance between successive
        sampled values times the multiplier supplied. Suggested by Neal in 
        section 4.4 (p. 721). At least one sample must be drawn before
        calling this method (because otherwise there is no basis for
        computing the average distance between successive sampled values).

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> ss.adaptNeal(1.0)
        Traceback (most recent call last):
            ...
        AssertionError: must draw at least one sample before calling adaptNeal()
        >>> samples = [ss.sample() for i in range(10)]
        >>> print ss.getSliceUnitWidth()
        1.0
        >>> ss.adaptNeal(1.0)
        >>> print ss.getSliceUnitWidth()
        8.42916765004

        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling adaptNeal()'
        SliceSamplerBase.adaptNeal(self, multiplier)

    def getSliceUnitWidth(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns current slice unit width. A slice begins with a horizontal
        interval of this width randomly positioned around last sampled value.
        The slice is increased to the right and left by adding more units of
        this width until the density curve is crossed on each side.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> samples = [ss.sample() for i in range(10)]
        >>> print ss.getSliceUnitWidth()
        1.0
        >>> ss.adaptSimple(1.0)
        >>> print ss.getSliceUnitWidth()
        24.3843066307

        """
        return SliceSamplerBase.getSliceUnitWidth(self)

    def getMinX(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the minimum sampled value since last call to the
        resetDiagnostics() function. At least one sample must be drawn before
        calling this method, otherwise an exception is raised.

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> ss.getMinX()
        Traceback (most recent call last):
            ...
        AssertionError: must draw at least one sample before calling getMinX()
        >>> samples = [ss.sample() for i in range(10)]
        >>> print ss.getMinX()
        0.35589724917
        >>> print min(samples)
        0.35589724917
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getMinX()'
        return SliceSamplerBase.getMinX(self)

    def getMaxX(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the maximum sampled value since last call to the
        resetDiagnostics() function. At least one sample must be drawn before
        calling this method, otherwise an exception is raised.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> ss.getMaxX()
        Traceback (most recent call last):
            ...
        AssertionError: must draw at least one sample before calling getMaxX()
        >>> samples = [ss.sample() for i in range(10)]
        >>> print ss.getMaxX()
        40.4682348007
        >>> print max(samples)
        40.4682348007

        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getMaxX()'
        return SliceSamplerBase.getMaxX(self)

    def getOrigLeftEdgeOfSlice(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        For last sample taken, returns the original left edge of the slice.
        The slice is created by randomly positioning an interval of length
        slice unit width around the position of the last sampled point. Then
        the slice is extended by adding more intervals of length slice unit
        width to the left and right ends of the original unit until the
        slice edges are above the density function at their distal endpoints.
        This function returns the left edge of the slice at this point before
        any samples are attempted. The slice is shortened whenever there is
        a failed sampling attempt (point chosen is above the density). To get
        the left edge of the final slice (after a successful sampling), use
        the function getLeftEdgeOfSlice().

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(50.0, 50.0))
        >>> print round(ss.sample(), 12)
        0.591999204028
        >>> print round(ss.getOrigLeftEdgeOfSlice(), 12)
        -0.515610599787
        >>> print round(ss.getLeftEdgeOfSlice(), 12)
        -0.266608357228

        In the above example, it is clear that at least one sampling attempt
        failed because the original left edge of the slice is more extreme
        than the final left edge of the slice, indicating that the slice was
        cropped on the left end after a point (-0.266608357228) was rejected.
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getOrigLeftEdgeOfSlice()'
        return SliceSamplerBase.getOrigLeftEdgeOfSlice(self)

    def getOrigRightEdgeOfSlice(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        For last sample taken, returns the original right edge of the slice.
        The slice is created by randomly positioning an interval of length
        slice unit width around the position of the last sampled point. Then
        the slice is extended by adding more intervals of length slice unit
        width to the left and right ends of the original unit until the
        slice edges are above the density function at their distal endpoints.
        This function returns the right edge of the slice at this point before
        any samples are attempted. The slice is shortened whenever there is
        a failed sampling attempt (point chosen is above the density). To get
        the right edge of the final slice (after a successful sampling), use
        the function getRightEdgeOfSlice().

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(50.0, 50.0))
        >>> print round(ss.sample(), 12)
        0.591999204028
        >>> print round(ss.getOrigRightEdgeOfSlice(), 12)
        1.48438940021
        >>> print round(ss.getRightEdgeOfSlice(), 12)
        1.48438940021

        Unlike the example for getOrigLeftEdgeOfSlice(), in this case the
        original and final right slice edges are identical, meaning that no
        failed sampling attempts occurred on the right end of the slice.
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getOrigRightEdgeOfSlice()'
        return SliceSamplerBase.getOrigRightEdgeOfSlice(self)

    def getLeftEdgeOfSlice(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        For last sample taken, returns the final left edge of the slice.
        The slice is created by randomly positioning an interval of length
        slice unit width around the position of the last sampled point. Then
        the slice is extended by adding more intervals of length slice unit
        width to the left and right ends of the original unit until the
        slice edges are above the density function at their distal endpoints.
        The slice is cropped whenever there is a failed sampling attempt
        (point chosen is above the density). This function returns the left
        edge of the final slice (after a successful sampling). To find the
        original slice width before any cropping, use the function
        getLeftEdgeOfSlice().

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(50.0, 5.0))
        >>> print round(ss.sample(), 12)
        0.591999204028
        >>> print round(ss.getOrigLeftEdgeOfSlice(), 12)
        -0.515610599787
        >>> print round(ss.getLeftEdgeOfSlice(), 12)
        -0.266608357228

        It is clear in the above example that at least one attempted sample
        failed because the point -0.266608357228 was used to crop the slice
        on the left end.
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getLeftEdgeOfSlice()'
        return SliceSamplerBase.getLeftEdgeOfSlice(self)

    def getRightEdgeOfSlice(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        For last sample taken, returns the final right edge of the slice.
        The slice is created by randomly positioning an interval of length
        slice unit width around the position of the last sampled point. Then
        the slice is extended by adding more intervals of length slice unit
        width to the left and right ends of the original unit until the
        slice edges are above the density function at their distal endpoints.
        The slice is cropped whenever there is a failed sampling attempt
        (point chosen is above the density). This function returns the right
        edge of the final slice (after a successful sampling). To find the
        original slice width before any cropping, use the function
        getRightEdgeOfSlice().

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(5.0, 50.0))
        >>> print round(ss.sample(), 12)
        0.087950013611
        >>> print round(ss.getOrigRightEdgeOfSlice(), 12)
        0.484389400213
        >>> print round(ss.getRightEdgeOfSlice(), 12)
        0.202508966415

        It is clear in the above example that at least one attempted sample
        failed because the point 0.202508966415 was used to crop the slice
        on the right end.
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getRightEdgeOfSlice()'
        return SliceSamplerBase.getRightEdgeOfSlice(self)

    def getSliceYValue(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the height (on the log scale) of the slice for the last sample
        taken.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(5.0, 50.0))
        >>> print round(ss.sample(), 12)
        0.087950013611
        >>> print ss.getSliceYValue()
        1.6559680552

        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getSliceYValue()'
        return SliceSamplerBase.getSliceYValue(self)

    def getNumFuncEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns number of function evaluations required since object was
        constructed, or since the last call to resetDiagnostics(). Divide by
        the value returned by getNumSamples() to compute average number of
        function evaluations required for each sample drawn.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(10.0))
        >>> samples = [ss.sample() for i in range(10)]
        >>> evals = ss.getNumFuncEvals()
        >>> print round(evals/10.0, 12)
        8.0
        
        """
        return SliceSamplerBase.getNumFuncEvals(self)

    def getNumFailedSamples(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the total number of failed samples since the object was
        instantiated, or since the last call of resetDiagnostics(). Each time
        a sample is attempted but fails (because the sampled value is in an
        extreme part of the slice that lies above the density curve), that
        value is used to crop the slice. Divide by the result of calling
        getNumSamples() to compute the average number of failures per sample
        drawn.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(10.0))
        >>> print ss.getNumFailedSamples()
        0
        >>> samples = [ss.sample() for i in range(10)]
        >>> failures = ss.getNumFailedSamples()
        >>> print failures
        40
        >>> print round(failures/10.0, 12)
        4.0

        """
        return SliceSamplerBase.getNumFailedSamples(self)

    def getNumUnitsRequired(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the total number of slice units required to bracket all slices
        since object was constructed (or last call to resetDiagnostics()).
        Divide by the result of calling getNumSamples() to compute the average
        number of units required per sample drawn.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.2))
        >>> print ss.getNumUnitsRequired()
        0
        >>> samples = [ss.sample() for i in range(10)]
        >>> nunits = ss.getNumUnitsRequired()
        >>> print nunits
        85
        >>> print round(nunits/10.0, 12)
        8.5

        """
        return SliceSamplerBase.getNumUnitsRequired(self)

    def getNumSamples(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of samples drawn since the last call to
        resetDiagnostics(). This number can be used to calculate means of
        various quantities, such as the mean cropped slice interval width
        used by adaptSimple().

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> samples = [ss.sample() for i in range(10)]
        >>> print len(samples)
        10
        >>> print ss.getNumSamples()
        10
        
        """
        return SliceSamplerBase.getNumSamples(self)

    def resetDiagnostics(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets all diagnostic counters to 0. These counters are used by
        getNumSamples(), getNumUnitsRequired(), getNumFuncEvals(),
        getNumFailedSamples(), getMinX() and getMaxX(). It also affects the
        functions that set the slice unit width based on diagnostics gleaned
        from previous sampling, i.e. adaptSimple(), adaptNeal() and
        adaptYConditional().
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Exponential(0.1))
        >>> samples = [ss.sample() for i in range(10)]
        >>> print ss.getNumSamples()
        10
        >>> ss.resetDiagnostics()
        >>> print ss.getNumSamples()
        0

        """
        SliceSamplerBase.resetDiagnostics(self)

    def adaptYConditional(self, from_ends, multiplier):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function provides a way to continuously tune the slice unit
        width so that it is close to optimal for the distribution being
        sampled. Unlike the other methods, adaptSimple() and adaptNeal(),
        the slice unit width is not set just once, but instead is recalulated
        for each sample drawn based on the height (y) of the slice. The
        purpose of adaptYConditional() is to parameterize the function used
        to compute the slice unit width given y. The width of the density is
        computed using bisection at two points, one a fraction from_ends from
        the bottom and the other a fraction from_ends from the top of the
        density at its estimated mode (I say estimated here, because the
        position of the mode is simply the previously-sampled point having
        the highest density, which could innacurate to a degree depending on
        how much sampling has been done since the object was constructed (or
        since the last call to resetDiagnostics()). These two widths define a
        slope that can be used to crudely approximate the width of the density
        at any height. The parameter multiplier is stored for use by the
        calcW() function.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(20.0, 20.0))
        >>> samples = [ss.sample() for i in range(20)]
        >>> print ss.getNumFailedSamples()
        38
        >>> ss.adaptYConditional(0.1, 1.0)
        >>> ss.resetDiagnostics()
        >>> samples = [ss.sample() for i in range(20)]
        >>> print ss.getNumFailedSamples()
        14

        Note the fact that y-conditional adaptation improved the efficiency
        of the slice sampler in this case: the number of failed sample
        attempts out of 20 samplings was reduced. The improvement is somewhat
        deceptive, however. Most of the improvement comes from simply reducing
        the size of the slice unit width: the default value of 1.0 is way too
        large for this narrow Beta(20,20) density function.

        """
        SliceSamplerBase.adaptYConditional(self, from_ends, multiplier)

    def calcW(self, y):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Meaningful only if adaptYConditional() has been previously called.
        calcW() returns an approximation to the width of the density at height
        y, scaled by the multiplier previously supplied to the function
        adaptYConditional(). Note: y should not be on the log scale. You will
        probably never need to call this function explicitly except for
        debugging purposes.
        
        >>> from phycas.ProbDist import *
        >>> import math
        >>> ss = SliceSampler(Lot(73915), Beta(20.0, 20.0))
        >>> print ss.getSliceUnitWidth()
        1.0
        >>> samples = [ss.sample() for i in range(20)]
        >>> ss.adaptYConditional(0.1, 1.0)
        >>> samples = [ss.sample() for i in range(20)]
        >>> lny = ss.getLnDensityAtMode()
        >>> y = math.exp(lny)
        >>> w = ss.calcW(y/2.0) # approximate width halfway up density at mode
        >>> print round(w, 12)
        0.21875
        >>> w = ss.calcW(y/4.0) # approximate width one-fourth way up density
        >>> print round(w, 12)
        0.3017578125

        """
        return SliceSamplerBase.calcW(self, y)

    def getMode(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns current estimate of position of mode. It assumes that the
        density is unimodal; if density is actually multimodal, it will
        estimate the position of the highest mode. The estimation process is
        very simple: the estimated mode is that value having the highest
        density of all previously sampled values. If few samples have been
        drawn since the object was constructed (or since the resetDiagnostics()
        function was last called), then the accuracy of the estimate will be
        correspondingly diminished. Here is an example involving a highly
        variable distribution in which the true mode is located at 0.5.

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(2.0, 2.0))
        >>> first = ss.sample()
        >>> print round(ss.getMode(), 12)
        0.591999204028
        >>> next_ten = [ss.sample() for i in range(10)]
        >>> print round(ss.getMode(), 12)
        0.503134657448
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getMode()'
        return SliceSamplerBase.getMode(self)

    def getLastSampledXValue(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the last value sampled. The number of samples drawn must be
        greater than 0 in order for there to be a last sampled value, so
        attempting to call this function before calling sample will trip an
        assert.

        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(2.0, 2.0))
        >>> first = ss.sample()
        >>> print "%.5f" % first
        0.59200
        >>> lastx = ss.getLastSampledXValue()
        >>> print "%.5f" % lastx
        0.59200
        
        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getLastSampledXValue()'
        return SliceSamplerBase.getLastSampledXValue(self)

    def getLnDensityAtMode(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the natural logarithm of the density function at the estimated
        mode. The estimated mode can be obtained by calling getMode(). Here is
        the same example used for getMode(), but instead of asking for the
        mode, we ask for the log-density at the mode. Drawing one sample does
        not give us a very good estimate of the mode (0.591999204028 vs true
        mode 0.5). After drawing 10 more samples, our mode estimate has
        improved (0.503134657448) and thus the height of the density at that
        estimated mode is higher.
        
        >>> from phycas.ProbDist import *
        >>> ss = SliceSampler(Lot(73915), Beta(2.0, 2.0))
        >>> first = ss.sample()
        >>> print round(ss.getLnDensityAtMode(), 12)
        0.371023326913
        >>> next_ten = [ss.sample() for i in range(10)]
        >>> print round(ss.getLnDensityAtMode(), 12)
        0.405425803026

        """
        assert self.getNumSamples() > 0, 'must draw at least one sample before calling getLnDensityAtMode()'
        return SliceSamplerBase.getLnDensityAtMode(self)
