import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.conversions import *
from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *
from _Lot import *
from _StopWatch import *
from _SliceSampler import *
from _NormalDist import *
from _LognormalDist import *
from _GammaDist import *
from _ExponentialDist import *
from _InverseGammaDist import *
from _UniformDist import *
from _ImproperUniformDist import *
from _BinomialDist import *
from _BernoulliDist import *
from _BetaDist import *
from _BetaPrimeDist import *
from _DirichletDist import *
from _MVNormalDist import *
from _RelRateDist import *
from _SquareMatrix import *
from _RectangularMatrix import *

#print 'importing probdist...'

def testExamples(verbose = False):
    import doctest
    a = [0,0]

    if verbose: print '...testing examples in file _BernoulliDist.py'
    r = doctest.testfile('_BernoulliDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _BetaDist.py'
    r = doctest.testfile('_BetaDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _BetaPrimeDist.py'
    r = doctest.testfile('_BetaPrimeDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _BinomialDist.py'
    r = doctest.testfile('_BinomialDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _DirichletDist.py'
    r = doctest.testfile('_DirichletDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _RelRateDist.py'
    r = doctest.testfile('_RelRateDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _TreeLengthDist.py'
    r = doctest.testfile('_TreeLengthDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _MVNormalDist.py'
    r = doctest.testfile('_MVNormalDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _ExponentialDist.py'
    r = doctest.testfile('_ExponentialDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _NormalDist.py'
    r = doctest.testfile('_NormalDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _LognormalDist.py'
    r = doctest.testfile('_LognormalDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _GammaDist.py'
    r = doctest.testfile('_GammaDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _InverseGammaDist.py'
    r = doctest.testfile('_InverseGammaDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _Lot.py'
    r = doctest.testfile('_Lot.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _StopWatch.py'
    r = doctest.testfile('_StopWatch.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _SliceSampler.py'
    r = doctest.testfile('_SliceSampler.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _UniformDist.py'
    r = doctest.testfile('_UniformDist.py')
    a[0] += r[0] ; a[1] += r[1] 

    if verbose: print '...testing examples in file _ImproperUniformDist.py'
    r = doctest.testfile('_ImproperUniformDist.py')
    a[0] += r[0] ; a[1] += r[1] 
    return tuple(a)
    


