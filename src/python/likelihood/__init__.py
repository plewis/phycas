import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.conversions import *
from _LikelihoodExt import *
from _TreeLikelihood import *
from _Model import *
from _MCMCChainManager import *
from _SimData import *
from _TopoPriorCalculator import *
from _QMatrix import *
from _TreeLengthDist import *

#print 'importing likelihood...'

def testExamples(verbose = False):
    import doctest
    a = [0,0]

    if verbose: print '...testing examples in file _MCMCChainManager.py'
    r = doctest.testfile('_MCMCChainManager.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _Model.py'
    r = doctest.testfile('_Model.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _QMatrix.py'
    r = doctest.testfile('_QMatrix.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _SimData.py'
    r = doctest.testfile('_SimData.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _TopoPriorCalculator.py'
    r = doctest.testfile('_TopoPriorCalculator.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _TreeLikelihood.py'
    r = doctest.testfile('_TreeLikelihood.py')
    a[0] += r[0] ; a[1] += r[1]
    return tuple(a)

