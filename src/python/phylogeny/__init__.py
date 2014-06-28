import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.conversions import *
from _PhylogenyExt import *
from _Tree import *
from _TreeNode import *
from _TreeManip import *
from _Split import *

#print 'importing phylogeny...'

def testExamples(verbose = False):
    import doctest
    a = [0,0]

    if verbose: print '...testing examples in file _Tree.py'
    r = doctest.testfile('_Tree.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _TreeManip.py'
    r = doctest.testfile('_TreeManip.py')
    a[0] += r[0] ; a[1] += r[1]

    if verbose: print '...testing examples in file _Split.py'
    r = doctest.testfile('_Split.py')
    a[0] += r[0] ; a[1] += r[1]
    return tuple(a)
