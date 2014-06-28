import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.conversions import *
from _ReadNexusExt import *
from _NexusReader import *

#print 'importing readnexus...'

def testExamples(verbose = False):
    import doctest

    if verbose: print '...testing examples in file _NexusReader.py'
    return doctest.testfile('_NexusReader.py')
