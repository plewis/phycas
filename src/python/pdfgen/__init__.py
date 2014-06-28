import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

from phycas.conversions import *
#from _PDFGen import *
from _PDFGenerator import *

def testExamples(verbose = False):
    import doctest
    if verbose: print '...testing examples in file _PDFGenerator.py'
    return doctest.testfile('_PDFGenerator.py')
