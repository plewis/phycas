import warnings
warnings.filterwarnings('ignore', '.*second conversion method ignored.*', RuntimeWarning)

#from phycas.conversions import *
from _DataMatrixExt import *
from _DataMatrix import *

#print 'importing datamatrix...'

def testExamples():
    import doctest
    return doctest.testfile('_DataMatrix.py')
