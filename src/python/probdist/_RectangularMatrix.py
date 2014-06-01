from _ProbDistExt import *

class RectangularMatrix(RectangularMatrixBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a rectangular matrix of floating point values (underlying
    C++ implementation stores these as doubles).
        
    """
    def __init__(self, nrows, ncols, value = 0.0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Create a rectangular matrix of size nrows by ncols containing value 
        in every cell.
        
        """
        RectangularMatrixBase.__init__(self, nrows, ncols, value)
        
    def getNRows(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of rows in the matrix.
        
        """
        return RectangularMatrixBase.getNRows(self)
        
    def getNCols(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of columns in the matrix.
        
        """
        return RectangularMatrixBase.getNCols(self)
        
    def getDimensions(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple containing the number of rows and number of columns in
        th matrix.
        
        """
        return (RectangularMatrixBase.getNRows(self),RectangularMatrixBase.getNCols(self))
        
    def getElement(self, i, j):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns (i,j)th element of the rectangular matrix.
        
        """
        return RectangularMatrixBase.getElement(self, i, j)
        
    def setElement(self, i, j, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets (i,j)th element of the rectangular matrix to value v.
        
        """
        RectangularMatrixBase.setElement(self, i, j, v)
        
    def getRow(self, i):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns ith row of the rectangular matrix.
        
        """
        return RectangularMatrixBase.getRow(self, i)
        
    def setRow(self, i, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets ith row of the rectangular matrix to the values in v.
        
        """
        RectangularMatrixBase.setRow(self, i, v)
        
    def getMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns rectangular matrix in the form of a two-dimensional list.
        
        """
        nrows, ncols = self.getDimensions()
        v = RectangularMatrixBase.getMatrix(self)
        m = []
        k = 0
        for i in range(nrows):
            tmp = []
            for j in range(ncols):
                tmp.append(v[k])
                k += 1
            m.append(tmp)
        return m
        
    def setMatrixFromFlattenedList(self, nrows, ncols, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Replaces existing or creates a new square matrix using the supplied
        unidimensional list or tuple v. The supplied list v is expected to 
        have length equal to nrows*ncols.
        
        """
        RectangularMatrixBase.setMatrix(self, nrows, ncols, v)
        
    def setMatrix(self, m):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Replaces existing or creates a new square matrix using the supplied
        two-dimensional list or tuple m.
        
        """
        nrows = len(m)
        ncols = len(m[0])
        v = []
        for row in m:
            for col in row:
                v.append(col)
        RectangularMatrixBase.setMatrix(self, nrows, ncold, v)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Represents matrix as string.
        
        """
        s = RectangularMatrixBase.__repr__(self)
        return s
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns tuple of length ncols representing the mean over rows.
        
        """
        return RectangularMatrixBase.getMean(self)
        
    def getVarCovMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns SquareMatrix representing the variance-covariance matrix,
        treating rows as independent observations and columns as variables.
        
        """
        return RectangularMatrixBase.getVarCovMatrix(self)
        