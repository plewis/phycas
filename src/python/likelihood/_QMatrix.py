from _LikelihoodExt import *

class QMatrix(QMatrixBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    The QMatrix class encapsulates the instantaneous rate matrix
    representing a substitution model. It is most useful for models that
    do not have closed-form expressions representing transition
    probabilities. The models that can be represented by QMatrix are
    those having n*(n-1)/2 relative rates for n states; that is, the
    relative rate of a change from state i to state j equals the relative
    rate of a change from state j to state i. Models can additionally
    have unequal state frequencies. The class will be described in the
    sequel using the GTR model as an example. The six relative rates of
    the GTR model are symbolized r_AC, r_AG, r_AT, r_CG, r_CT and r_GT
    (or more generically as r_ij), and the four base frequencies will be
    symbolized pi_A, pi_C, pi_G and pi_T (or pi_i).

    Ordinarily, the off-diagonal element for row i, column j, of the Q
    matrix for the GTR model is set to pi_j*r_ij and the diagonal elements
    are set equal to the negative sum of the other elements in the same
    row. To obtain the transition probability matrix P, Q is diagonalized
    yielding a 4x4 matrix V of eigenvectors and a one-dimensional array
    of eigenvalues. If the eigenvalues are used to form a diagonal matrix
    L, the matrix Q could be reconstituted as follows:

    Q = V L V'

    where V' is the transpose of the eigenvector matrix. The P matrix is
    computed as follows:

    P = exp(Q*t) = V exp(L*t) V'
    
    where exp(L*t) is obtained by replacing each eigenvalue lambda_i on
    the diagonal with exp(lambda_i*t). The value t is the edge length
    scaled by dividing by the sum of the off diagonal elements of the
    Q matrix.
    
    One difficulty is that Q is not a symmetric matrix if the nucleotide
    frequencies are unequal. The diagonalization is simplified by using
    a modified Q matrix (call it Z) in which each off-diagonal element
    is replaced by r_ij*sqrt(pi_i)*sqrt(pi_j), making Z symmetric. (The
    diagonal elements of Z remain the same as those of Q.) This is
    equivalent to

    Z = pi^(1/2) Q pi^(-1/2)

    The change probability matrix P will then be obtained as

    P = pi^(-1/2) exp(Z*t) pi^(1/2). This P matrix is identical to
    exp(Q*t).
     
    """

    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The constructor takes no arguments and initializes the object to
        represent the Jukes-Cantor model.
        
        """
        QMatrixBase.__init__(self)
        
    def getDimension(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of rows (columns) of the stored Q matrix. This
        equals the number of states in the model.
        
        """
        return QMatrixBase.getDimension(self)
        
    def setRelativeRates(self, relrates):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the relative rates used to construct the underlying Q matrix.
        For a k-state model, relrates should contain k*(k-1)/2 relative
        rates. For the GTR model, the expected order is rAC, rAG, rAT, rCG,
        rCT and rGT. In addition to relative rates, state frequencies are
        needed to fully specify the Q matrix. The state frequencies are
        supplied using the setStateFreqs function.
        
        """
        QMatrixBase.setRelativeRates(self, relrates)
        
    def setStateFreqs(self, freqs):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the state frequencies that, in combination with relative rates,
        specify the Q matrix. The supplied freqs should contain one frequency
        for every state. The supplied frequencies will be normalized so that
        they sum to 1.0. In addition to state frequencies, relative rates are
        needed to fully specify the Q matrix. The relative rates are supplied
        using the setRelativeRates function.
        
        """
        QMatrixBase.setStateFreqs(self, freqs)

    # Uncomment if numarray re-instated        
    #def getPMatrix(self, edgelen):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Returns a numarray representing the matrix of transition
    #    probabilities for an edge length equal to the supplied value edgelen.
    #    The following example sets up an HKY model with kappa = 4 and base
    #    frequencies piA = 0.1, piC = 0.2, piG = 0.3 and piT = 0.4, then
    #    calculates and prints out the transition probability matrix for an
    #    edgelength equal to 0.1.
    #    
    #    >>> import math
    #    >>> import numarray
    #    >>> from phycas import *
    #    >>> 
    #    >>> qmatrix = likelihood.QMatrix()
    #    >>> 
    #    >>> rr = [1.0,4.0,1.0,1.0,4.0,1.0]
    #    >>> qmatrix.setRelativeRates(rr)
    #    >>> 
    #    >>> pi = [0.1, 0.2, 0.3, 0.4]
    #    >>> qmatrix.setStateFreqs(pi)
    #    >>> 
    #    >>> print qmatrix.getPMatrix(0.1)
    #    [[ 0.87734732  0.01417824  0.08011797  0.02835647]
    #     [ 0.00708912  0.86650039  0.02126735  0.10514313]
    #     [ 0.02670599  0.01417824  0.9307593   0.02835647]
    #     [ 0.00708912  0.05257157  0.02126735  0.91907196]]
    #
    #    """
    #    return QMatrixBase.getPMatrix(self, edgelen)
        
    # Comment out this version if numarray is re-introduced
    def printSquareMatrix(self, matrix):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Utility function used to interpret a vector as a square matrix and
        print it out.
        
        """
        total_len = len(matrix)
        import math
        row_len = int(math.sqrt(total_len))
        assert row_len*row_len == total_len, 'Attempting to print a matrix that is not square'
        k = 0
        for i in range(row_len):
            for j in range(row_len):
                print '% .8f ' % matrix[k],
                k += 1
            print
        
    # Comment out if numarray re-instated        
    def getPMatrix(self, edgelen):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a vector representing the matrix of transition probabilities
        for an edge length equal to the supplied value edgelen. The following
        example sets up an HKY model with kappa = 4 and base frequencies
        piA = 0.1, piC = 0.2, piG = 0.3 and piT = 0.4, then calculates and
        prints out the transition probability matrix for an edgelength equal
        to 0.1.
        
        >>> import math
        >>> from phycas import *
        >>> 
        >>> qmatrix = likelihood.QMatrix()
        >>> 
        >>> rr = [1.0,4.0,1.0,1.0,4.0,1.0]
        >>> qmatrix.setRelativeRates(rr)
        >>> 
        >>> pi = [0.1, 0.2, 0.3, 0.4]
        >>> qmatrix.setStateFreqs(pi)
        >>> 
        >>> qmatrix.printSquareMatrix(qmatrix.getPMatrix(0.1))
         0.87734732   0.01417824   0.08011797   0.02835647 
         0.00708912   0.86650039   0.02126735   0.10514313 
         0.02670599   0.01417824   0.93075930   0.02835647 
         0.00708912   0.05257157   0.02126735   0.91907196 
        
        """
        return QMatrixBase.getPMatrix(self, edgelen)
        
    def getQMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a numarray representing the Q matrix.
        
        """
        return QMatrixBase.getQMatrix(self)
        
    def getEigenVectors(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a numarray representing the matrix.
        
        """
        return QMatrixBase.getEigenVectors(self)
        
    def getEigenValues(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a numarray representing the Q matrix.
        
        """
        return QMatrixBase.getEigenValues(self)
        
