from _PhylogenyExt import *

class TreeManip(TreeManipBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    A manipulator of Tree objects. Can be used to create a tree de novo,
    or rearrange the topology of an existing tree.

    """
    def __init__(self, t):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Constructs a TreeManip object that operates on the supplied tree t.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> d = ProbDist.Exponential(10.0)
        >>> tm.starTree(5, d)
        >>> print t.walkPreorder()
        (0) -> [5] -> (4) -> (3) -> (2) -> (1)

        """
        TreeManipBase.__init__(self, t)

#    def buildTreeFromSplitVector(self, split_vect, edge_len_dist):
#         #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
#         """
#         Creates a tree having the splits defined in split_vect, which should
#         be a list or tuple of string representations of splits. All splits in
#         split_vect should be compatible. Later splits not compatible with
#         earlier ones already in the tree will be ignored.
#         
#         >>> from phycas import *
#         >>> t = Phylogeny.Tree()
#         >>> tm = Phylogeny.TreeManip(t)
#         >>> v = ['--****', '--**--', '----**']
#         >>> tm.buildTreeFromSplitVector(v, ProbDist.Exponential(10.0))
#         >>> print t.walkPreorder(2)
#         (0) -> [6] -> (1) -> [4294967295] -> [4294967295] -> (3) -> (2) -> [4294967295] -> (5) -> (4)
# 
#         """
#        TreeManipBase.buildTreeFromSplitVector(self, split_vect, edge_len_dist)

    def starTree(self, num_tips, edge_len_dist = None):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a star tree having num_tips tips and edge lengths drawn at
        random from the supplied edge_len_dist probability distribution.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> d = ProbDist.Exponential(10.0)
        >>> tm.starTree(5, d)
        >>> print t.walkPreorder()
        (0) -> [5] -> (4) -> (3) -> (2) -> (1)

        """
        return TreeManipBase.starTree(self, num_tips, edge_len_dist)

    def equiprobTree(self, num_tips, rng, internal_dist = None, external_dist = None):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a random tree from a discrete uniform distribution. For 
        example, for 6 taxa there are 105 distinct unrooted tree topologies,
        and if this function were to be called many times, it would generate
        each of these 105 tree topologies with probability 1/105 = 0.0095238.
        The probability distribution internal_dist is used to set internal
        edge lengths, whereas external_dist is used to set the lengths of 
        edges associated with tip nodes. If neither internal_dist nor 
        external_dist are specified, every edge in the tree will have length 
        1.0. If only internal_dist is specified, external_dist will be made 
        equal to internal_dist. Likewise, if only external_dist is specified,
        then internal_dist will be made equal to external_dist.
        
        >>> from phycas import *
        >>> r = ProbDist.Lot(13579)
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.equiprobTree(5,r)
        >>> print t.walkPreorder()
        (0) -> [7] -> [5] -> [6] -> (1) -> (3) -> (2) -> (4)
        >>> print t.makeNewick()
        (1:1.00000,((2:1.00000,4:1.00000):1.00000,3:1.00000):1.00000,5:1.00000)
        
        In the example above, the edge lengths are all 1.0 because no edge
        length distribution was supplied to the equiprobTree function. Below
        is the same example except that both the internal and external edge 
        length distributions have been defined to be Gamma distributions with 
        mean 0.1 and variance 0.05.
        
        >>> from phycas import *
        >>> r = ProbDist.Lot(94593)
        >>> d = ProbDist.Gamma(1.0, 1.0)
        >>> d.setMeanAndVariance(0.1, 0.001)
        >>> print 'mean = %.5f, variance = %.5f' % (d.getMean(), d.getVar())
        mean = 0.10000, variance = 0.00100
        >>> d.setLot(r)
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.equiprobTree(5, r, d)
        >>> print t.makeNewick()
        (1:0.06106,(2:0.06113,3:0.03849):0.05451,(4:0.13251,5:0.08943):0.11240)
        
        """
        return TreeManipBase.equiprobTree(self, num_tips, rng, internal_dist, external_dist)
        
    def equiprobTreeUsingTreeLengthPrior(self, num_tips, rng, tree_length_dist = None):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a random tree from a discrete uniform distribution. For 
        example, for 6 taxa there are 105 distinct unrooted tree topologies,
        and if this function were to be called many times, it would generate
        each of these 105 tree topologies with probability 1/105 = 0.0095238.
        The probability distribution tree_length_dist is used to set edge 
        lengths. This should be an object of type TreeLengthDist, which
        implements the joint tree length/edge length prior introduced by 
        Rannala, Zhu and Yang (2011. MBE 29:325-335).
                
        """
        return TreeManipBase.equiprobTreeUsingTreeLengthPrior(self, num_tips, rng, tree_length_dist)
        
    def yuleTree(self, num_tips, rng, lambd):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Yule tree having num_tips tips, speciation rate lambd, and
        a topology and edge lengths determined by the supplied random number
        generator rng. num_tips should include the tip used to root the
        resulting tree.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> rng = ProbDist.Lot(13579)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.yuleTree(5, rng, 1.5)
        >>> print t.walkPreorder()
        (0) -> [5] -> [6] -> (1) -> (3) -> [7] -> (2) -> (4)
        >>> print t.makeNewick()
        (1:0.10962,(2:0.29077,4:0.29077):0.00340,(3:0.22779,5:0.22779):0.06638)

        """
        import phycas.ProbDist
        dist = phycas.ProbDist.Exponential(lambd)
        dist.setLot(rng)
        TreeManipBase.randomTree(self, num_tips, rng, dist, True)

    def randTree(self, num_tips, rng, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a random tree topology having num_tips tips, with edge lengths
        determined by the supplied probability distribution edge_len_dist.
        num_tips should include the tip used to root the resulting tree.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> rng = ProbDist.Lot(13579)
        >>> dist = ProbDist.Exponential(0.5)
        >>> dist.setLot(rng)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.randTree(5, rng, dist)
        >>> print t.walkPreorder()
        (0) -> [5] -> [6] -> [7] -> (1) -> (4) -> (3) -> (2)
        >>> print t.makeNewick()
        (1:0.90722,((2:0.22471,5:0.43341):0.36279,4:3.41679):0.56683,3:0.02039)

        """
        TreeManipBase.randomTree(self, num_tips, rng, edge_len_dist, False)

    def setRandomEdgeLengths(self, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets edge lengths using random draws from the probability distribution
        edge_len_dist.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> t.buildFromString('(1,(2,5),(3,4))')
        >>> rng = ProbDist.Lot(13579)
        >>> dist = ProbDist.Exponential(0.5)
        >>> dist.setLot(rng)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.setRandomEdgeLengths(dist)
        >>> print t.walkPreorder()
        1 -> [1002] -> [1000] -> 2 -> 5 -> [1001] -> 3 -> 4
        >>> print t.makeNewick()
        (1:0.32887,(2:0.02039,5:1.18233):0.90722,(3:3.41679,4:0.72478):0.56683)

        """
        TreeManipBase.setRandomEdgeLens(self, edge_len_dist)
        
    def setRandomInternalExternalEdgeLengths(self, internal_dist, external_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets edge lengths using random draws from the probability distribution
        internal_dist (for internal edges) or external_dist (for edges 
        subtending tips).
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> t.buildFromString('(1,(2,5),(3,4))')
        >>> rng = ProbDist.Lot(13579)
        >>> edist = ProbDist.Exponential(1.0)
        >>> edist.setLot(rng)
        >>> idist = ProbDist.Exponential(10.0)
        >>> idist.setLot(rng)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.setRandomInternalExternalEdgeLengths(idist, edist)
        >>> print t.walkPreorder()
        1 -> [1002] -> [1000] -> 2 -> 5 -> [1001] -> 3 -> 4
        >>> print t.makeNewick()
        (1:0.01644,(2:0.01019,5:0.59116):0.04536,(3:1.70839,4:0.36239):0.02834)

        """
        TreeManipBase.setRandomInternalExternalEdgeLens(self, internal_dist, external_dist)

    def deleteRandomInternalEdge(self, rng):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Deletes a random internal edge using the supplied random number 
        generator object. 

        """
        TreeManipBase.deleteRandomInternalEdge(self, rng)

    def addRightChild(self, parent_node, new_rightmost_child):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds new_rightmost_child to parent_node. 

        """
        return TreeManipBase.addRightChild(self, parent_node, new_rightmost_child)


