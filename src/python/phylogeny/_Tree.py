from _PhylogenyExt import *
from phycas.readnexus._NexusReader import TreeDescription

class Tree(TreeBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates the notion of a phylogenetic tree. This class has only
    methods necessary for representing the tree, copying trees from other
    trees, manipulating the tree, and representing the tree graphically.
    It does not perform any specialized activities, such as computing its
    own likelihood: these sorts of activities are left to member
    functions of other classes.

    """
    def __init__(self, builder=None, newick=None, zero_based=False, taxon_labels=None, rooted=False):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls clear() to initialize data members.
        If a `builder` is passed in, then the builder.buildTree() function
        will be called to initialize the tree.
        If no builder is passed in, but `newick` is specified then the tree 
        definition comes from the newick string and `zero_based` is used to 
        indicate whether the taxa indices start at 0 or 1.

        >>> from phycas.phylogeny import *
        >>> t1 = Tree()
        >>> print t1.getNNodes()
        0

        """
        TreeBase.__init__(self)
        self.setRootedness(rooted)
        if builder:
            builder.buildTree(self)
        elif newick:
            self.buildFromString(newick, zero_based)
        self.taxon_labels = taxon_labels
        if self.getFirstPreorder():
            if self.taxon_labels:
                self.rectifyNumbers(self.taxon_labels)
            else:
                self.taxon_labels = []
                someLabels = False
                for tip in self.getFullTipList():
                    n = tip.getNodeName()
                    if not n:
                        n = None
                    else:
                        someLabels = True
                    self.taxon_labels.append(n)
                if not someLabels:
                    self.taxon_labels = []

    def buildTree(self, tree=None):
        """Calls buildFromString or other appropriate method to construct a tree 
        in the the variable `tree`
        Returns the `tree` instance."""
        if tree is None:
            return phylogeny.Tree(newick=self.newick, taxon_labels=list(self.taxon_labels))
        # makeNumberedNewick returns a one-based number
        tree.buildFromString(self.makeNumberedNewick(), False) 
        return tree
        
    def __str__(self):
        return self.newick
        
    def __iter__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes in preorder sequence. The tip node named "a" appears
        first in the example below because the tree is rooted at the first
        tip node for "unrooted" trees.

        >>> from phycas.phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree:
        ...     print nd.getNodeName(),
        ... 
        a z b y c x d e

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            yield nd
            nd = nd.getNextPreorder()
        
    def nodesWithEdges(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes that have edges, which is all nodes in an unrooted tree
        except the root node. The nodes are visited in preorder sequence.

        >>> from phycas.phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree.nodesWithEdges():
        ...     print nd.getNodeName(),
        ... 
        z b y c x d e

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            if not nd.isRoot():
                yield nd
            nd = nd.getNextPreorder()
        
    def iterTipNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes that correspond to tip nodes, where a tip node is any
        node (including the root node) that is of degree one (i.e. has only
        one edge attached to it). The nodes are visited in preorder sequence.

        >>> from phycas.phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree.iterTipNodes():
        ...     print nd.getNodeName(),
        ... 
        a b c d e

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            if nd.isTip():
                yield nd
            nd = nd.getNextPreorder()
        
    def iterInternalNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes that correspond to internal nodes, where an internal
        node is any node that has degree greater than two (i.e. has at least
        three edges attached to it). The nodes are visited in preorder
        sequence.

        >>> from phycas.phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree.iterInternalNodes():
        ...     print nd.getNodeName(),
        ... 
        z y x

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            if nd.isInternal():
                yield nd
            nd = nd.getNextPreorder()
        
    def buildFromString(self, newick, zero_based_tips = False):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Builds a tree from a newick tree description. If edge lengths are
        specified in the tree description, they must be specified for all
        edges. No checking of node names is performed as the tree is created.
        Before returning, roots tree at first tip node. Raises an XPhylogeny
        exception if a problem is encountered.

        >>> from phycas.phylogeny import *
        >>> t2 = Tree()
        >>> t2.buildFromString('(c,d,(a, b))')
        >>> print t2.walkPreorder()
        c -> [1001] -> d -> [1000] -> a -> b

        """
        TreeBase.buildFromString(self, newick, zero_based_tips)

    def clear(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns object to just-constructed state.

        >>> from phycas.phylogeny import *
        >>> t3 = Tree()
        >>> t3.buildFromString('(c,d,(a, b))')
        >>> t3.getNNodes()
        6
        >>> t3.clear()
        >>> t3.getNNodes()
        0

        """
        TreeBase.clear(self)

    def findTipByName(self, tipname):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the node number of a tip node named tipname, if one can be
        found. Otherwise, returns None.

        """
        n = TreeBase.findTipByName(self, tipname)
        if n < self.getNTips():
            return n
        else:
            return None

    def getNInternals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of internal nodes currently composing the tree
        (i.e. getNInternals(). The number of internal nodes equals the number
        of nodes in the tree having degree 2 or higher. Its value is 1 for a
        (rooted or unrooted) star tree, n-2 for a fully-resolved unrooted
        tree, and n-1 for a fully-resolved rooted tree. Calls
        refreshNodeCounts() if node counts have been invalidated.

        >>> from phycas.phylogeny import *
        >>> t4 = Tree()
        >>> t4.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t4.getNNodes()
        6
        >>> print t4.getNInternals()
        2

        """
        return TreeBase.getNInternals(self)

    def getNTips(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of leaf nodes currently composing the tree (i.e.
        getNLeaves(). The number of leaves equals the number of degree-1 nodes
        in the tree. Note that one of these degree-1 nodes is the root node.
        The number of leaves equals the number of taxa for unrooted trees.
        For rooted trees of n taxa, there are n+1 leaves because the root node
        is a leaf node but not a taxon. Calls refreshNodeCounts() if node
        counts have been invalidated.

        >>> from phycas.phylogeny import *
        >>> t5 = Tree()
        >>> t5.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t5.getNNodes()
        6
        >>> print t5.getNTips()
        4

        """
        return TreeBase.getNTips(self)

    def getNNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the total number of nodes currently composing the tree (i.e.
        getNLeaves() plus getNInternals(). For unrooted trees of n taxa, one
        of the tip nodes serves as the root node, so this number can range
        from n+1 (star tree) to 2n-2 (fully-resolved). For rooted trees, the
        range is n+2 (star tree) to 2n-1 (fully-resolved) because the root
        node exists but is not identical to one of the tips. Calls
        refreshNodeCounts() if node counts have been invalidated.

        >>> from phycas.phylogeny import *
        >>> t6 = Tree()
        >>> t6.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t6.getNNodes()
        6

        """
        return TreeBase.getNNodes(self)

    def isRooted(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if tree is rooted, or False if the root is actually a
        tip.

        >>> from phycas.phylogeny import *
        >>> t7 = Tree()
        >>> t7.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t7.isRooted()
        False

        """
        return TreeBase.isRooted(self)

    #def hasEdgeLens(self):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Returns True if edge lengths have been specified for the tree, False
    #    otherwise.
    #
    #    """
    #    return TreeBase.hasEdgeLens(self)

    def hasEdgeLens(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if tree has edge lengths defined, False if edge lengths
        have not yet been set.

        >>> from phycas.phylogeny import *
        >>> t8 = Tree()
        >>> t8.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t8.hasEdgeLens()
        False
        >>> t8.buildFromString('(fish:0.2,shark:0.25,(bird:0.1, mammal:0.15):0.09)')
        >>> print t8.hasEdgeLens()
        True
        
        """
        return TreeBase.hasEdgeLens(self)

    def walkPreorder(self, verbosity = 0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Walks through the tree in preorder fashion (visits parents before
        descendants, and leftmost child before siblings), building up a string
        showing the path taken. Each unnamed internal node is represented by
        the node number (0..getNInternals() - 1) in square brackets, and each
        unnamed leaf node is represented by the node number (0..getNLeaves()
        - 1) in parentheses. If a name has been provided for a node, that
        name is used instead of the node number, both for internal and leaf
        nodes.

        >>> from phycas.phylogeny import *
        >>> t9 = Tree()
        >>> t9.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t9.walkPreorder()
        fish -> [1001] -> shark -> [1000] -> bird -> mammal

        """
        return TreeBase.debugWalkTree(self, True, verbosity)

    def walkPostorder(self, verbosity = 0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Walks through the tree in postorder fashion (visiting children before
        parents and siblings on right before those on the left, building up a
        string showing the path taken. Each unnamed internal node is represented by
        the node number (0..getNInternals() - 1) in square brackets, and each
        unnamed leaf node is represented by the node number (0..getNLeaves()
        - 1) in parentheses. If a name has been provided for a node, that
        name is used instead of the node number, both for internal and leaf
        nodes.

        >>> from phycas.phylogeny import *
        >>> tA = Tree()
        >>> tA.buildFromString('(fish,shark,(bird, mammal))')
        >>> print tA.walkPostorder()
        mammal -> bird -> [1000] -> shark -> [1001] -> fish

        """
        return TreeBase.debugWalkTree(self, False, verbosity)

    def setRootNode(self, root_tip):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets root of tree to supplied root_tip, which is assumed to be a tip.
        """
        TreeBase.rerootAtThisTip(self, root_tip)
    
    def rerootAtTip(self, num):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reroots the tree at the leaf node numbered num. An XPhylogeny
        exception is raised if a leaf node having number num cannot be found.

        >>> from phycas.phylogeny import *
        >>> tC = Tree()
        >>> tC.buildFromString('((a,b),c,(d,e))')
        >>> print tC.walkPreorder()
        a -> [1000] -> b -> [1002] -> c -> [1001] -> d -> e
        >>> tC.rerootAtTip(5)
        Traceback (most recent call last):
            ...
        Exception: there is no tip node having number 5
        >>> tC.rerootAtTip(4)
        >>> print tC.walkPreorder()
        e -> [1001] -> d -> [1002] -> c -> [1000] -> b -> a

        """
        return TreeBase.rerootAtTip(self, num)

    def edgeLenSum(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sums all edge lengths in the tree. Raises an Exception if edge lengths
        were never provided for the tree.

        >>> from phycas.phylogeny import *
        >>> tB = Tree()
        >>> tB.buildFromString('(fish,shark,(bird, mammal))')
        >>> print tB.edgeLenSum()
        Traceback (most recent call last):
            ...
        Exception: no edge lengths were specified for this tree
        >>> tB.buildFromString('(fish:0.2,shark:0.25,(bird:0.1, mammal:0.15):0.09)')
        >>> print round(tB.edgeLenSum(), 3)
        0.79

        """
        return TreeBase.edgeLenSum(self)

    def edgeLens(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list of all edge lengths.
        
        Need an example here.

        """
        return TreeBase.edgeLens(self)

    def makeNewick(self, float_precision = 5):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string representing the tree in newick format (i.e. using
        nested parentheses). If a node is named, the name will be used in the
        tree description. Tip nodes are represented by node numbers if they
        do not have a name. If edge lengths are present, they will be shown.

        >>> from phycas.phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a:0.11,b:0.12)x:0.10,c:0.21,(d:0.31,e:0.32)y:0.30)z')
        >>> print t.makeNewick()
        (a:0.11000,b:0.12000,(c:0.21000,(d:0.31000,e:0.32000)y:0.30000)z:0.10000)x

        """
        return TreeBase.makeNewick(self, float_precision)

    def makeNumberedNewick(self, float_precision = 5):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a 1-based newick representing the tree in newick format (i.e. 
        using nested parentheses). If a node is named, the name will be used in the
        tree description. Tip nodes are represented by node numbers if they
        do not have a name. If edge lengths are present, they will be shown.

        >>> from phycas.phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a:0.11,b:0.12)x:0.10,c:0.21,(d:0.31,e:0.32)y:0.30)z')
        >>> print t.makeNumberedNewick()
        (1:0.11000,2:0.12000,(3:0.21000,(4:0.31000,5:0.32000):0.30000):0.10000)

        """
        return TreeBase.makeNumberedNewick(self, float_precision)
        
    def rectifyNumbers(self, taxon_names):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets tip node numbers in the tree to the index of the tip's name in
        the supplied object taxon_names, which should be a list or tuple of
        taxon names.

        >>> from phycas.phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a,b),c,(d,e))')
        >>> print t.walkPreorder(verbosity=1)
        a (0) -> ? [1000] -> b (1) -> ? [1002] -> c (2) -> ? [1001] -> d (3) -> e (4)
        >>> t.rectifyNumbers(['e','d','c','b','a'])
        >>> print t.walkPreorder(verbosity=1)
        a (4) -> ? [1000] -> b (3) -> ? [1002] -> c (2) -> ? [1001] -> d (1) -> e (0)

        """
        self.taxon_labels = taxon_names
        return TreeBase.rectifyNumbers(self, taxon_names)

    def rectifyNames(self, taxon_names):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets tip node names in the tree to the names in the supplied object
        taxon_names, which should be a list or tuple of taxon names. The node
        numbers are treated as indices into the taxon_names list.

        >>> from phycas.phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a,b),c,(d,e))')
        >>> print t.walkPreorder(verbosity=1)
        a (0) -> ? [1000] -> b (1) -> ? [1002] -> c (2) -> ? [1001] -> d (3) -> e (4)
        >>> t.rectifyNames(['e','d','c','b','a'])
        >>> print t.walkPreorder(verbosity=1)
        e (0) -> ? [1000] -> d (1) -> ? [1002] -> c (2) -> ? [1001] -> b (3) -> a (4)

        """
        self.taxon_labels = taxon_names
        return TreeBase.rectifyNames(self, taxon_names)

    def unselectAllNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Turns the selected attribute off for all nodes in the tree.

        """
        return TreeBase.unselectAllNodes(self)

    def getFirstPreorder(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the root node. Use the TreeNode getNextPreorder to walk
        through the other nodes in the tree in preorder fashion.

        """
        return TreeBase.getFirstPreorder(self)

    def getFirstPostorder(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the node at the upper right corner of the tree (the last
        node in the preorder sequence, which is also the first node in the
        postorder sequence). Use the TreeNode getNextPostorder to walk
        through the other nodes in the tree in postorder fashion.

        """
        return TreeBase.getFirstPostorder(self)

    def ladderizeRight(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Rearranges tree so that for each internal node, the child with the
        most descendants is on the right.

        """
        return TreeBase.ladderize(self, True)

    def ladderizeLeft(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Rearranges tree so that for each internal node, the child with the
        most descendants is on the left.

        """
        return TreeBase.ladderize(self, False)
        
    def getTips(self):
        """Returns the tips in preorder traversal order"""
        return [i for i in self.iterTipNodes()]

    def getFullTipList(self):
        """Returns the list of nodes or None up to the index corresponding to 
        the highest node index.
        """
        d = {}
        max_n = 0
        for nd in self.iterTipNodes():
            n = nd.getNodeNumber()
            d[n] = nd
            max_n = max(max_n, n)
        l = [None] * (max_n + 1)
        for k, v in d.iteritems():
            l[k] = v
        return l

    rooted = property(isRooted)
    newick = property(makeNewick)
