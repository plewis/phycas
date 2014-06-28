import copy, os, sys, math, random, copy
from phycas.utilities.PhycasCommand import str_value_for_user
from phycas.utilities.CommonFunctions import CommonFunctions
from phycas.utilities.io import TreeCollection
from phycas import phylogeny, Newick

class TreeSimulator(CommonFunctions, TreeCollection):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    A collection of simulated trees.

    The n_trees attribute controls the size of the collection.  If it is set to
    0 then the collection will be bottomless and trees will not be stored (thus
    indexing the collection is not supported).

    If n_trees is > 0, then the trees will be simulated on demand, but stored
    so that the same tree will be returned with subsequent indexing with the same
    number.
    """

    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes TreeSummarizer object by assigning supplied phycas object
        to a data member variable.

        """
        CommonFunctions.__init__(self, opts)
        self.r = self._getLot()
        eld = self.edgelen_dist
        self._current = 0
        self.trees = []
        if eld is None:
            if self.distribution.lower() != "yule":
                raise ValueError("edgelen_dist must be set if the tree simulator's distribution is not Yule")
            if self.speciation_rate is None:
                raise ValueError("edgelen_dist or the speciation_rate must be set to generate trees from the Yule process")
        else:
            eld.setLot(self.r)
        self.stdout.debugging("Generating %s" % str(self))

    def __str__(self):
        nt = self.n_trees
        n = nt > 0 and "%d " % nt or ""
        eld = self.edgelen_dist
        if self.distribution.lower() != "yule":
            t = "from the equiprobable distribution"
        else:
            t = "from the Yule process"
        if eld is None:
            y = "with speciation rate = %f" % self.speciation_rate
        else:
            y = "with edge lengths drawn from %s" % str_value_for_user(eld)
        return "collection of %stree(s) simulated by %s %s" % (n, t, y)

    def setActiveTaxonLabels(self, tl):
        if not self.taxon_labels:
            self.taxon_labels = tl
        if tl:
            self.active_taxa = [i.lower() for i in tl]
        else:
            self.active_taxa = None
        self._checkActiveTaxa()

    def writeTree(self, tree, name="", rooted=None):
        raise ValueError("TreeSimulator collection of trees cannot be used as an output destination")

    def finish(self):
        raise ValueError("TreeSimulator collection of trees cannot be used as an output destination")

    def __iter__(self):
        self._current = 0
        return self

    def _simulateTree(self):
        if (not self.taxon_labels) and self.n_taxa < 1 and (not self.newick):
            raise ValueError("TreeSimulator cannot be created with an empty list of taxon_labels and n_taxa set to 0")
        tl = self.taxon_labels
        if tl:
            t = phylogeny.Tree(taxon_labels=self.taxon_labels)
            ntax = len(self.taxon_labels)
        else:
            t = phylogeny.Tree()
            ntax = self.n_taxa
        tm = phylogeny.TreeManip(t)
        eld = self.edgelen_dist
        self.stdout.debugging("RandomTree._simulateTree seed = %d" %self.r.getSeed())
        n = self.newick
        if n:
            if not isinstance(n, Newick):
                n = Newick(n)
            n.buildTree(t)
            tm.setRandomEdgeLengths(self.edgelen_dist)
        elif self.distribution.lower() == "yule":
            if eld is None:
                isp = self.speciation_rate
                if isp is None or isp <= 0.0:
                    raise ValueError("edgelen_dist or the speciation_rate must be set to generate trees from the Yule process")
                eld = Exponential(isp)
                eld.setLot(self.r)
                self.stdout.debugging("Calling TreeManip.yuleTree")
                tm.yuleTree(ntax, self.r, eld)
            else:
                self.stdout.debugging("Calling TreeManip.randTree")
                tm.randTree(ntax, self.r, eld)
        else:
            self.stdout.debugging("Calling TreeManip.equiprobTree")
            tm.equiprobTree(ntax, self.r, eld)
        return t

    def next(self):
        n = self._current
        nt = self.n_trees
        if n >= nt and nt > 0:
            raise StopIteration()
        t = None
        if n < len(self.trees):
            t = self.trees[n]
        if t is None:
            t = self._simulateTree()
        self._current += 1
        if nt > 0:
            self._pad_to(self._current)
            self.trees[n] = t
        return t

    def __getitem__(self, k):
        nt = self.n_trees
        if nt < 1:
            raise TypeError("TreeSimulator without n_trees cannot return trees by index")
        if isinstance(k, slice):
            inds = k.indices(nt)
            l = []
            for i in range(inds[0], inds[1], inds[2]):
                l.append(self.__getitem__(i))
            return l
        if k < 0:
            k += nt
        if k > nt and nt > 0:
            raise IndexError("TreeSimulator index out of range")
        self._pad_to(k+1)
        t = self.trees[k]
        if t is None:
            t = self._simulateTree()
            self.trees[k] = t
        return t

    def _pad_to(self, k):
        diff_len = k - len(self.trees)
        if diff_len > 0:
            self.trees.extend([None]*diff_len)


    def __len__(self):
        if self.n_trees < 1:
            raise TypeError("TreeSimulator without n_trees has no len()")
        return self.n_trees

    def getNTrees(self):
        return self.opts.n_trees

    def setNTrees(self, x):
        self.opts.n_trees = x

    def getTaxa(self):
        return self.opts.taxon_labels

    def setTaxa(self, x):
        self.opts.taxon_labels = x

    def getNTaxa(self):
        return self.opts.n_taxa

    def setNTaxa(self, x):
        self.opts.n_taxa = x

    def getDistribution(self):
        return self.opts.distribution

    def setDistribution(self, x):
        self.opts.distribution = x

    def getEdgeLenDist(self):
        return self.opts.edgelen_dist

    def setEdgeLenDist(self, x):
        self.opts.edgelen_dist = x

    def getSpeciationRate(self):
        return self.opts.speciation_rate

    def setSpeciationRate(self, x):
        self.opts.speciation_rate = x

    def getNewick(self):
        return self.opts.newick

    def setNewick(self, x):
        self.opts.newick = x

    distribution = property(getDistribution, setDistribution)
    taxon_labels = property(getTaxa, setTaxa)
    edgelen_dist = property(getEdgeLenDist, setEdgeLenDist)
    speciation_rate = property(getSpeciationRate, setSpeciationRate)
    n_trees = property(getNTrees, setNTrees)
    n_taxa = property(getNTaxa, setNTaxa)
    newick = property(getNewick, setNewick)

