import os, sys, subprocess
from phycas.readnexus._NexusReader import FileFormats
from phycas.utilities.CommonFunctions import getDefaultOutputter
from phycas import phylogeny, Newick
_phycas_dir = None

def getPhycasDir():
    "Returns the absolute path to the directory that is the top of the phycas package"
    global _phycas_dir
    if _phycas_dir is None:
        _phycas_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    return _phycas_dir

def getPhycasTestDir():
    return os.path.join(getPhycasDir(), "tests")

def getPhycasTestDataDir():
    return os.path.join(getPhycasTestDir(), "data")

def getPhycasTestData(filen):
    """Takes a string `filen` that represents the name of one of the files in
    Phycas' test suite.
    The function returns the full path to the test file."""
    return os.path.join(getPhycasTestDataDir(), filen)

def _runRegressionTests(out):
    d = getPhycasTestDir()
    r = os.path.join(d, "runall.py")
    spawnPython(r)
    r = os.path.join(d, "doctestall.py")
    spawnPython(r)
    out.info("\nAll tests passed.")

def spawnPython(f):
    if not os.path.exists(f):
        raise RuntimeError('The python script "%s" does not exist' % f)
    retcode = subprocess.call([sys.executable, f])
    if retcode < 0:
        raise RuntimeError("python execution of %s was terminated by signal %s" % (f, str(-retcode)))
    elif retcode > 0:
        raise RuntimeError("python execution of %s was failed with code %s" % (f, str(retcode)))

def _readFileSanityCheck(filepath, format=FileFormats.NEXUS, out=None):
    if not format == FileFormats.NEXUS:
        if out is None:
            out = getDefaultOutputter()
        out.phycassert(format == FileFormats.NEXUS, "Currently only the NEXUS format is supported")
    if not os.path.exists(filepath):
        raise ValueError('The file "%s" does not exist' % filepath)

def readData(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a data matrix (or None if there is no data) from `filepath`

    Currently only supports NEXUS and only returns the last data matrix, but
    this will be generalized to read other formats and return the
    supermatrix of all data matrices in the file."""
    _readFileSanityCheck(filepath, format, out)
    from phycas.readnexus import NexusReader
    reader = NexusReader()
    reader.readFile(filepath)
    return reader.getLastDiscreteMatrix(True)

def readTrees(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a (taxon_list, list of TreeDescription objects) (or an empty list if there are no trees) from `filepath`

    Currently only supports NEXUS and only returns the last data matrix, but
    this will be generalized to read other formats and return the
    supermatrix of all data matrices in the file."""

    _readFileSanityCheck(filepath, format, out)
    from phycas.readnexus import NexusReader
    reader = NexusReader()
    reader.readFile(filepath)
    return reader.taxa, reader.getTrees()

class TreeCollection(object):

    def __init__(self, **kwargs):
        self.init(**kwargs)

    def init(self, **kwargs):
        """Valid keyword arguments include filename, format, newick, taxon_labels, trees and title.
        filename should be a string and should represent the path to a tree file.
        format is relevant only if filename is specified. The values for format are
        in the FileFormats class (but currently include only FileFormats.NEXUS).
        newick is a newick string representation of a tree topology (with or without
        edge lengths).
        taxon_labels provide names for the numbers used in a newick description (alternatively,
        taxon labels can be provided directly in the newick string)
        trees is an iterable collection of phylogeny.Tree objects.
        title is a label for the collection.
        """
        self.reset()
        self.title = kwargs.get("title")
        self.filename = kwargs.get("filename")
        if self.filename is None:
            self.filename = kwargs.get("file")
        if self.filename:
            self.format = kwargs.get("format", FileFormats.NEXUS)
            self.trees = []
            self.taxon_labels = None
        else:
            self.format = None
            self.taxon_labels = kwargs.get("taxon_labels")
            newick = kwargs.get("newick")
            if newick:
                if not isinstance(newick, Newick):
                    newick = Newick(newick)
                tree = phylogeny.Tree()
                # by default the tree is read as 1-based (as in NEXUS trees)
                newick.buildTree(tree)
                self.trees = [tree]
            else:
                self.trees = kwargs.get("trees", [])
        self.active_taxa = self.taxon_labels

    def reset(self):
        self.title = None
        self.filename = None
        self.format = None
        self.trees = []
        self.taxon_labels = None
        self.active_taxa = None
        self._needToRelabel = False

    def __iter__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns an iterator over the trees (this will trigger the reading of
        the input file if the object was initialized from a string.
        Stores tree descriptions in self.stored_tree_defs and taxon labels in
        self.taxon_labels.
        """
        if not self.trees:
            if not self.filename:
                return iter([])
            self.taxon_labels, tree_descriptions = readTrees(self.filename, self.format)
            # POL/DLS: Originally above was just:
            #            self.taxon_labels, self.trees = readTrees(self.filename, self.format)
            #          but readTrees really returns tree descriptions.  We added the loop below to convert the
            #          tree description list into a list of trees, which seems to be the intended purpose of self.trees
            self.trees = []
            for td in tree_descriptions:
                n = td.getNewickFromC()
                t = phylogeny.Tree(newick=n, rooted=td.rooted)
                self.trees.append(t)
            self._checkActiveTaxa()
        if self._needToRelabel:
            print "before rectify: ", self.trees[0].debugListTree()
            self.trees = [i.rectifyNumbers(self.active_taxa) for i in self.trees]
            print "after rectify: ", self.trees[0].debugListTree()
            self.taxon_labels = self.active_taxa
            self._needToRelabel = False
        return iter(self.trees)

    def __str__(self):
        if self.title:
            return self.title
        s, t = "", ""
        d = self.trees
        if d:
            s = "Collection of %d trees in memory. " % len(self.trees)
        if self.filename:
            t =  "Trees from the file %s" % repr(self.filename)
        elif d:
            t = "\nThe id of python Tree Collection %d" % id(self.trees)

        if s or t:
            return s + t
        return "None"

    def __deepcopy__(self, memo):
        #trees are expensive, so we don't make a deepcopy
        c = memo.get(self)
        if c is not None:
            return c
        c = TreeCollection()
        self._copyInternals(c)
        memo[self] = c
        return c

    def _copyInternals(self, c):
        c.init(**self.__dict__)

    def writeTree(self, tree, name="", rooted=None):
        if not self.trees:
            self.trees = []
        self.trees.append(tree)

    def finish(self):
        pass

    def setActiveTaxonLabels(self, tl):
        if tl:
            self.active_taxa = [i.lower() for i in tl]
        else:
            self.active_taxa = None
        self._checkActiveTaxa()

    def _checkActiveTaxa(self):
        al = self.active_taxa
        if al is None:
            self._needToRelabel = False
            return True
        if self.taxon_labels:
            tl = [i.lower() for i in self.taxon_labels]
            if al != tl:
                for a in al:
                    if a not in tl:
                        raise ValueError("Taxon %s is not found in the TreeCollection's taxon_labels (%s)" % (a, "\n".join(tl)))
                self._needToRelabel = True
            else:
                self._needToRelabel = False
        return True


class DataSource(object):
    def __init__(self, **kwargs):
        """`arg` can be a string or an iterator containing trees.
        If `arg` is a string, then the `format` keyword argument can be used to
        specify the file format (the default is NEXUS).
        If the trees are passed in as a list then a `taxon_labels` keyword
        argument is expected.
        """
        self.title = kwargs.get("title")
        self.filename = kwargs.get("filename")
        self.data_obj = kwargs.get("matrix") or kwargs.get("data_obj")
        self.taxon_labels = kwargs.get("taxon_labels")
        self.format = kwargs.get("format", FileFormats.NEXUS)

    def _reset(self):
        self.__init__()

    def _setEqualTo(self, other):
        "Copy the state from other (currently other must be a DataSource object)"
        self._reset()
        self.__init__(**other.__dict__)

    def getMatrix(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a DataMatrix object (or None if the DataSource is empty)"
        """
        if self.data_obj is None:
            if not self.filename:
                return None
            self.data_obj = readData(self.filename, self.format)
        return self.data_obj

    def __str__(self):
        if self.title:
            return self.title
        s, t = "", ""
        d = self.data_obj
        if d:
            s = "Character matrix in memory (id %d) with %d characters for %d taxa.  " % (id(d), d.n_char, d.n_tax)
        if self.filename:
            t =  "Characters from the file %s" % self.filename
        if s or t:
            return s + t
        return "None"

    def __deepcopy__(self, memo):
        #trees are expensive, so we don't make a deepcopy
        if self.data_obj:
            return DataSource(matrix=self.data_obj, taxon_labels=self.taxon_labels, title=self.title)
        elif self.filename:
            return DataSource(filename=self.filename, format=self.format, title=self.title)
        return DataSource(taxon_labels=self.taxon_labels, title=self.title)

