from _ReadNexusExt import *
from phycas.DataMatrix._DataMatrix import DataMatrix
import re

_ROUND_TRIP_EVERY_NEXUS_READ = False

class FileFormats:
    NEXUS, FASTA, PHYLIP, PDF, RAW_TEXT = range(5)
    _names = ["NEXUS", "FASTA", "PHYLIP", "PDF", "RAW_TEXT"]
    def to_str(v):
        if v < 0 or v > len(FileFormats._names):
            raise VauleError("Invalid file format code (%s) specified" % str(v))
        return FileFormats._names[v]
    to_str = staticmethod(to_str)


class TreeDescription(object):
    def __init__(self, nxs_full_tree_description):
        self.ftd = nxs_full_tree_description

    def getNameFromC(self):
        return self.ftd.getName()

    def getNewickFromC(self):
        return self.ftd.getNewick()
    makeNewick = getNewickFromC

    def getIsRootedFromC(self):
        return self.ftd.isRooted()

    def buildTree(self, tree=None):
        """Calls buildFromString or other appropriate method to construct a tree 
        in the the variable `tree`
        Returns the `tree` instance."""
        if tree is None:
            tree = Phylogeny.Tree()
        tree.buildFromString(self.newick, False) # tree descriptions from NCL 2.1 are 1-based not 0-based
        return tree
    name = property(getNameFromC)
    newick = property(getNewickFromC)
    rooted = property(getIsRootedFromC)

_tokenBreakers = re.compile(r"[(){}_\"-\[\]/\\,;:=*`+<>' \t\n]")
_quoteDemandingTokenBreakers = re.compile('[-\\]\\[/\\,;:=*`+<>\'\\t\\n_(){}\\"]')
def escapeNexus(s):
    global _tokenBreakers, _quoteDemandingTokenBreakers
    ls = len(s)
    if ls > 1:
        if len(_tokenBreakers.split(s)) > 1:
            if len(_quoteDemandingTokenBreakers.split(s)) > 1:
                return "'%s'" % "''".join(s.split("'"))
            return s.replace(" ", "_")
        return s
    if ls == 0:
        return "''"
    if (s in "'_[") or (s.isspace()):
        if s == "'":
            return "''''"
        return "'%s'" % s
    return s

class NexusWriter(object):
    NO_BLOCK, TAXA_BLOCK, CHARACTERS_BLOCK, TREES_BLOCK = range(4)
    def __init__(self, outstream, preexisting=False, taxa=[]):
        self.outstream = outstream
        self._preexisting = preexisting
        self._started_writing = False
        self._in_block = NexusWriter.NO_BLOCK
        self._written = []
        self.taxa = taxa
    def _writeTaxa(self):
        assert(self.outstream)
        if (not self._preexisting) and (not self._started_writing):
            self.outstream.write("#NEXUS\n")
        self._started_writing = True
        if self.taxa:
            taxa = self.taxa
            ntax = len(self.taxa)
            labels = [escapeNexus(t) for t in taxa]
            self.outstream.write("""BEGIN TAXA;
  Dimensions ntax = %d;
  TaxLabels
    %s ;
END;
""" % (ntax, "\n    ".join(labels)))
            self._written.append(NexusWriter.TAXA_BLOCK)

    def writeTree(self, tree, name="", rooted=None):
        # tree is of type Tree
        assert(self.outstream)
        assert(self.taxa)
        if not self._in_block == NexusWriter.TREES_BLOCK:
            if not NexusWriter.TAXA_BLOCK in self._written:
                self._writeTaxa()
            taxa = self.taxa
            ntax = len(self.taxa)
            labels = [escapeNexus(t) for t in taxa]
            t_arg = ",\n    ".join(["%d %s" % (i+1, n) for i,n in enumerate(labels)])
            self.outstream.write("BEGIN TREES;\n  Translate\n    %s;\n" % t_arg)
            self._in_block = NexusWriter.TREES_BLOCK
        if rooted is None:
            rooted = tree.rooted
        c = rooted and "R" or "U"
        n = name or tree.name
        self.outstream.write("  Tree %s = [&%s] %s;\n" % (escapeNexus(n), c, tree.newick))
        self._written.append(NexusWriter.TREES_BLOCK)

    def writeCharacters(self, data_matrix):
        assert(self.outstream)
        assert(self.taxa)
        if not NexusWriter.TAXA_BLOCK in self._written:
            self._writeTaxa()
        self.finish_block()
        out.write("Begin Characters;\n Dimensions nchar = %d;\n " % data_matrix.n_char)
        out.write(data_matrix.getNEXUSFormatCommand())
        out.write("\n Matrix\n")
        assert(len(et) == dm.n_tax)
        for ind, name in enumerate(et):
            out.write("%-20s "% name)
            r = data_matrix.getRow(ind)
            for c in r:
                out.write(data_matrix.getCodeToSymbol(c))
            out.write("\n")
        out.write(";\nEnd;\n\n")
        self._written.append(NexusWriter.CHARACTERS_BLOCK)
        
    def finish_block(self):
        if not self.outstream:
            return
        if self._in_block != NexusWriter.NO_BLOCK:
            self.outstream.write("END;\n")
            self._in_block = NexusWriter.NO_BLOCK
            
    finish = finish_block

class CharacterSuperMatrix(object): # Kansas2010Todo
    def __init__(self, raw_supermatrix, taxa):
        self.raw_supermatrix = raw_supermatrix
        n = raw_supermatrix.getNumMatrices()
        self._raw_matrix_list = [raw_supermatrix.getMatrix(i) for i in range(n)]
        self.matrix_list = [DataMatrix(i, taxa) for i in self._raw_matrix_list]
    def getNChar(self):
        return sum([i.getNChar() for i in self.matrix_list])
    def getNTax(self):
        return sum([i.getNTax() for i in self.matrix_list])
    def getNSubsets(self):
        return len(self.matrix_list)
    def __getitem__(self, k):
        return self.matrix_list[k]
    def getTaxLabels(self):
        for i in self.matrix_list[1:]:
            phycassert(self.matrix_list[0].taxa == i.taxa)
        return self.matrix_list[0].taxa
    n_char = property(getNChar)
    n_tax = property(getNTax)
    n_subsets = property(getNSubsets)
    taxa = property(getTaxLabels)
    
class NexusReader(NexusReaderBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Need to write.
    
    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        >>> from phycas import *
        >>> reader = ReadNexus.NexusReader()
        >>> reader.readFile(getPhycasTestData('nyldna4.nex'))
        >>> print reader.getLastDiscreteMatrix().n_char
        3080

        """
        NexusReaderBase.__init__(self, -1)
        self.taxa = None
        self._data_matrices = None

    def readFile(self, fn):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to              write.

        """
        self.taxa = None
        self._data_matrices = None
        NexusReaderBase.readFile(self, fn)
        self.taxa = NexusReaderBase.getTaxLabels(self)
        if _ROUND_TRIP_EVERY_NEXUS_READ:
            import sys
            sys.stdout.write("==========================DEBUG OUTPUT FROM %s  ================\n" % (__file__))
            self.writeNEXUS(sys.stdout)
            sys.stdout.write("========================== END DEBUG OUTPUT FROM %s  ================\n" % (__file__))
        
        
    def getErrorMessage(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        return NexusReaderBase.getErrorMessage(self)
        
    def getTrees(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Need to write.

        """
        return [TreeDescription(t) for t in NexusReaderBase.getTrees(self)]
    
    def getLastDiscreteMatrix(self, gaps_to_missing=True):
        if self._data_matrices is None:
            o = getLastNexusDiscreteMatrix(self, gaps_to_missing)
            if o is not None:
                self._data_matrices = CharacterSuperMatrix(o, self.taxa)
                return self._data_matrices
        return self._data_matrices

    def writeNEXUS(self, out, appending=False):
        nw = NexusWriter(out, appending, self.getTaxLabels())
        dm = self.getLastDiscreteMatrix()
        if dm:
            nw.writeCharacters(dm)
        tr = self.getTrees()
        for t in tr:
            nw.writeTree(t)

