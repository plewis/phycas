from phycas import P
from phycas.Utilities.io import DataSource, TreeCollection, _readFileSanityCheck
from phycas.ReadNexus import FileFormats
class PhyloPackage(object):
    """A container for holding characters and trees that refer to the same set
    of taxa. The object has three attributes of interest:
    taxon_labels (a list of taxon labels),
    characters (a DataSource object)
    trees (a TreeCollection object)
    """
    def __init__(self):
        self.taxon_labels = []
        self.characters = DataSource()
        self.trees = TreeCollection()

def readFile(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a (list of taxon labels, DataSource, TreeCollection) from the file `filepath`
    
    Currently only supports NEXUS and only returns the last data matrix, but
    this will be generalized to read other formats and return the 
    supermatrix of all data matrices in the file.
    
    A side effect of calling this function is the resetting of the 
        taxon_labels, 
        characters, and
        trees
    attributes of the global object P
    """
    global P
    x = readFileNoSideEffects(filepath, format, out)
    taxon_labels, characters, trees = x.taxon_labels, x.characters, x.trees
    P.set_taxon_labels(x.taxon_labels)
    P.set_characters(x.characters)
    P.set_trees(x.trees)
    return x

def readFileNoSideEffects(filepath, format=FileFormats.NEXUS, out=None):
    """Returns a (list of taxon labels, DataSource, TreeCollection) from the file `filepath`"""
    _readFileSanityCheck(filepath, format, out)
    from phycas.ReadNexus import NexusReader
    reader = NexusReader()
    reader.readFile(filepath)
    x = PhyloPackage()
    x.taxon_labels = reader.taxa
    matrix = reader.getLastDiscreteMatrix(True)
    x.characters = DataSource(matrix=matrix, taxon_labels=x.taxon_labels)
    x.trees = TreeCollection(trees=reader.getTrees(), taxon_labels=x.taxon_labels)
    return x
