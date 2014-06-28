from phycas import P
from phycas.utilities.PhycasCommand import *
from phycas.commands.RandomTreeImpl import TreeSimulator
from phycas.probdist import Exponential
import copy
class RandomTree(PhycasCommand):
    def __init__(self):
        args = tuple(
                    PhycasCommand._getRNGOptions() +
                   [("distribution", "Yule", "'Yule' or 'Equiprobable' -- The name of the tree generation process or distribution from which the trees will be drawn", EnumArgValidate(["yule", "equiprobable"])),
                    ("taxon_labels", P.taxon_labels , "The names of the taxa to simulate"),
                    ("edgelen_dist",  Exponential(10.0), "Used to generate edge lengths. This can be None if distribution = 'yule'; in this case the branch lengths from the Yule process will be used"),
                    ("speciation_rate",  None, "The rate of speciation that governs the branch lengths of the Yule tree simulation process. This is only used if edgelen_dist is None", FloatArgValidate(min=1.0e-10)),
                    ("n_trees", 0, "The number of trees to generate (if 0 is specified, a bottomless collection of trees is generated)", IntArgValidate(min=0)),
                    ("n_taxa", 0, "The number of taxa to generate (only used if the taxon_labels attribute is not specified)", IntArgValidate(min=0)),
                    ("newick", None, "Tree topology to simulate branch lengths on."),
                   ]
                   )
        o = PhycasCommandOutputOptions()
        PhycasCommand.__init__(self, args, "randomtree", "Produces a TreeCollection by simulating trees using the Yule Process", o)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        c = copy.deepcopy(self)
        return TreeSimulator(c)

