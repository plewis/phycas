from phycas import randomtree
from phycas.utilities.PhycasCommand import *
from phycas.commands.SumTImpl import TreeSummarizer
from phycas.commands.RefDistImpl import RefDistImpl
import copy

class RefDist(PhycasCommand):
    def __init__(self):
        args = (    ("trees",        None,   "A source of trees (list of trees or the name of the input tree file) to be summarized. Required for multi-tree-topology generalized stepping-stone analysis", TreeSourceValidate),
                    ("params",       None,   "Name of file containing sampled parameter values. Required.", FileExistsValidate),
                    ("skip",            1,      "Number of trees from the input list of trees to skip", IntArgValidate(min=0)),
                    ("epsilon",      0.01,   "Fraction of total sample size used to ensure 0 < split posterior < 1", ProbArgValidate()),
                    ("rooted",      False,  "Set to True if trees are rooted; otherwise, leave set to default value of False to assume trees are unrooted", BoolArgValidate),
                )
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["refdistfile"]

        t = TextOutputSpec(prefix='refdistfile', suffix=".txt", help_str="The reference distribution definition file that can be used as input to the stepping-stone (ss) command.")
        o.__dict__["refdistfile"] = t

        PhycasCommand.__init__(self, args, "refdist", "The refdist command is used to summarize a posterior sample of parameters and trees and use that summary to define a reference distribution that approximates the posterior distribution. This command should be used after an MCMC analysis and before a stepping-stone analysis.", o)

    def hidden():
        """
        Overrides the PhycasCommand.hidden method to keep SumT's name from being displayed
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return False

    hidden = staticmethod(hidden)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        c = copy.deepcopy(self)
        ref_dist = RefDistImpl(c)
        ref_dist.run()

