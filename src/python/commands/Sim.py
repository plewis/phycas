from phycas.Utilities.PhycasCommand import *
from phycas import model, randomtree
from phycas.ProbDist import Exponential
from phycas.Commands.SimImpl import SimImpl

class Sim(PhycasCommand):
    def __init__(self):
        args = tuple([
                ("file_name",       'simulated.nex',                          "Name of file in which to save simulated data'"),
                ("model",           model,                                    "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("taxon_labels",    ['taxon1', 'taxon2', 'taxon3', 'taxon4'], "Names to use for taxa in simulated data set (number of labels defined determines the number of taxa in the simulated dataset)"),
                ("nchar",           1000,                                     "Number of characters to generate", IntArgValidate(min=1)),
                ("edgelen_dist",    Exponential(10.0),                        "Distribution from which random edge lengths are drawn when edge lengths are not specified in tree_source"),
                ("tree_source",     None,                                     "TreeCollection that will provide the model tree for the simulation. Examples: tree_source = TreeCollection(filename='treefile.tre'); tree_source = TreeCollection(newick='(A:0.01,B:0.07,(C:0.03,D:0.11):0.07'). Can also be object returned by randomtree(); type randomtree.help for details", TreeSourceValidate),
                ] + PhycasCommand._getRNGOptions()
                )
        PhycasCommand.__init__(self, args, "sim", "Simulates DNA sequences according to the specified model.")

        # The data members added below should be hidden from the user because they are irrelevant to simulating data
        # They must be present, however, because they are referenced in the LikelihoodCore class, which is also used
        # for likelihood and Bayesian analyses.
        #
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["fix_edgelens"]   = False
        self.__dict__["uf_num_edges"]   = 50
        self.__dict__["data_source"]    = None

    def hidden():
        """
        Overrides the PhycasCommand.hidden method to keep Sim's name from being displayed
        in the list of classes displayed when users type help. Change the return value to
        False when it is ready to be advertised.
        """
        return True

    hidden = staticmethod(hidden)

    def __call__(self, **kwargs):
        self.set(**kwargs)
        simulate = SimImpl(self)
        simulate.run()
        return simulate

