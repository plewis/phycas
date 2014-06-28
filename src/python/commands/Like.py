from phycas.utilities.PhycasCommand import *
from phycas import model, randomtree, P
from phycas.probdist import Exponential
from phycas.commands.LikeImpl import LikeImpl

class Like(PhycasCommand):
    def __init__(self):
        args = (
                 PhycasCommand._getRNGOptions() +
                [
                ("data_source",             P.characters,               "The DataSource that provides the data, to be used in the MCMC analysis. Should be a DataSource object", DataSourceValidate),
                ("model",                   model,                      "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("tree_source",             randomtree(),               "TreeCollection that will provide the tree.", TreeSourceValidate),
                ("starting_edgelen_dist",   Exponential(10.0),          "Used to select the starting edge lengths when tree_source is 'random'"),
                ("store_site_likes",         False,                      "If True, site log-likelihoods will be stored and can be retrieved using the getSiteLikes() function"),
                ("uf_num_edges",              50,    "Number of edges to traverse before taking action to prevent underflow", IntArgValidate(min=1)),
                ]
                )
        PhycasCommand.__init__(self, args, "like", "Calculates the log-likelihood under the current model.")

        # The data members added below should be hidden from the user because they are irrelevant to
        # computing the likelihood. They must be present, however, because they are referenced in the
        # LikelihoodCore class, which is also used for Bayesian analyses.
        #
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["random_seed"] = 0
        self.__dict__["fix_edgelens"]   = False
        #self.__dict__["uf_num_edges"]   = 50
        self.__dict__["data_source"]    = 'file'
        self.__dict__["site_likes"]    = []
        self.__dict__["pattern_counts"]    = []
        self.__dict__["char_to_pattern"]    = []
        self.__dict__["site_uf"]    = []

        # The following data members are hidden from user but useful for developers in debugging
        self.__dict__["preorder_edgelens"]    = None    # can provide a vector of edgelens that will be substituted before likelihood is computed

    def getSiteLikes(self):
        '''
        Returns a list of pattern log-likelihoods. Note that sites with the same pattern are combined.
        Use getPatternCounts() to obtain a list of the number of sites for each pattern and use
        getCharToPattern() to obtain a list mapping individual sites to the pattern index (i.e.,
        if site_like is the list returned by this function, and char_to_pattern is the list returned
        by getCharToPattern(), then site_like[char_to_pattern[i]] holds the site log-likelihood for
        site i.
        '''
        return self.site_likes;

    def getPatternCounts(self):
        '''
        See documentation for getSiteLikes() function.
        '''
        return self.pattern_counts;

    def getCharToPattern(self):
        '''
        See documentation for getSiteLikes() function.
        '''
        return self.char_to_pattern;

    def getSiteUF(self):
        '''
        See documentation for getSiteLikes() function.
        '''
        return site_uf;

    def __call__(self, **kwargs):
        self.set(**kwargs)
        calclike = LikeImpl(self)
        return calclike.run()

