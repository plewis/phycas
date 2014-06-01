from phycas.Utilities.io import TreeCollection, DataSource
class DefaultData:
    "Bundle of references most loaded taxa, trees, and character matrices."
    _default_data = None
    def getInstance():
        if DefaultData._default_data is None:
            DefaultData._default_data = DefaultData()
        return DefaultData._default_data
    getInstance = staticmethod(getInstance)

    def __init__(self):
        self.__dict__["taxon_labels"] = []
        self.__dict__["characters"] = DataSource()
        self.__dict__["trees"] = TreeCollection()

    def __setattr__(self, name, value):
        raise TypeError("Assignment is not supported.")

    def set_taxon_labels(self, t):
        del self.taxon_labels[:]
        if t:
            self.taxon_labels.extend(t)
    def set_trees(self, t):
        c = self.trees
        c.reset()
        if t:
            t._copyInternals(c)
    def set_characters(self, d):
        c = self.__dict__["characters"]
        if d:
            c._setEqualTo(d)
        else:
            c._reset()

