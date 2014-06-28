from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions
from phycas import model

class Subset(object):
    '''An instance of this class is imported as subset in phycas.__init__

    Thus, to create a list of indices (for example, a list used to define a set
        of characters), one calls:

    subset(1,10,2)

        to get a list of [1, 3, 5, 7, 9]

    '''
    def __init__(self):
        self.start = 0
        self.stop = 0
        self.incr = 0

    def __call__(self, first_site, last_site, step_size = 1):
        """
        Creates a range using the supplied first site, last site and step size. If first_site and last_site are both zero,
        however, returns None.
        """
        if (first_site == 0) and (last_site == 0):
            return None
        else:
            self.start = int(first_site)
            self.stop = int(last_site) + 1
            self.incr = int(step_size)
            return range(self.start, self.stop, self.incr)

class Partition(PhycasCommand):
    """
    The Python class Partition handles the user interface for defining a partition.
    It has no exact counterpart on the C++ side. The class PartitionModel
    stores partitioning information on the C++ side.
    """
    def __init__(self):
        args =  (
                ('name',                    'mypart',   'Choose a name for this partitioning scheme'),
                ('fix_subset_relrates',     False,      'If True, the vector of subset relative rates will not be modified during the course of an MCMC analysis', BoolArgValidate),
                ('subset_relrates',         None,       'The vector of subset relative rates. If None, the vector of subset relative rates be set to a vector of the appropriate dimension consisting of all 1.0'),
                ('subset_relrates_prior',   None,       'The joint prior distribution for the relative rates of partition subsets. If specified, this should be a RelativeRateDistribution with dimension equal to the number of subsets. If None, a flat RelativeRateDistribution prior will be generated.'),
                )

        # Specify output options
        PhycasCommand.__init__(self, args, "partition", "Sets up a data partitioning scheme, allowing different models to be applied to different data subsets.")

        # set to True if data for at least one character is provided
        self.__dict__['is_data'] = False

        # a list containing (name, sitelist, model) tuples
        self.__dict__['subset'] = []

        # a list of integers representing the index of the model for each site
        self.__dict__['sitemodel'] = []

        # needed for using phycassert
        self.__dict__['cf'] = CommonFunctions(self)

    def hidden():
        """
        Overrides the PhycasCommand.hidden method to keep Partitions's name from being displayed
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return True

    hidden = staticmethod(hidden)

    def flatten(self, alist):
        #print 'flattend: alist =', alist

        flattened = []
        for s in alist:
            if s.__class__.__name__ == 'list' or s.__class__.__name__ == 'tuple':
                # an element of the supplied list is itself a list or tuple, so recurse
                s = self.flatten(s)
                flattened.extend(s)
            else:
                if s.__class__.__name__ == 'int':
                    self.cf.phycassert(s > 0, 'elements of site lists must be greater than 0; you supplied the value %d' % s)
                    flattened.append(s)
                else:
                    self.cf.phycassert(False, 'elements of site lists must be of type int; you supplied a value of type %s' % s.__class__)
        return flattened

    def getSiteModelVector(self):
        return self.sitemodel

    def getModels(self):
        """
        Creates and returns a list of just the models stored in subset.
        """
        return [m for n,s,m in self.subset]

    def noData(self):
        return not self.is_data

    def resetPartition(self):
        """
        Returns object to the state it was in when first constructed: subset and sitemodel are both empty lists
        and is_data is False.
        """
        self.is_data = False
        self.subset = []
        self.sitemodel = []

    def validate(self, nsites):
        """
        If user has not specified a partition, create a default partition now using the current model.
        """
        if nsites > 0:
            self.is_data = True
        else:
            self.is_data = False

        # Better check to make sure user did not define empty subsets
        if self.is_data:
            num_empty_subsets = 0
            for i,(n,s,m) in enumerate(self.subset):
                if s is None or len(s) == 0:
                    num_empty_subsets += 1
            self.cf.phycassert(num_empty_subsets == 0 or num_empty_subsets == len(self.subset), 'There should either be no empty subsets, or all subsets should be empty. In this case, %d of %d subsets were empty' % (num_empty_subsets, len(self.subset)))

        if self.is_data:
            if len(self.subset) < 1:
                self.sitemodel = []
                self.addSubset(range(1, nsites + 1), model, 'default')
        else:
            if self.subset is None or len(self.subset) == 0:
                self.addSubset(range(1,1), model, 'default')

    def addSubset(self, sites, model_for_sites, name = None):
        """
        Applies a model (model_for_sites) to a list of sites (sites). For example,

        model.type = 'hky'
        hky = model()
        partition.addSubset(subset(1,100,3), hky, 'first')

        """
        # the index of the new subset is simply the length of self.subset
        subset_index = len(self.subset)
        if name is None:
            name = 'subset%d' % (1 + subset_index)

        self.cf.output('Processing subset %s...' % name)

        # make a copy of the supplied list in order to flatten it (if necessary)
        if sites is None:
            # No sites were specified. User apparently wishes to just explore the prior for this subset
            sitelist = []
        else:
            sitelist = self.flatten(sites)
            sitelist.sort()

        # expand sitemodel list if last site in sorted `sites' list is larger
        # than the last site in self.sitemodel
        curr_size = len(self.sitemodel)
        needed_size = len(sitelist) > 0 and sitelist[-1] or 0
        if curr_size < needed_size:
            xtra = [-1]*(needed_size - curr_size)
            self.sitemodel.extend(xtra)

        # add sites in sitelist to sitemodel (the master list)
        for s in sitelist:
            assigned_index = self.sitemodel[s - 1]
            if assigned_index == subset_index:
                self.cf.phycassert(False, 'site %d has already been assigned a model by the current subset (%s)' % (s,name))
            elif assigned_index >= 0:
                assigned_subset = self.subset[assigned_index]
                self.cf.phycassert(False, 'site %d has already been assigned a model by subset %s' % (s,assigned_subset[0]))
            self.sitemodel[s - 1] = subset_index

        self.subset.append((name, sitelist, model_for_sites))

    def save(self):
        return copy.deepcopy(self)

    def __deepcopy__(self, memo):
        # memo is a dictionary used to avoid infinite recursion
        # memo keeps track of objects already copied
        c = memo.get(self)
        if c:
            return c

        # self not already copied, so copy it now
        new_partition = Partition()
        new_partition.name = copy.deepcopy(self.name, memo)
        memo[self] = new_partition
        return new_partition

    def handleAllMissingSites(self, all_missing):
        """
        Go through the specified list of sites with all missing data and revise
        numbers of sites in affected partition subsets accordingly.
        """
        for i,(name,sitelist,model) in enumerate(self.subset):
            print name,len(sitelist)
            for site in all_missing:
                if site in sitelist:
                    sitelist.remove(site)

    def partitionReport(self, reporter):
        """
        Outputs a summary of the partitioning scheme.
        """
        #cf = CommonFunctions(self)
        reporter.output('Partition report:')
        reporter.output('%12s %12s %12s   %s' % ('subset', 'size', 'model', 'name'))
        for i,(n,s,m) in enumerate(self.subset):
            model_descr = m.type.upper()
            if m.pinvar_model:
                model_descr += '+I'
            if m.num_rates > 1:
                model_descr += '+G'
            reporter.output('%12d %12d %12s   %s' % (i+1, len(s), model_descr, n))

    def getSubsetProportions(self):
        """
        Obtains a list of the number of sites in each subset from getSubsetSizes, then
        returns a normalized vector in which these subset sizes are divided by the
        total number of sites. If all subsets have size 0, then the proportion of each
        subset is set to 1/number_of_subsets.
        """
        proportions = []

        if self.getNumSubsets() > 0:
            size_vect = self.getSubsetSizes()

            n = sum(size_vect)
            if n == 0:
                num_subsets = len(size_vect)
                p = 1.0/float(num_subsets)
                proportions = [p]*num_subsets
            else:
                for s in size_vect:
                    proportions.append(float(s)/float(n))

        return proportions;

    def getSubsetSizes(self):
        """
        Extracts the length of the sitelist for each subset in self.subset, which
        is a list of tuples (name, sitelist, model).
        """
        size_vect = []
        for s in self.subset:
            size_vect.append(len(s[1]))
        return size_vect;

    def getNormalizedSubsetRelRates(self):
        sizes = self.getSubsetSizes()
        n = sum(sizes)
        x = sum([rr*float(sz)/float(n) for rr,sz in zip(self.subset_relrates,sizes)])
        normrr = [rr/x for rr in self.subset_relrates]
        return normrr

    def getNumSubsets(self):
        """
        Returns the length of the subset vector.
        """
        return len(self.subset);

    def __call__(self, **kwargs):
        """
        This function should perform some sanity checks, but right now it does nothing.
        """
        cf = CommonFunctions(self)
        cf.phycassert(not model.update_freqs_separately, 'update_freqs_separately is no longer allowed')
        cf.phycassert(not model.update_relrates_separately, 'update_freqs_separately is no longer allowed')
        #self.set(**kwargs)
        #return self.partitionReport()
