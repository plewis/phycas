
class PyDistributionBase(object):
    def __deepcopy__(self, memo):
        c = memo.get(self)
        if c is not None:
            return c
        c = self.clone()
        memo[self] = c
        return c
        


