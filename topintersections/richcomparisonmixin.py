class RichComparisonMixin(object):
    """
    Downloaded from <http://www.voidspace.org.uk/python/recipebook.shtml#comparison>.  Makes it easier to define all the rich comparison operators.
    """

    def __eq__(self, other):
        raise NotImplementedError("Equality not implemented")

    def __lt__(self, other):
        raise NotImplementedError("Less than not implemented")

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __gt__(self, other):
        return not (self.__lt__(other) or self.__eq__(other))

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)
    

    

