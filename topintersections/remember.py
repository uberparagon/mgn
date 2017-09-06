"""
Some automatically remembering decorators.
"""

def remember_convert_args(converter, storing_dict = dict()):
    """
    The idea is that if a function has already been called with the same arguments, you would want to just retrieve the computed value instead of recomputing it.
    
    This decorator allows you to specify a converter that will convert the arguments passed into your function into an appropriate format for use as a dictionary key.  You can also specify a dictionary yourself, so, for example, two functions could share a dictionary.
    
    This also has save and load features, but I think the current implementation of the program doesn't need them.
    """
    class remember_class:
        def __init__(self, f):
            self.func = f
            self.stored = storing_dict
            
        def __call__(self, *args):
            key = converter(*args)
            value = self.stored.get(key, None)
            if value != None:
                #print("retrieving stored value for input {0}, returning {1}".format(key,value))
                return value
            #print("computing value")
            value = self.func(*args)
            self.stored[key] = value
            return value
            
        def save(self, filename):
            with open(filname, "w") as f:
                pickle.dump(self.stored,f)
                
        def load(self, filename):
            with open(filename, "r") as f:
                self.stored = pick.load(f)
            
    def use_converter(func):
        return remember_class(func)
            
    return use_converter
    
