
def lazy_attr(compute_value_func):
    """
    Use as a decorator for a method that computes a value based on unchanging data.  It turns the function into a read-only property.  It is computed only when requested, and stores the value so it only needs to be computed once.
    
    Notice that since it replaces the function with a property, you need access it without ().
    
    You need to be using new style classes for this to work!
    
    .. NOTE::
        
        The doctests don't work for this, so I put ``tsage`` instead of ``sage`` so it doesn't try to run automated tests on them.  This also kills the cool highlighting in html, too.
    
    Examples::
        
        tsage: class c(object):
            def __init__(self, x):
                self.x = x
                
            @lazy_attr
            def x_square(self):
                print "computing the square of x"
                return self.x**2
                
                
    The first time we ask for x_square, it calls the x_square method.
    
    ::
    
        tsage: t = c(5)
        tsage: t.x_square
        computing the square of x
        25
    
    In future requests, it just retrieves the stored value.
    
    ::
    
        tsage: t.x_square
        25
    """
    stored_value_name = "_" + compute_value_func.__name__ + "_stored_value"
    
    #def get_lazy_attr(self):
    #    if not self.__dict__.has_key(stored_value_name):
    #        self.__setattr__(stored_value_name, compute_value_func(self))
    #    return self.__getattribute__(stored_value_name)   
        
    def get_lazy_attr(self):
        value = self.__dict__.get(stored_value_name)
        if value == None:
            value = compute_value_func(self)
            self.__setattr__(stored_value_name, value)
        return value
    
    lazy_attr_doc = "\n\n\tThis uses the ``lazy_attr`` decorator, which turns the function into a read-only property.  It is computed only when requested, and stores the value so it only needs to be computed once."
    
    if compute_value_func.__doc__ != None:
        doc = compute_value_func.__doc__ +  lazy_attr_doc
    else:
        doc = lazy_attr_doc
        
    return property(get_lazy_attr, doc = doc)