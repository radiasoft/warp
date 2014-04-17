class Noparens:
    """Creates a class representation of functions that allow
    the function to be called without the added '()' when there are no
    arguments. Note that in that case, an extra blank line is output.

    """
    def __init__(self, func):
        self.func = func
        self.__doc__ = func.__doc__

    def __repr__(self):
        self.func()
        return ''
    __str__ = __repr__

    def __call__(self, *k, **kw):
        return apply(self.func, k, kw)
