class NestedAttributes(object):
    """ Container with nested attribute interface

        Each item in the container is identified by a name
        "something.this.that". The attributes can be accessed via the "."
        operator, as within a normal python object hierarchy.

        This container is not (easily) writable.

        The input is a dictionay, as follows:

        >>> nested = NestedAttributes({                 \
                ('this', 'that'): 'first object',       \
                ('this', 'there'): 'second object',     \
                ('anther', 'something'): 'third object' \
            })
        >>> nested.this.that
        'first object'
        >>> nested.this.there
        'second object'
        >>> nested.anther.something
        'third object'

        Hierarchies where an item is both a branch and a leaf are illegal:

        >>> nested = NestedAttributes({                    \
                ('this', 'that'): 'first object',          \
                ('this', 'that', 'there'): 'other object', \
            })
        ValueError: Ambiguous object hierarchy

        It is not clear whether `nested.this.that` should return the string
        'first object', or an object with an attribute 'there'.

        The keys must be tuples.
    """
    def __init__(self, chain):
        super(NestedAttributes, self).__init__()

        # Sanity check
        if any(not isinstance(k, tuple) for k in chain):
            raise TypeError("Keys of input dictionary should be tuples")
        if not set(chain).isdisjoint(k[:-1] for k in chain):
            raise ValueError("Ambiguous object hierarchy")

        self._chain = chain
        """ Operator chain indicating location and name of c function """

    def __dir__(self):
        return list(set([u[0] for u in self._chain.keys()]))

    def __getattr__(self, name):
        chain = {k[1:]: v for k, v in self._chain.items() if k[0] == name}
        if len(chain) == 0:
            raise AttributeError(name)
        return NestedAttributes(chain) if len(chain) > 1 \
            else next(iter(chain.values()))
