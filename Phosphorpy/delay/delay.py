"""
Script for the delay system
"""


class Delay:
    _calls = []

    def __init__(self):
        """
        The class delay provides a system run the code more efficient by managing the memory in a optimised way
        and by executing query simultaneously.
        """
        pass

    def __call__(self, *args, **kwargs):
        self._calls.append((args[0], args[1:], kwargs))

    def execute(self):
        for func, arg, kwarg in self._calls:
            pass
