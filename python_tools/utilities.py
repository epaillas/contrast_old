class Utilities:

    def __init__(self):
        return

    @staticmethod
    def next_pow_two(n):
        '''
        Returns the largest power of two
        smaller than a given positive integer.
        '''
        i = 1
        while i < n:
            i = i << 1
        i = i >> 1
        return i