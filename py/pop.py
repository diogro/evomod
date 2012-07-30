import numpy as np


class Individual:
    """class for Individuals in a Population"""
    def __init__(self, n, p):
        self.n = n # number of
        self.p = p
        self.B = np.zeros(n, p)
        
