import pandas as pd
import numpy as np

def read_pop(filename, nind, gen_size=1000, phen_size=10):
    ind = []
    raw = pd.read_csv(filename, header=None)

    for i in range(nind):
        ind_line = i*(gen_size + phen_size*gen_size)
        y = np.array(raw[ind_line:  ind_line + gen_size])
        b = np.array(raw[ind_line + gen_size:  ind_line + gen_size +phen_size*gen_size]).reshape((phen_size,gen_size))
        ind += [(y, b)]

    return ind
