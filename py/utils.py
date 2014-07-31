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

def calc_fitness(inds, omega, theta, amb):

    fit_inds = np.zeros(len(inds))

    for i, ind in enumerate(inds):
        z = np.dot(ind[1], ind[0]) + np.random.normal(0, amb, ind[1].shape[0]).reshape((10,1))
        delta_s = z - theta.reshape((10,1))
        fit_inds[i] = np.exp(-np.dot(delta_s.T, np.linalg.solve(omega, delta_s)))

    return fit_inds

