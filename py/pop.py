#! /usr/bin/env python
"""Usage: ./pop.py [options]

Options:
   -h --help              show this
   -l loci                number of loci [default: 200]
   -p traits              number of traits [default: 10]
   -m mu                  genetic mutation rate [default: 5e-4]
   -b mu_b                ontogenetic mutation rate [default: 1e-4]
   -n ne                  population size [default: 25000]
   -s sigma               mutation size [default: 0.2]
   -e amb                 enviromental noise [default: 0.8]
   -v omega_var           selection variance [default: 1.0]
   -o omega_mat           selection correlation matrix [default: omega.csv]
   -t time                number of generations [default: 1]
   -d delta_S             change in optimal per generation [default: 0.0]
"""

import numpy as np
import os
import errno
import multiprocessing as mp
import functools as f
#import matplotlib.pylab as mpl
from docopt import docopt


class Individual:
    """class for Individuals in a Population"""
    def __init__(self, m, p, amb, mu, mu_b, sigma):
        self.m = m                # number of loci
        self.p = p                # number of traits
        self.amb = amb            # environmental noise
        self.mu = mu              # genetic effects mutation rate
        self.mu_b = mu_b          # ontogenetic effects mutation rate
        self.sigma = sigma        # genetic mutation variance

    def generate(self):
        """creates an ind with Individual parameters"""
        b = np.concatenate((np.ones(self.m * self.p),   # binary ontogenetic
                           np.zeros(self.m * self.p)))  # matrix
        np.random.shuffle(b)
        b = b.reshape((self.p, 2 * self.m))
        y = np.random.normal(0, self.sigma, size=2 * self.m)  # gene vector
        x = np.dot(b, y)                       # additive effects vector
        z = x + np.random.normal(0, self.amb,
                                 self.p)       # phenotipic values vector
        return {'y': y, 'x': x, 'z': z, 'b': b}

    def mutate(self, ind):
        """Mutates an ind with Individual parameters"""
        mutation_number_y = np.random.binomial(2 * self.m, self.mu)
        if (mutation_number_y > 0):
            mutation_vector = np.concatenate((
                np.random.normal(self.sigma, size=mutation_number_y),
                np.zeros(2 * self.m - mutation_number_y)))
            np.random.shuffle(mutation_vector)
            ind['y'] = ind['y'] + mutation_vector
        mutation_number_b = np.random.binomial(2 * self.m * self.p,
                                               self.mu_b)
        if (mutation_number_b > 0):
            mutation_mask = np.concatenate((
                np.ones(mutation_number_b),
                np.zeros(2 * self.m * self.p - mutation_number_b)))
            np.random.shuffle(mutation_mask)
            mutation_mask = mutation_mask.reshape((self.p, 2 * self.m))
            ind['b'] = (ind['b'] + mutation_mask) % 2
        return ind

    def fitness(self, ind, omega, teta):
        """calculates ind fitness from population"""
        delta_s = ind['z'] - teta
        return np.exp(-np.dot(delta_s, np.linalg.solve(omega, delta_s)))

    def cross(self, ind_1, ind_2):
        """Generates a new individual by crossing ind_1 and ind_2"""
        b = np.zeros((self.p, 2 * self.m),
                     dtype=float)              # binary ontogenetic matrix
        y = np.zeros(2 * self.m, dtype=float)  # gene vector
        alele = np.random.randint(0, 2, size=2 * self.m)
        for locus in xrange(self.m):
            y[2 * locus] = ind_1['y'][2 * locus + alele[2 * locus]]
            y[2 * locus + 1] = ind_2['y'][2 * locus + alele[2 * locus + 1]]
            b[:, 2 * locus] = ind_1['b'][:, 2 * locus + alele[2 * locus]]
            b[:, 2 * locus + 1] = ind_2['b'][:, 2 * locus + alele[2 * locus + 1]]
        x = np.dot(b, y)                       # additive effects vector
        z = x + np.random.normal(0, self.amb,
                                 self.p)       # phenotipic values vector
        return {'y': y, 'x': x, 'z': z, 'b': b}


def mut(p, i):
    return p.indmod.mutate(p.pop[i])

def fit(p, i):
    fitness = p.indmod.fitness(p.pop[i],
                               p.omega,
                               p.teta)
    if (np.isnan(fitness) or np.isinf(fitness)):
        fitness = 0.0
    return fitness

def cross(p, c):
    return p.indmod.cross(p.pop[c[0]], p.pop[c[1]])

class Population:
    """class for population of Individuals"""
    def __init__(self, n_e, teta, delta_s, omega, indmod):
        self.indmod = indmod
        self.n_e = n_e
        self.teta = teta
        self.omega = omega
        self.current_gen = 0
        self.modules = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
        self.pop = [indmod.generate() for k in xrange(n_e)]
        self.fitness = np.ones(n_e) / n_e
        self.pop_name = ('./dats/Ne.' + str(n_e) + '-mp.' + str(indmod.m) + '_'
                         + str(indmod.p) + '-mu_muB.' + str(indmod.mu) + '_'
                         + str(indmod.mu_b) + '-Delta_S.' + str(delta_s))
        self.out_files = self.set_outfile()

    def mutate(self):
        """mutates every individual of population"""
        #TODO Obviously parallel

        self.pop = pool.map(f.partial(mut, self), range(0, self.n_e))

    def pupdate_fitness(self):
        self.fitness = np.array(pool.map(f.partial(fit, self), range(0,
            self.n_e)))
        fitness_total = self.fitness.sum()
        if (fitness_total == 0.0):
            self.fitness = np.ones(self.n_e) / self.n_e
        else:
            self.fitness /= fitness_total
 
    def update_fitness(self):
        """calculates the fitness of every individual of population"""
        #TODO Obviously parallel
        for k in xrange(self.n_e):
            self.fitness[k] = self.indmod.fitness(self.pop[k],
                                                  self.omega,
                                                  self.teta)
            if (np.isnan(self.fitness[k]) or np.isinf(self.fitness[k])):
                self.fitness[k] = 0.0
        fitness_total = self.fitness.sum()
        if (fitness_total == 0.0):
            self.fitness = np.ones(self.n_e) / self.n_e
        else:
            self.fitness /= fitness_total

    def next_generation(self, delta_s):
        """creates next generation by mutating then crossing with probability
        proportional do fitness"""
        self.mutate()
        self.update_fitness()
        sires = np.random.choice(self.n_e, size=self.n_e,
                                 p=self.fitness, replace=True)
        dames = np.random.choice(self.n_e, size=self.n_e,
                                 p=self.fitness, replace=True)

	pairs = zip(sires, dames)
        new_pop = []
        #TODO Obviously parallel
        #for k in xrange(self.n_e):
        #    new_pop.append(self.indmod.cross(self.pop[sires[k]],
        #                                     self.pop[dames[k]]))

	new_pop = pool.map(f.partial(cross, self), pairs)
        self.pop = new_pop
        self.current_gen += 1
        self.teta += delta_s

    def moments(self):
        """Calculate covariance and correlation matrices,
        trait, genotipic and ontogenetic means"""
        zs = np.array([ind['z'] for ind in self.pop])
        xs = np.array([ind['x'] for ind in self.pop])
        ys = np.array([ind['y'] for ind in self.pop])
        bs = np.array([ind['b'] for ind in self.pop])
        ymean = ys.mean(axis=0)
        zmean = zs.mean(axis=0)
        xmean = xs.mean(axis=0)
        ymean = ys.mean(axis=0)
        bmean = bs.mean(axis=0)
        phenotipic = np.cov(zs, rowvar=0, bias=1)
        genetic = np.cov(xs, rowvar=0, bias=1)
        heridability = (genetic[np.diag_indices_from(genetic)] /
                        phenotipic[np.diag_indices_from(phenotipic)])
        corr_phenotipic = np.corrcoef(zs, rowvar=0, bias=1)
        corr_genetic = np.corrcoef(xs, rowvar=0, bias=1)
        avgP = avg_ratio(corr_phenotipic, self.modules)
        avgG = avg_ratio(corr_genetic, self.modules)
        return {'y.mean': ymean,
                'b.mean': bmean,
                'z.mean': zmean,
                'x.mean': xmean,
                'P': phenotipic,
                'G': genetic,
                'h2': heridability,
                'avgP': avgP,
                'avgG': avgG,
                'corrP': corr_phenotipic,
                'corrG': corr_genetic}

    def print_moments(self):
        mats = self.moments()
        gen = self.current_gen
        mat_print(self.teta, 'vector', self.out_files['teta'], gen)
        mat_print(mats['corrG'], 'tri', self.out_files['corrG'], gen)
        mat_print(mats['corrP'], 'tri', self.out_files['corrP'], gen)
        mat_print(mats['G'], 'diag', self.out_files['varG'], gen)
        mat_print(mats['P'], 'diag', self.out_files['varP'], gen)
        mat_print(mats['h2'], 'vector', self.out_files['varH'], gen)
        mat_print(mats['avgP'], 'vector', self.out_files['avgP'], gen)
        mat_print(mats['avgG'], 'vector', self.out_files['avgG'], gen)
        mat_print(mats['h2'], 'vector', self.out_files['varH'], gen)
        mat_print(mats['y.mean'], 'vector', self.out_files['y.mean'], gen)
        mat_print(mats['z.mean'], 'vector', self.out_files['z.mean'], gen)
        mat_print(mats['x.mean'], 'vector', self.out_files['x.mean'], gen)
        mat_print(mats['b.mean'], 'total', self.out_files['b.mean'], gen)

    def set_outfile(self):
        try:
            os.makedirs('./dats/')
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        varG = open(self.pop_name + '-varG.dat', "w")
        varP = open(self.pop_name + '-varP.dat', "w")
        varH = open(self.pop_name + '-varH.dat', "w")
        corrP = open(self.pop_name + '-corrP.dat', "w")
        corrG = open(self.pop_name + '-corrG.dat', "w")
        avgP = open(self.pop_name + '-avgP.dat', "w")
        avgG = open(self.pop_name + '-avgG.dat', "w")
        z_traits = open(self.pop_name + '-traits_z.dat', "w")
        y_traits = open(self.pop_name + '-traits_y.dat', "w")
        x_traits = open(self.pop_name + '-traits_x.dat', "w")
        b_matrix = open(self.pop_name + '-b.dat', "w")
        teta = open(self.pop_name + '-teta.dat', "w")
        return {'varG': varG,
                'varP': varP,
                'varH': varH,
                'corrP': corrP,
                'corrG': corrG,
                'avgP': avgP,
                'avgG': avgG,
                'teta': teta,
                'b.mean': b_matrix,
                'z.mean': z_traits,
                'y.mean': y_traits,
                'x.mean': x_traits}


def mat_print(matrix, out_format, output, generation):
    s = str(generation) + ' '
    if (out_format == 'vector'):
        n = len(matrix)
        for i in xrange(n):
                s += str(matrix[i]) + ' '
    else:
        n, m = matrix.shape

    if (out_format == 'total'):
        for i in xrange(n):
            for j in xrange(m):
                s += str(matrix[i, j]) + ' '
    if (out_format == 'tri'):
        for i in xrange(1, n):
            for j in xrange(i):
                s += str(matrix[i, j]) + ' '
    if (out_format == 'diag'):
        for i in xrange(n):
                s += str(matrix[i, i]) + ' '
    s += '\n'
    output.write(s)
    output.flush()


def avg_ratio(matrix, modules):
    """Calculates average correlation given between and within modules
    given a correlation matrix and a list of modules"""
    n, p = matrix.shape
    mask = np.zeros((n, n), int)
    avgs = []
    k = 0  # number of modules
    for traits in modules:
        k = k + 1
        for i in traits:
            for j in traits:
                if i == j:
                    mask[i, j] = -17  # any negative integer
                else:
                    mask[i, j] = k
    for module in xrange(k + 1):
        avgs.append(matrix[mask == module].mean())
    # avg ratio, avg+, avg-, avg by module
    return [np.mean(avgs[1:]) / avgs[0]] + [np.mean(avgs[1:])] + avgs


pool = mp.Pool(4)
pool.map(id, range(mp.cpu_count()))

def main(options):
    teta_init = np.zeros(int(options['-p']))
    delta_teta = np.ones(int(options['-p'])) * float(options['-d'])
    omega = np.genfromtxt(options['-o'])
    i = Individual(int(options['-l']),
                   int(options['-p']),
                   float(options['-e']),
                   float(options['-m']),
                   float(options['-b']),
                   float(options['-s']))
    p = Population(int(options['-n']),
                   teta_init,
                   float(options['-d']),
                   omega,
                   i)
    for generation in xrange(int(options['-t'])):
        print generation
        p.next_generation(delta_teta)
        p.print_moments()

if __name__ == '__main__':
    options = docopt(__doc__)
    main(options)
