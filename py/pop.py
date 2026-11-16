#!/usr/bin/env python3
"""Individual-based simulation of the evolution of modularity.

Python port of the C model in ../c-gsl (Melo & Marroig 2015, PNAS).

Usage: pop.py [options]

Options:
   -h --help              show this
   -l loci                number of loci [default: 200]
   -p traits              number of traits [default: 10]
   -m mu                  genetic mutation rate [default: 5e-4]
   -b mu_b                ontogenetic mutation rate [default: 1e-4]
   -n ne                  population size [default: 2500]
   -s sigma               mutation size [default: 0.2]
   -e amb                 enviromental noise (standard deviation) [default: 0.8]
   -v omega_var           selection variance, used when no -o matrix is given [default: 1.0]
   -o omega_mat           selection covariance matrix file [default: none]
   -t time                number of generations [default: 1]
   -d delta_S             change in optimal per generation [default: 0.0]
   -r seed                random seed [default: none]
"""

import os

import numpy as np
from docopt import docopt


class Individual:
    """Parameters and operations shared by all Individuals in a Population.

    An individual is a dict with
      y: allelic effects, length 2m (alleles 2k and 2k+1 belong to locus k)
      b: binary ontogenetic matrix, p x 2m
      x: additive values, b . y
      z: phenotype, x + environmental noise
    """

    def __init__(self, m, p, amb, mu, mu_b, sigma, rng):
        self.m = m                # number of loci
        self.p = p                # number of traits
        self.amb = amb            # environmental noise
        self.mu = mu              # genetic effects mutation rate
        self.mu_b = mu_b          # ontogenetic effects mutation rate
        self.sigma = sigma        # genetic mutation standard deviation
        self.rng = rng
        self.elocus = np.arange(0, 2 * self.m, 2)
        self.olocus = np.arange(1, 2 * self.m, 2)

    def phenotype(self, b, y):
        x = b @ y                                            # additive effects
        z = x + self.rng.normal(0, self.amb, self.p)         # phenotype
        return {'y': y, 'x': x, 'z': z, 'b': b}

    def generate(self):
        """creates an ind with Individual parameters"""
        b = np.concatenate((np.ones(self.m * self.p),        # binary ontogenetic
                            np.zeros(self.m * self.p)))      # matrix
        self.rng.shuffle(b)
        b = b.reshape((self.p, 2 * self.m))
        y = self.rng.normal(0, self.sigma, size=2 * self.m)  # gene vector
        return self.phenotype(b, y)

    def mutate(self, ind):
        """Mutates an ind in place (x and z are recomputed at crossing)"""
        n_y = self.rng.binomial(2 * self.m, self.mu)
        if n_y > 0:
            pos = self.rng.choice(2 * self.m, n_y, replace=False)
            ind['y'] = ind['y'].copy()
            ind['y'][pos] += self.rng.normal(0, self.sigma, size=n_y)
        n_b = self.rng.binomial(2 * self.m * self.p, self.mu_b)
        if n_b > 0:
            pos = self.rng.choice(2 * self.m * self.p, n_b, replace=False)
            b = ind['b'].copy().ravel()
            b[pos] = 1 - b[pos]
            ind['b'] = b.reshape((self.p, 2 * self.m))
        return ind

    def cross(self, ind_1, ind_2):
        """Generates a new individual by crossing ind_1 and ind_2"""
        b = np.zeros((self.p, 2 * self.m))                   # binary ontogenetic matrix
        y = np.zeros(2 * self.m)                             # gene vector
        alele = self.rng.integers(0, 2, size=2 * self.m)

        y[self.elocus] = ind_1['y'][self.elocus + alele[self.elocus]]
        y[self.olocus] = ind_2['y'][self.elocus + alele[self.olocus]]
        b[:, self.elocus] = ind_1['b'][:, self.elocus + alele[self.elocus]]
        b[:, self.olocus] = ind_2['b'][:, self.elocus + alele[self.olocus]]
        return self.phenotype(b, y)


class Population:
    """class for population of Individuals"""

    def __init__(self, n_e, teta, delta_s, omega, indmod, out_dir='./dats/'):
        self.indmod = indmod
        self.rng = indmod.rng
        self.n_e = n_e
        self.teta = np.array(teta, dtype=float)
        self.omega = omega
        self.current_gen = 0
        half = indmod.p // 2
        self.modules = [list(range(half)), list(range(half, indmod.p))]
        self.pop = [indmod.generate() for _ in range(n_e)]
        self.fitness = np.ones(n_e) / n_e
        self.pop_name = os.path.join(
            out_dir,
            'Ne.' + str(n_e) + '-mp.' + str(indmod.m) + '_' + str(indmod.p)
            + '-mu_muB.' + str(indmod.mu) + '_' + str(indmod.mu_b)
            + '-Delta_S.' + str(delta_s))
        os.makedirs(out_dir, exist_ok=True)
        self.out_files = self.set_outfile()

    def mutate(self):
        """mutates every individual of population"""
        self.pop = [self.indmod.mutate(ind) for ind in self.pop]

    def update_fitness(self):
        """Gaussian fitness, w = exp(-1/2 (z - theta)' Omega^-1 (z - theta)),
        normalised to sum to one"""
        dz = np.array([ind['z'] for ind in self.pop]) - self.teta
        fitness = np.exp(-0.5 * np.einsum('ij,ji->i', dz, np.linalg.solve(self.omega, dz.T)))
        fitness[~np.isfinite(fitness)] = 0.0
        total = fitness.sum()
        if total == 0.0:
            self.fitness = np.ones(self.n_e) / self.n_e
        else:
            self.fitness = fitness / total

    def next_generation(self, delta_s):
        """creates next generation by mutating then crossing with probability
        proportional do fitness"""
        self.update_fitness()
        self.mutate()
        sires = self.rng.choice(self.n_e, size=self.n_e, p=self.fitness, replace=True)
        dames = self.rng.choice(self.n_e, size=self.n_e, p=self.fitness, replace=True)
        self.pop = [self.indmod.cross(self.pop[s], self.pop[d])
                    for s, d in zip(sires, dames)]
        self.current_gen += 1
        self.teta += delta_s

    def moments(self):
        """Calculate covariance and correlation matrices,
        trait, genotipic and ontogenetic means"""
        zs = np.array([ind['z'] for ind in self.pop])
        xs = np.array([ind['x'] for ind in self.pop])
        ys = np.array([ind['y'] for ind in self.pop])
        bs = np.array([ind['b'] for ind in self.pop])
        phenotipic = np.cov(zs, rowvar=False, bias=True)
        genetic = np.cov(xs, rowvar=False, bias=True)
        corr_phenotipic = np.corrcoef(zs, rowvar=False)
        corr_genetic = np.corrcoef(xs, rowvar=False)
        return {'y.mean': ys.mean(axis=0),
                'b.mean': bs.mean(axis=0),
                'z.mean': zs.mean(axis=0),
                'x.mean': xs.mean(axis=0),
                'P': phenotipic,
                'G': genetic,
                'h2': np.diag(genetic) / np.diag(phenotipic),
                'avgP': avg_ratio(corr_phenotipic, self.modules),
                'avgG': avg_ratio(corr_genetic, self.modules),
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
        mat_print(mats['y.mean'], 'vector', self.out_files['y.mean'], gen)
        mat_print(mats['z.mean'], 'vector', self.out_files['z.mean'], gen)
        mat_print(mats['x.mean'], 'vector', self.out_files['x.mean'], gen)
        mat_print(mats['b.mean'], 'total', self.out_files['b.mean'], gen)

    def set_outfile(self):
        suffixes = {'varG': '-varG.dat',
                    'varP': '-varP.dat',
                    'varH': '-varH.dat',
                    'corrP': '-corrP.dat',
                    'corrG': '-corrG.dat',
                    'avgP': '-avgP.dat',
                    'avgG': '-avgG.dat',
                    'teta': '-teta.dat',
                    'b.mean': '-b.dat',
                    'z.mean': '-traits_z.dat',
                    'y.mean': '-traits_y.dat',
                    'x.mean': '-traits_x.dat'}
        return {k: open(self.pop_name + s, 'w') for k, s in suffixes.items()}

    def close(self):
        for f in self.out_files.values():
            f.close()


def mat_print(matrix, out_format, output, generation):
    matrix = np.asarray(matrix)
    if out_format == 'vector':
        values = matrix
    elif out_format == 'total':
        values = matrix.ravel()
    elif out_format == 'tri':
        values = matrix[np.tril_indices_from(matrix, k=-1)]
    elif out_format == 'diag':
        values = np.diag(matrix)
    else:
        raise ValueError(out_format)
    output.write(str(generation) + ' ' + ''.join(str(v) + ' ' for v in values) + '\n')
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
    for module in range(k + 1):
        avgs.append(matrix[mask == module].mean())
    # avg ratio, avg+, avg-, avg by module
    return [np.mean(avgs[1:]) / avgs[0]] + [np.mean(avgs[1:])] + avgs


def main(options):
    p = int(options['-p'])
    seed = None if options['-r'] == 'none' else int(options['-r'])
    rng = np.random.default_rng(seed)
    teta_init = np.zeros(p)
    delta_teta = np.ones(p) * float(options['-d'])
    if options['-o'] == 'none':
        omega = float(options['-v']) * np.eye(p)
    else:
        omega = np.genfromtxt(options['-o'])
    i = Individual(int(options['-l']),
                   p,
                   float(options['-e']),
                   float(options['-m']),
                   float(options['-b']),
                   float(options['-s']),
                   rng)
    pop = Population(int(options['-n']),
                     teta_init,
                     float(options['-d']),
                     omega,
                     i)
    for generation in range(int(options['-t'])):
        print(generation)
        pop.next_generation(delta_teta)
        pop.print_moments()
    pop.close()


if __name__ == '__main__':
    main(docopt(__doc__))
