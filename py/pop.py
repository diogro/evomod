#! /usr/bin/env python
"""Usage: ./pop.py [options]

Options:
   -h --help                show this
   -m loci                  number of loci [default: 20]
   -p traits                number of traits [default: 10]
   -mu mu                   genetic mutation rate [default: 5e-4]
   -mu_b mu_b               ontogenetic mutation rate [default: 1e-4]
   -ne ne                   population size [default: 2500]
   -s sigma                 mutation size [default: 0.2]
   -e amb                   enviromental noise [default: 0.8]
   -o omega                 selection variance [default: 1.0]
   -omega_mat file_name     selection correlation matrix [default: omega.csv]
   -t time                  number of generations [default: 1]
   -ds delta_S              change in optimal per generation [default: 0.0]
"""

import numpy as np
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
        b = np.ones((self.p, 2 * self.m),
                    dtype=float)              # binary ontogenetic matrix
        y = np.zeros(2 * self.m, dtype=float)  # gene vector
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

    def fitness(self, ind, omega, teta):
        """calculates ind fitness from population"""
        delta_s = ind['z'] - teta
        return np.exp(-np.dot(delta_s, np.linalg.solve(omega, delta_s)))

    def cross(self, ind_1, ind_2):
        """Generates a new individual by crossing ind_1 and ind_2"""
        b = np.zeros((self.p, 2 * self.m),
                     dtype=float)              # binary ontogenetic matrix
        y = np.zeros(2 * self.m, dtype=float)  # gene vector
        for locus in range(self.m):
            alele_1 = np.random.randint(0, 2)
            alele_2 = np.random.randint(0, 2)
            y[2 * locus] = ind_1['y'][2 * locus + alele_1]
            y[2 * locus + 1] = ind_2['y'][2 * locus + alele_2]
            b[:, 2 * locus] = ind_1['b'][:, 2 * locus + alele_1]
            b[:, 2 * locus + 1] = ind_2['b'][:, 2 * locus + alele_2]
        x = np.dot(b, y)                       # additive effects vector
        z = x + np.random.normal(0, self.amb,
                                 self.p)       # phenotipic values vector
        return {'y': y, 'x': x, 'z': z, 'b': b}


class Population:
    """class for population of Individuals"""
    def __init__(self, n_e, teta, omega, indmod):
        self.indmod = indmod
        self.n_e = n_e
        self.teta = teta
        self.omega = omega
        self.current_gen = 0
        self.pop = [indmod.generate() for k in range(n_e)]
        self.fitness = np.ones(n_e) / n_e

    def mutate(self):
        """mutates every individual of population"""
        #TODO Obviously parallel
        for k in range(self.n_e):
            self.indmod.mutate(self.pop[k])

    def update_fitness(self):
        """calculates the fitness of every individual of population"""
        #TODO Obviously parallel
        for k in range(self.n_e):
            self.fitness[k] = self.indmod.fitness(self.pop[k],
                                                  self.omega,
                                                  self.teta)
            if (np.isnan(self.fitness[k]) or np.isinf(self.fitness[k])):
                self.fitness[k] = 0.0
        fitness_total = self.fitness.sum()
        if (fitness_total == 0.0):
            self.fitness = np.ones(self.indmod.p) / self.indmod.p
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
        new_pop = []
        #TODO Obviously parallel
        for k in range(self.n_e):
            new_pop.append(self.indmod.cross(self.pop[sires[k]],
                                             self.pop[dames[k]]))
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
        phenotipic = np.cov(zs.transpose())
        genetic = np.cov(xs.transpose())
        outer_diagonal = phenotipic[np.diag_indices_from(phenotipic)]
        outer_diagonal = np.sqrt(outer_diagonal[:, np.newaxis] *
                                 outer_diagonal)
        corr_phenotipic = phenotipic / outer_diagonal
        outer_diagonal = genetic[np.diag_indices_from(genetic)]
        outer_diagonal = np.sqrt(outer_diagonal[:, np.newaxis] *
                                 outer_diagonal)
        corr_genetic = genetic / outer_diagonal
        return {'y.mean': ymean,
                'b.mean': bmean,
                'z.mean': zmean,
                'x.mean': xmean,
                'P': phenotipic,
                'G': genetic,
                'corrP': corr_phenotipic,
                'corrG': corr_genetic}


def main(options):
    teta = np.ones(int(options['-p']))
    omega = np.genfromtxt(options['-omega_mat'])
    i = Individual(int(options['-m']),
                   int(options['-p']),
                   float(options['-e']),
                   float(options['-mu']),
                   float(options['-mu_b']),
                   float(options['-s']))
    p = Population(int(options['-ne']),
                   teta,
                   omega,
                   i)
    for generation in range(int(options['-t'])):
        print generation
        p.next_generation(teta * float(options['-ds']))

if __name__ == '__main__':
    options = docopt(__doc__)
    main(options)
