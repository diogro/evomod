import numpy as np
import random as rand



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
        """creates an ind dict with Individual parameters"""
        b = np.zeros(2 * self.m, self.p)   # binary ontogenetic matrix
        y = np.zeros(2 * self.m)      # gene vector
        x = b*y                  # additive effects vector
        z = x                    # phenotipic values vector TODO incluir ruido
                                   #ambiental com variancia amb
        return {'y': y,'x': x,'z': z,'b': b}


    def mutate (self, ind): 
        """Mutates an ind dict with Individual parameters"""
        for i in range(2 * self.m):
            if (rand.random() < self.mu):
                ind['y'] = ind['y'] + rand.normalvariate(0, self.sigma)
        for i in range(2 * self.m):
            for j in range(2 * self.p):
                if (rand.random() < self.mu_b):
                    ind['b'][i, j] = 1 - ind['b'][i, j]


    def fitness(self, ind, omega, teta):
        """calculates ind fitness from population"""
        delta_s = ind['z'] - teta
        return np.exp(-np.dot(delta_s, np.linalg.solve(omega, delta_s))) 

        
class Population:
    """class for population of Individuals"""
    def __init__(self, n_e, teta, omega, indmod):
        self.indmod = indmod
        self.n_e = n_e
        self.teta = teta
        self.omega = omega
        self.pop = [ indmod.generate() for k in range(n_e) ]
        self.fitness = np.ones(n_e)


    def mutate(self):
        """mutates every individual of popualation"""
        for k in range(self.n_e):
            self.indmod.mutate(self.pop[k])
    

    def update_fitness(self):
        """calculates the fitness of every individual of population"""
        for k in range(self.n_e):
            self.fitness[k] = self.indmod.fitness(self.pop[k], self.omega, self.teta)

   
