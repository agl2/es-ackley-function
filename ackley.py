import random
import numpy as np
from numpy.random import randn

class Ackley:
    def __init__(self, N=30, C1=20, C2=.2, C3=2*np.pi):
        self.N = N
        self.c1 = C1
        self.c2 = C2
        self.c3 = C3

    def f_x(self, x):
        part1 = -1. * self.c1 * np.exp(
            -1. * self.c2 * np.sqrt((1./self.N) * sum(map(lambda nb: nb**2, x)))
            )
        part2 = -1. * np.exp(
            (1./self.N) * \
            sum(map(lambda nb: np.cos(self.c3 * nb), x))
            )
        return part1 + part2 + self.c1 + np.exp(1)

class EvolutionStrategy:
    def __init__(self, generations=5000, population_size=5, mutation_step=0.1, adjust_mutation_constant=0.8):
        self.generations = generations
        self.population_size = population_size
        self.population = []
        self.f_ackley = Ackley().f_x
        self.cromossome = None
        self.mutation_step = mutation_step
        self.adjust_mutation_constant = adjust_mutation_constant
        self.success_rate = .2
        self.num_mutations = 0
        self.num_successful_mutations = 0
        self.verbose = 0

    def print_cromossome(self):
        for i in range(self.population_size):
            print self.population[i]

    def init_cromossome(self):
        for i in range(self.population_size):
            self.cromossome = 30*np.random.random(30)-15
            self.population.append(self.cromossome)

    def get_mutation_vector(self):
        return np.random.normal(0, self.mutation_step, 30)

    def get_success_probability(self):
        return self.num_successful_mutations / float(self.num_mutations)


    def fitness(self, cromossome):
        return -1.*abs(self.f_ackley(cromossome))

    def adjust_mutation_step(self):
        ps = self.get_success_probability()
        if self.verbose == 1:
            print "ps: %.4f" % ps
        if ps > self.success_rate:
            self.mutation_step /= self.adjust_mutation_constant
        elif ps < self.success_rate:
            self.mutation_step *= self.adjust_mutation_constant
        if self.verbose == 1:
            print "mutation_step: %.4f" % self.mutation_step

    def apply_mutation(self):
        cromossome_prime = self.cromossome + self.get_mutation_vector()
        self.num_mutations += 1
        if self.fitness(self.cromossome) < self.fitness(cromossome_prime):
            self.cromossome = cromossome_prime
            self.num_successful_mutations += 1
        self.adjust_mutation_step()
        
    def parent_selection(self, k_parents = 2):
        return random.sample(self.population, k_parents)

    def run(self, verbose=0):
        self.verbose = verbose
        self.init_cromossome()
        gen = 0
        history = [(self.cromossome, self.f_ackley(self.cromossome))]
        if self.verbose == 1:
            print "gen: %d" % gen
            self.print_cromossome()
            print "Ackley(x): %.5f" % self.f_ackley(self.cromossome)
        while gen < self.generations:
            gen += 1
            #aqui aplicar a recombinacao
            self.apply_mutation()
            if self.verbose == 1:
                print "gen: %d" % gen
                self.print_cromossome()
                print "Ackley(x): %.5f" % self.f_ackley(self.cromossome)
            history.append((self.cromossome, self.f_ackley(self.cromossome)))
        return history

es = EvolutionStrategy()
es.run(verbose = 1)