import random
import numpy as np
from numpy.random import randn
import copy

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
    def __init__(self, generations=10000, population_size=20, sons_per_iter = 100, mutation_step=0.1, adjust_mutation_constant=0.8):
        self.generations = generations
        self.population_size = population_size
        self.sons_per_iter = sons_per_iter
        self.population = []
        self.f_ackley = Ackley().f_x
        self.mutation_step = mutation_step
        self.adjust_mutation_constant = adjust_mutation_constant
        self.success_rate = .2
        self.num_mutations = 0
        self.num_successful_mutations = 0
        self.verbose = 0

    def print_population(self):
        for i in range(len(self.population)):
            print self.population[i]

    def init_population(self):
        for i in range(self.population_size):
            cromossome = 30*np.random.random(30)-15
            self.population.append(cromossome)

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
        for i in range(len(self.population)):
            cromossome_prime = self.population[i] + self.get_mutation_vector()
            self.num_mutations += 1
            if self.fitness(self.population[i]) < self.fitness(cromossome_prime):
                self.population[i] = cromossome_prime
                self.num_successful_mutations += 1
            self.adjust_mutation_step()
        
    def parent_selection(self, k_parents = 2):
        return random.sample(self.population, k_parents)

    def apply_recombination(self):
        new_population = []
        for i in range(self.sons_per_iter):
            parents = parent_selection()
            cromossome_son = np.array(30*[0])
            for cromossome in parents:
                cromossome_son += cromossome
            cromossome_son = cromossome_son/len(parents)
            new_population.append(cromossome_son)

        self.population = new_population

    def apply_selection(self):
        self.population.sort(key=lambda cromossome : self.fitness(cromossome))	
        self.population = self.population[:self.population_size]

    def run(self, verbose=0):
        self.verbose = verbose
        self.init_population()
        self.apply_selection()
        gen = 0
        history = [(self.population[0], self.f_ackley(self.population[0]))]
        if self.verbose == 1:
            print "gen: %d" % gen
            self.print_population()
            print "Ackley(x): %.5f" % self.f_ackley(self.population[0])
        while gen < self.generations:
            if self.verbose == 2:
                print "gen: %d" % gen
            gen += 1
            self.apply_recombination
            self.apply_mutation()
            self.apply_selection()
            if self.verbose == 1:
                print "gen: %d" % gen
                self.print_population()
                print "Ackley(x): %.5f" % self.f_ackley(self.population[0])
            history.append((self.population[0],  self.f_ackley(self.population[0])))
        return history

es = EvolutionStrategy()
hs = es.run(verbose = 2)
print "First: " + str(hs[0])
print "Best: " + str(hs[-1])
