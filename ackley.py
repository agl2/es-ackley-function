# -*- coding: cp1252 -*-
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

class Chromossome:
    def __init__(self, genes = np.array(30*[0]), mutation_step = 1):
        self.genes = genes
        self.mutation_step = mutation_step
        self.num_mutations = 0
        self.num_successful_mutations = 0
        self.f_ackley = Ackley().f_x

    def get_mutation_vector(self, n_genes = 30):
        return np.random.normal(0, self.mutation_step, n_genes)

    def get_success_probability(self):
        if(self.num_mutations == 0):
            return 0;
        return float(self.num_successful_mutations) / float(self.num_mutations)

    def adjust_mutation_step(self, mutation_constant, success_rate):
        ps = self.get_success_probability()       
        if ps > success_rate:
            self.mutation_step /= mutation_constant
        elif ps < success_rate:
            if(self.mutation_step > 0.001):
                self.mutation_step *= mutation_constant
                

    def mutation_one_fifth(self, mutation_constant):
        genes_prime = self.genes + self.get_mutation_vector()
        self.num_mutations += 1
        if self.fitness() < self.fitness_(genes_prime):
            self.genes = genes_prime
            self.num_successful_mutations += 1

    #Essa mutacao nao funciona, normalmente aumenta o fitness    
    def mutation_delta_exp (self, delta_mutation, n_genes = 30):
        self.mutation_step *= np.exp(np.random.normal(0, delta_mutation))
        self.genes += np.random.normal(0, self.mutation_step, n_genes)

    def fitness(self):
        return -1.*abs(self.f_ackley(self.genes))

    def fitness_(self, genes):
        return -1.*abs(self.f_ackley(genes))

class prec2float(float):
    def __repr__(self):
        return "%0.2f" % self

class EvolutionStrategy:
    def __init__(self, generations=10000, population_size=30, sons_per_iter = 100, mutation_step=1, mutation_constant=0.8, k_parents = 2, delta_mutation = 0.1826):
        self.generations = generations
        self.population_size = population_size
        self.sons_per_iter = sons_per_iter
        self.population = []
        self.mutation_step = mutation_step
        self.mutation_constant = mutation_constant
        self.success_rate = .2
        self.verbose = 0
        self.k_parents = k_parents
        self.f_ackley = Ackley().f_x
        self.count_adjust = 0
        self.delta_mutation = delta_mutation

    def print_population(self):
        for i in range(len(self.population)):
            print "Gene " + str(i) + ": " + str(map(prec2float, self.population[i].genes))
            print "Passo de mutação: " + str(prec2float(self.population[i].mutation_step))
            print "Taxa de sucesso: " + str(self.population[i].get_success_probability())

    def init_population(self):
        for i in range(self.population_size):
            self.population.append(Chromossome(genes = 30*np.random.random(30)-15, mutation_step = self.mutation_step))

    def apply_mutation(self):
        for i in range(len(self.population)):
            #self.population[i].mutation_delta_exp(delta_mutation = self.delta_mutation)
            self.population[i].mutation_one_fifth(self.mutation_constant)
            if(self.count_adjust > 500):
                self.population[i].adjust_mutation_step(self.mutation_constant, self.success_rate)
                self.count_adjust = 0
            else:
                self.count_adjust += 1

    def parent_selection(self, k_parents = 2):
        return random.sample(self.population, k_parents)

    def apply_recombination(self):
        new_population = []
        for i in range(self.sons_per_iter):
            parents = self.parent_selection(k_parents = self.k_parents)
            genes_son = np.array(30*[0.])
            mutation_step_son = 0
            for chromossome in parents:
                genes_son += chromossome.genes
                mutation_step_son += chromossome.mutation_step
            genes_son /= len(parents)
            mutation_step_son /= len(parents)           
            new_population.append(Chromossome(genes = genes_son, mutation_step = mutation_step_son))

        for chromossome in new_population:
            self.population.append(chromossome)
        #if(self.sons_per_iter > 0):
            #self.population = new_population
       

    def apply_selection(self):
        self.population.sort(key=lambda chromossome : chromossome.fitness(), reverse=True)	
        self.population = self.population[:self.population_size]

    def run(self, verbose=0):
        self.verbose = verbose
        self.init_population()
        self.apply_selection()
        gen = 0
        history = [(self.population[0].genes, self.f_ackley(self.population[0].genes))]
        if self.verbose == 1:
            print "=========================================================="
            print "GERAÇÃO: %d" % gen
            #self.print_population()
            print "Ackley(x): %.5f" % self.f_ackley(self.population[0].genes)
        while gen < self.generations:
            gen += 1
            self.apply_mutation()
            self.apply_recombination()
            self.apply_selection()
            if self.verbose == 1:
                print "==========================================================="
                print "GERAÇÃO: %d" % gen
                #self.print_population()
                print "Ackley(x) of best: %.4f" % self.f_ackley(self.population[0].genes)
                print "Ackley(x) of worst: %.4f" % self.f_ackley(self.population[-1].genes)
                print "Population length: %d" % len(self.population)
            history.append((self.population[0].genes,  self.f_ackley(self.population[0].genes)))
        return history

es = EvolutionStrategy()
hs = es.run(verbose = 0)
print "First: " + str(hs[0])
print "Best: " + str(hs[-1])
