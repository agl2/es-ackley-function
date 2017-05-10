# -*- coding: cp1252 -*-
import random
import numpy as np
from numpy.random import randn
import copy
import Gnuplot

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
            self.mutation_step *= mutation_constant
            if(self.mutation_step < 0.001):
                self.mutation_step = 0.001
                

    def mutation_one_fifth(self, mutation_constant):
        genes_prime = self.genes + self.get_mutation_vector()
        self.num_mutations += 1
        if self.fitness() < self.fitness_(genes_prime):
            self.genes = genes_prime
            self.num_successful_mutations += 1
            
    def mutation_delta_exp (self, delta_mutation, n_genes = 30):
        new_mutation_step = self.mutation_step*np.exp(np.random.normal(0, delta_mutation))
        
        if(new_mutation_step < 0.001):
            new_mutation_step = 0.001
            
        new_genes = self.genes + np.random.normal(0, self.mutation_step, n_genes)
        
        if(self.fitness() < self.fitness_(new_genes)):
           self.mutation_step = new_mutation_step
           self.genes = new_genes

    def fitness(self):
        return -1.*abs(self.f_ackley(self.genes))

    def fitness_(self, genes):
        return -1.*abs(self.f_ackley(genes))

class prec2float(float):
    def __repr__(self):
        return "%0.4f" % self

class EvolutionStrategy:
    def __init__(self, generations=1000, population_size=30, sons_per_iter = 200, mutation_step=1, mutation_constant=0.8, k_parents = 2, delta_mutation = 1, \
                 iter_to_adjust = 30, elitist_suvivor = True, mutation_type = "delta_exp", recombination_type="random"):
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
        self.iter_to_adjust = iter_to_adjust
        
        if(mutation_type == "one_fifth"):
            self.elitist_suvivor = True
        else:
            self.elitist_suvivor = elitist_suvivor
        
        self.mutation_type = mutation_type
        self.recombination_type = recombination_type

    def print_chromossome(self, chromossome):
        print "\tGenes: " + str(map(prec2float, chromossome.genes))
        print "\tPasso de mutação: " + str(prec2float(chromossome.mutation_step))
        print "\tAckley Value: %.4f" % self.f_ackley(chromossome.genes)
        if(self.mutation_type == "one_fifth"):
            print "\tTaxa de sucesso: " + str(chromossome.get_success_probability())

    def print_population(self):
        for i in range(len(self.population)):
            print "Indivíduo: " + str(i)
            self.print_chromossome(population[i])

    def init_population(self):
        for i in range(self.population_size):
            self.population.append(Chromossome(genes = 30*np.random.random(30)-15, mutation_step = self.mutation_step))

    def apply_mutation(self):
        for i in range(len(self.population)):
            if(self.mutation_type == "delta_exp"):
                self.population[i].mutation_delta_exp(delta_mutation = self.delta_mutation)
            elif(self.mutation_type == "one_fifth"):
                self.population[i].mutation_one_fifth(mutation_constant = self.mutation_constant)
                if(self.count_adjust > self.iter_to_adjust):
                    self.population[i].adjust_mutation_step(self.mutation_constant, self.success_rate)
                    self.count_adjust = 0
                else:
                    self.count_adjust += 1

    def parent_selection(self, k_parents = 2):
        return random.sample(self.population, k_parents)

    def apply_recombination(self):
        if(self.recombination_type == "mean"):
            self.mean_recombination()
        elif(self.recombination_type == "random"):
            self.random_recombination()
    
    def mean_recombination(self):
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


        if(self.elitist_suvivor):
            for chromossome in new_population:
                self.population.append(chromossome)
        else:
            self.population = new_population

    def random_recombination(self):
        new_population = []
        for i in range(self.sons_per_iter):
            parents = self.parent_selection(k_parents = self.k_parents)
            genes_son = []
            for i in range(len(parents[0].genes)):
                parent_select = np.random.randint(0,self.k_parents)
                genes_son.append(parents[parent_select].genes[i])
                
            parent_select = np.random.randint(0,self.k_parents)
            mutation_step_son = parents[parent_select].mutation_step          
            new_population.append(Chromossome(genes = genes_son, mutation_step = mutation_step_son))
        if(self.elitist_suvivor):
            for chromossome in new_population:
                self.population.append(chromossome)
        else:
            self.population = new_population

            

    def apply_selection(self):
        self.population.sort(key=lambda chromossome : chromossome.fitness(), reverse=True)	
        self.population = self.population[:self.population_size]

    def run(self, verbose=0):
        self.verbose = verbose
        self.init_population()
        self.apply_selection()
        gen = 0
        history = [self.population[0]]
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
            history.append(self.population[0])
        return history

if __name__ == '__main__':
    es = EvolutionStrategy()
    hs = es.run(verbose = 0)
    print "First: "
    es.print_chromossome(hs[0])
    print "Best: "
    es.print_chromossome(hs[-1])

    gplt = Gnuplot.Gnuplot(debug=1)

    gplt.reset()
    maximo = []
    passo = []
    for chrm in hs:
        maximo.append(chrm.fitness())
        passo.append(chrm.mutation_step)
       
    #mediaData = Gnuplot.Data(media, title='Fitness Médio')
    #minData = Gnuplot.Data(minimo, title='Fitness Mínimo')
    maxData = Gnuplot.Data(maximo, title='Fitness Máximo')
    passoData = Gnuplot.Data(passo, title='Passo do cromossomo')

    gplt('set data style lines')
    gplt.xlabel('Geracao')
    gplt.ylabel('Fitness')
        
    gplt.plot(maxData)
    raw_input('Please press return to continue...\n')

    gplt.reset()
    gplt('set data style lines')
    gplt.xlabel('Geracao')
    gplt.ylabel('Passo')
        
    gplt.plot(passoData)
    raw_input('Please press return to continue...\n')
