# -*- coding: cp1252 -*-
import random
import numpy as np
from numpy.random import randn
import copy
import Gnuplot
import os,sys

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
    def __init__(self, genes = np.array(30*[0]), mutation_step = 1, num_mutations = 0, num_successful_mutations = 0):
        self.genes = genes
        self.mutation_step = mutation_step
        self.num_mutations = num_mutations
        self.num_successful_mutations = num_successful_mutations
        self.f_ackley = Ackley().f_x

    def get_mutation_vector(self, n_genes = 30):
        return np.random.normal(0, self.mutation_step, n_genes)

    def get_success_probability(self):
        if(self.num_mutations == 0):
            return 0;
        return float(self.num_successful_mutations) / float(self.num_mutations)

    def adjust_mutation_step(self, mutation_constant, success_rate):
        ps = self.get_success_probability()
        self.num_mutations = 0
        self.num_successful_mutations = 0
        if ps > success_rate:
            self.mutation_step /= mutation_constant
        elif ps < success_rate:
            self.mutation_step *= mutation_constant
            if(self.mutation_step  <  1e-15):
                self.mutation_step  = 1e-15
                

    def mutation_one_fifth(self, mutation_constant):
        genes_prime = self.genes + self.get_mutation_vector()
        self.num_mutations += 1
        if self.fitness() < self.fitness_(genes_prime):
            self.genes = genes_prime
            self.num_successful_mutations += 1
            
            
    def mutation_delta_exp (self, delta_mutation, n_genes = 30):
        new_mutation_step = self.mutation_step*np.exp(np.random.normal(0, delta_mutation))
        

        if(new_mutation_step <  1e-15):
            new_mutation_step = 1e-15
        elif(new_mutation_step > 1e10):
            new_mutation_step = 1e10
            
        new_genes = self.genes + np.random.normal(0, self.mutation_step, n_genes)
        self.mutation_step = new_mutation_step
        self.genes = new_genes

    def fitness(self):
        return -1.*abs(self.f_ackley(self.genes))

    def fitness_(self, genes):
        return -1.*abs(self.f_ackley(genes))

class sciNot(float):
    def __repr__(self):
        return "%0.4e" % self

class EvolutionStrategy:
    def __init__(self, generations=10000, population_size=30, sons_per_iter = 200, init_mutation_step=1000, mutation_constant=0.95, delta_mutation = 1, \
                 iter_to_adjust = 5, elitist_suvivor = False, mutation_type = "delta_exp", recombination_type="random", parent_sel = "global"):
        self.generations = generations
        self.population_size = population_size
        self.sons_per_iter = sons_per_iter
        self.population = []
        self.init_mutation_step = init_mutation_step
        self.mutation_constant = mutation_constant
        self.success_rate = .2
        self.verbose = 0
        self.parent_sel = parent_sel
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
        print "\tGenes: " + str(map(sciNot, chromossome.genes))
        print "\tPasso de mutacao: " + str(sciNot(chromossome.mutation_step))
        print "\tAckley Value: " +  str(sciNot(self.f_ackley(chromossome.genes)))
        if(self.mutation_type == "one_fifth"):
            print "\tTaxa de sucesso: " + str(chromossome.get_success_probability())

    def print_population(self):
        for i in range(len(self.population)):
            print "Individuo: " + str(i)
            self.print_chromossome(self.population[i])

    def init_population(self):
        self.population = []
        self.count_adjust = 0
        for i in range(self.population_size):
            self.population.append(Chromossome(genes = 30*np.random.random(30)-15, mutation_step = self.init_mutation_step*np.random.random()))

    def apply_mutation(self):
        if(self.mutation_type == "one_fifth"):
             self.count_adjust += 1
             
        for i in range(len(self.population)):

            if(self.mutation_type == "delta_exp"):
                self.population[i].mutation_delta_exp(delta_mutation = self.delta_mutation)

            elif(self.mutation_type == "one_fifth"):
                self.population[i].mutation_one_fifth(mutation_constant = self.mutation_constant)
                if(self.count_adjust == self.iter_to_adjust):
                    self.population[i].adjust_mutation_step(self.mutation_constant, self.success_rate)

        if(self.count_adjust == self.iter_to_adjust):
            self.count_adjust = 0
    #applyMutation


    def parent_selection(self):
        return random.sample(self.population, 2)

    def apply_recombination(self):
        if(self.recombination_type == "mean"):
            self.mean_recombination()
        elif(self.recombination_type == "random"):
            self.random_recombination()
                   
    def mean_recombination(self):
        new_population = []
        for sons in range(self.sons_per_iter):
            if(self.parent_sel == "local"):
                parents = self.parent_selection()
            genes_son = []
            for i in range(30):
                if(self.parent_sel == "global"):
                    parents = self.parent_selection()
                    
                genes_son.append((parents[0].genes[i] + parents[1].genes[i])/2)

            
            if(self.parent_sel == "global"):
                parents = self.parent_selection()                
            if(self.mutation_type == "delta_exp"):
                mutation_step_son = (parents[0].mutation_step + parents[1].mutation_step)/2
                new_population.append(Chromossome(genes = genes_son, mutation_step = mutation_step_son))
            else:
                mutation_step_son = (parents[0].mutation_step + parents[1].mutation_step)/2
                num_mutations_son = parents[0].num_mutations
                num_successful_mutations_son = (parents[0].num_successful_mutations + parents[1].num_successful_mutations)/2
                new_population.append(Chromossome(genes = genes_son, mutation_step = mutation_step_son, \
                                                  num_successful_mutations = num_successful_mutations_son, \
                                                  num_mutations = num_mutations_son
                                                  ))
                
        if(self.elitist_suvivor):
            for chromossome in new_population:
                self.population.append(chromossome)
        else:
            self.population = new_population

               
    def random_recombination(self):
        new_population = []
        for sons in range(self.sons_per_iter):
            if(self.parent_sel == "local"):
                parents = self.parent_selection()
            genes_son = []
            for i in range(30):
                if(self.parent_sel == "global"):
                    parents = self.parent_selection()
                parent_select = np.random.randint(0,2)
                genes_son.append(parents[parent_select].genes[i])
            #endFor
            
            if(self.parent_sel == "global"):
                parents = self.parent_selection()                
            if(self.mutation_type == "delta_exp"):
                parent_select = np.random.randint(0,2)
                mutation_step_son = parents[parent_select].mutation_step
                new_population.append(Chromossome(genes = genes_son, mutation_step = mutation_step_son))
            else:
                mutation_step_son = parents[parent_select].mutation_step
                num_mutations_son = parents[0].num_mutations
                num_successful_mutations_son = parents[parent_select].num_successful_mutations
                new_population.append(Chromossome(genes = genes_son, mutation_step = mutation_step_son, \
                                                  num_successful_mutations = num_successful_mutations_son, \
                                                  num_mutations = num_mutations_son
                                                  ))
        #endFor
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
        bests = [self.population[0]]
        means = [np.mean( [cr.fitness() for cr in self.population] )]
        variances = [np.var( [cr.fitness() for cr in self.population] )]
        mutation_st_mean = [np.mean( [cr.mutation_step for cr in self.population] )]
        mutation_st_var = [np.var( [cr.mutation_step for cr in self.population] )]
        
        if self.verbose == 1:
            print "=========================================================="
            print "GERAÇÃO: %d" % gen
            #self.print_population()
            print "Ackley(x): %.5f" % self.f_ackley(self.population[0].genes)
        while gen < self.generations:
            gen += 1
            self.apply_recombination()
            self.apply_mutation()
            self.apply_selection()
            if self.verbose == 1:
                print "==========================================================="
                print "GERAÇÃO: %d" % gen
                #self.print_population()
                print "Ackley(x) of best: %.4f" % self.f_ackley(self.population[0].genes)
                print "Ackley(x) of worst: %.4f" % self.f_ackley(self.population[-1].genes)
                print "Population length: %d" % len(self.population)
            bests.append(self.population[0])
            means.append(np.mean( [cr.fitness() for cr in self.population] ))
            variances.append(np.var( [cr.fitness() for cr in self.population] ))
            mutation_st_mean.append(np.mean( [cr.mutation_step for cr in self.population] ))
            mutation_st_var.append(np.var( [cr.mutation_step for cr in self.population] ))
        return (bests, means, variances, mutation_st_mean, mutation_st_var)

def main():
    # 
    #
    # recombination_type:   random                  ;               mean
    # parent_sel:           global                  ;               local
    # mutation_type:        one_fifth               ;               delta_exp 
    # elitist_suvivor:      True (lambda + delta)   ;               False (lambda, delta)
    
    
    #es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 200, mutation_step=1, mutation_constant=0.95, delta_mutation = 50, \
    #             iter_to_adjust = 5, elitist_suvivor = False, mutation_type = "one_fifth", recombination_type="random", parent_sel = "local")
    #es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 20, mutation_step=1, mutation_constant=0.95, delta_mutation = 0.18, \
    #    iter_to_adjust = 5, elitist_suvivor = True, mutation_type = "one_fifth", recombination_type="random", parent_sel = "local")    

    es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 200,\
                           init_mutation_step=100, mutation_constant=0.95, delta_mutation = 1.8, \
                           iter_to_adjust = 5, elitist_suvivor = False, mutation_type = "delta_exp",\
                           recombination_type="mean", parent_sel = "global")

    f = open('out.txt', 'w')
    sys.stdout = f   
    doTest(es, 'graph')

 
    
def doTest(es, out_name):
    ms = np.array((es.generations+1)*[0.])
    vs = np.array((es.generations+1)*[0.])
    msm = np.array((es.generations+1)*[0.])
    msv = np.array((es.generations+1)*[0.]) 
    for i in range(10):
        (hs, ms_temp, vs_temp, msm_temp, msv_temp) = es.run(verbose = 0)
        print "Melhor Inicial: " + str(i+1)
        es.print_chromossome(hs[0])
        print "Melhor Final: " + str(i+1)
        es.print_chromossome(hs[-1])

        ms += ms_temp
        vs += vs_temp
        msm += msm_temp
        msv += msv_temp

    ms /= 10.
    vs /= 10.
    msm /= 10.
    msv /= 10.



    gplt = Gnuplot.Gnuplot(debug=1)

    maximo = []
    for chrm in hs:
        maximo.append(chrm.fitness())

    maxData = Gnuplot.Data(maximo, title='-|F(X|)')
    mediaData = Gnuplot.Data(ms, title='Média')
    varData = Gnuplot.Data(vs, title='Variancia')
    mediaPassoData = Gnuplot.Data(msm, title='Média passo mutação')
    varPassoData = Gnuplot.Data(msv, title='Variância passo mutação')
       
    maxDataLog = Gnuplot.Data(map(np.log10, map(abs,maximo)), title='log10 |F(X|')
    mediaDataLog = Gnuplot.Data(map(np.log10, map(abs,ms)), title='log10 |Média|')
    varDataLog = Gnuplot.Data(map(np.log10, map(abs,vs)), title='log10 |Variância|')
    mediaPassoDataLog = Gnuplot.Data(map(np.log10, map(abs,msm)), title='log10 Média passo mutação')
    varPassoDataLog = Gnuplot.Data(map(np.log10, map(abs,msv)), title='log10 Var passo mutação')


    gplt('set data style lines')
    gplt.xlabel('Geracao')
    gplt.ylabel('Fitness log10')
    gplt('set terminal png size 1080,720 enhanced font "Helvetica,20"')
    gplt('set output "logGraph/' + out_name + '1_log.png"')
    gplt.plot(mediaDataLog, varDataLog)

    gplt.xlabel('Geracao')
    gplt.ylabel('Passo mutacao log10')
    gplt('set terminal png size 1080,720 enhanced font "Helvetica,20"')
    gplt('set output "logGraph/' + out_name + '2_log.png"')
    gplt.plot(mediaPassoDataLog, varPassoDataLog)

    gplt('set data style lines')
    gplt.xlabel('Geracao')
    gplt.ylabel('Fitness')
    gplt('set terminal png size 1080,720 enhanced font "Helvetica,20"')
    gplt('set output "graph/' + out_name + '1.png"')
    gplt.plot(mediaData, varData)

    gplt('set data style lines')
    gplt.xlabel('Geracao')
    gplt.ylabel('Passo Mutacao')
    gplt('set terminal png size 1080,720 enhanced font "Helvetica,20"')
    gplt('set output "graph/' + out_name + '2.png"')
    gplt.plot(mediaPassoData, varPassoData) 

def test():

    for init_mutation_step in [10]:
        for delta_mutation in [0.018]:
            for elitist_suvivor in [False]:
                for  recombination_type in ['mean']:
                    for parent_sel in ['local', 'global']:
                        output_name = 'delta_exp' + '_dm' + str(delta_mutation) + \
                                    '_ms' + str(init_mutation_step) + '_el' + str(elitist_suvivor) + '_' + recombination_type + '_' + parent_sel
                        f = open('data/' + output_name + '.txt','w')
                        sys.stdout = f
                           
                        es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 200,\
                        delta_mutation = delta_mutation, elitist_suvivor = elitist_suvivor, mutation_type = "delta_exp",\
                        init_mutation_step = init_mutation_step, recombination_type=recombination_type, parent_sel = parent_sel)

                        doTest(es, output_name)
                        
    for init_mutation_step in [100]:
        for delta_mutation in [0.018]:
            for elitist_suvivor in [True, False]:
                for  recombination_type in ['random', 'mean']:
                    for parent_sel in ['local', 'global']:
                        output_name = 'delta_exp' + '_dm' + str(delta_mutation) + \
                                    '_ms' + str(init_mutation_step) + '_el' + str(elitist_suvivor) + '_' + recombination_type + '_' + parent_sel
                        f = open('data/' + output_name + '.txt','w')
                        sys.stdout = f
                           
                        es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 200,\
                        delta_mutation = delta_mutation, elitist_suvivor = elitist_suvivor, mutation_type = "delta_exp",\
                        init_mutation_step = init_mutation_step, recombination_type=recombination_type, parent_sel = parent_sel)

                        doTest(es, output_name)
                        
    for init_mutation_step in [10, 100]:
        for delta_mutation in [0.18, 1.8]:
            for elitist_suvivor in [True, False]:
                for  recombination_type in ['random', 'mean']:
                    for parent_sel in ['local', 'global']:
                        output_name = 'delta_exp' + '_dm' + str(delta_mutation) + \
                                    '_ms' + str(init_mutation_step) + '_el' + str(elitist_suvivor) + '_' + recombination_type + '_' + parent_sel
                        f = open('data/' + output_name + '.txt','w')
                        sys.stdout = f
                           
                        es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 200,\
                        delta_mutation = delta_mutation, elitist_suvivor = elitist_suvivor, mutation_type = "delta_exp",\
                        init_mutation_step = init_mutation_step, recombination_type=recombination_type, parent_sel = parent_sel)

                        doTest(es, output_name)



def test2():
    for init_mutation_step in [10, 100]:
        for mutation_constant in [0.8, 0.9, 0.98]:
            for  recombination_type in ['random', 'mean']:
                for parent_sel in ['local', 'global']:
                    output_name = 'one_fifth' + '_spi' + '_ms' + str(init_mutation_step) + \
                        '_mc' + str(mutation_constant) + '_' + recombination_type + '_' + parent_sel

                    f = open('data/' + output_name + '.txt','w')
                    sys.stdout = f
                           
                    es = EvolutionStrategy(generations=1000, population_size=30, sons_per_iter = 200,\
                        init_mutation_step = init_mutation_step, elitist_suvivor = True, mutation_type = "one_fifth",\
                        recombination_type=recombination_type, parent_sel = parent_sel, iter_to_adjust = 5)

                    doTest(es, output_name)
                        
    
if __name__ == '__main__':
    main()



    
    
