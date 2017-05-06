# -*- coding: cp1252 -*-
import random
import Gnuplot
import math
import numpy as np


class Cromossomo:

    def __init__(self, genes, pai1, pai2, geracao, gerador):
        self.genes = genes
        self.fitness = fitness(genes)
        self.pai1 = pai1
        self.pai2 = pai2
        self.geracao = geracao
        self.gerador = gerador
        
        
def iniciar_populacao():
    final = [30*np.random.random(30)-15 for i in range(100)]
    for i in range(len(final)):
        final[i] = Cromossomo(final[i], None, None, 0, 'I')
    return final        
        
def vetor_mutacao(passo):
    return np.random.normal(0, passo, 30)

def probabilidade_sucesso(m_sucesso, num_mutacoes):
    return m_sucesso/num_mutacoes

class Estrategia:# tem que iniciar a estrategia evolutiva
    def __init__(self, passo_mutacao, ajuste):#passo_mutacao = 0.1, ajuste = 0.8
        self.passo_mutacao = passo_mutacao
        self.ajuste = ajuste
        self.verbose = 0
        self.num_mutacoes = 0
        self.m_sucesso = 0
        self.taxa_sucesso = 0.2
        
    
    def ajustar_passo(self):
        ps = probabilidade_sucesso(self.m_sucesso, self.num_mutacoes)
        if ps > self.taxa_sucesso:
            self.passo_mutacao /= self.ajuste
        elif ps < self.taxa_sucesso:
            self.passo_mutacao *= self.ajuste

#o calculo do fitness eh feito pelo valor da funcao de ackley
def fitness(genes):
    sum1 = 0.0
    sum2 = 0.0    
    for g in genes:
        sum1 += g**2.0
        sum2 += math.cos(2.0*math.pi*g)
    n = float(len(genes))
    return -20.0*math.exp(-0.2*math.sqrt(sum1/n)) - math.exp(sum2/n) + 20.0 + math.e


#a populacao eh iniciada com um vetor que tem valores de float entre -15 e 15

    


#ao chamar essa funcao tem que contar o numero de mutacoes e o numero de mutacoes bem sucessedidas
def mutacao(gene, passo):
    gene_final = gene + vetor_mutacao(passo)
    if (gene.fitness) < (fitness(gene_final)):# se a mutacao foi bem sucessedida retorna 1
        gene.genes = gene_final
        return 1
    else:#se nao retorna 0
        return 0

   
def main():
    
    filho1 = []
    filho2 = []
    i = 0

    #inicia a populacao com 100 individuos e os ordena de acordo com o fitness
    populacao = iniciar_populacao()
    #ordena os pais de acordo com o fitness
    populacao.sort(key=lambda p : p.fitness)
    print(populacao[0].fitness)
    print(populacao[-1].fitness)

    media = []
    maximo = []
    minimo = []
    
    #a partir daqui tem que ser feito o laco para poder tentar encontrar a solucao para o problema
    #while(i <= 100000):
    
       


if __name__ == '__main__':
    main()
