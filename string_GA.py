# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:11:57 2019

@author: Pran Regu
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

import numpy as np
import random
import time
import sys

import Multipoint as co
import string_mutate as mu
import Sc_fn2 as sc

# Read input file for intial population
def read_file(file_name):
  smiles_list = []
  with open(file_name,'r') as file:
    for smiles in file:
      smiles = smiles.replace('\n', '')
      smiles_list.append(smiles)

  return smiles_list

def make_initial_population(population_size,file_name):
  smiles_list = read_file(file_name)
  population = []
  for _ in range(population_size):
    smiles = random.choice(smiles_list)
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    string = co.smiles2string(smiles)
    if string:
      population.append(string)
    
  return population

# Scales down scores to a range of 0-1, where a higher number corresponds to a better probability of reproduction
def calculate_normalized_fitness(scores):
  sum_scores = len(scores)-sum(np.absolute(scores))
  sum_scores = sum_scores[0]
  normalized_fitness = [(1-abs(score[0]))/sum_scores for score in scores]
   
  return normalized_fitness

# Creates a mating pool of the specified population size
def make_mating_pool(population,fitness,mating_pool_size):
  mating_pool = []
  for i in range(mating_pool_size):
  	mating_pool.append(np.random.choice(population, p=fitness))

  return mating_pool
 
# Performs crossover and mutation to produce a new molecule
def reproduce(mating_pool,population_size,mutation_rate):
  new_population = []
  while len(new_population) < population_size:
    parent_A = random.choice(mating_pool)
    parent_B = random.choice(mating_pool)
    new_child = co.crossover(parent_A,parent_B)
    if new_child != None:
	    new_child = mu.mutate(new_child,mutation_rate)
	    if new_child != None:
	    	new_population.append(new_child)

  return new_population

# Compares the list of parents and new molecules and chooses the best from both put together
def sanitize(population,scores,population_size):
    prune_population = True 
    if prune_population:
      smiles_list = []
      population_tuples = []
      for score, string in zip(scores,population):
          canonical_smiles = Chem.MolToSmiles(co.string2mol(string))
          if canonical_smiles not in smiles_list:
              smiles_list.append(canonical_smiles)
              population_tuples.append((score, string))
    else:
      population_tuples = list(zip(scores,population))


    population_tuples = sorted(population_tuples, key=lambda x: x[0], reverse=True)[:population_size]
    new_population = [t[1] for t in population_tuples]
    new_scores = [t[0] for t in population_tuples]
    f = open("Result.out", "a")
    print(new_population, file=f)
    print(new_scores,file=f)   
    f.close()
    print(new_population)
    print(new_scores)
    return new_population, new_scores

def GA(args):
  population_size, file_name, generations, mating_pool_size, mutation_rate, \
  scoring_args, max_score = args
 
  maxscore = []  
  population = make_initial_population(population_size,file_name)
  f = open("Result.out", "a")
  print(population, file=f)
  f.close()
  scores = sc.calculate_scores(population,scoring_args, initial)
  fitness = calculate_normalized_fitness(scores)

  for generation in range(generations):
    mating_pool = make_mating_pool(population,fitness,mating_pool_size)
    new_population = reproduce(mating_pool,population_size,mutation_rate)
    f = open("Result.out", "a")
    print(new_population, file=f)
    f.close()
    new_scores = sc.calculate_scores(new_population,scoring_args, generation)
    population, scores = sanitize(population+new_population, scores+new_scores, population_size)  
    fitness = calculate_normalized_fitness(scores)
    
    maxscore.append(scores[0])
    if scores[0][0] <= max_score:
      break

  return (scores, population, generation+1, fitness, maxscore)


# Test case
if __name__ == "__main__":
    Compound = 'C1=CC=CC=C1C2=CC=C(C3=CC=CC=C3)C=C2'
    target = Chem.MolFromSmiles(Compound)
    population_size = 10 
    mating_pool_size = 10
    generations = 15
    mutation_rate = 0.9
    co.average_size = 39.15
    co.size_stdev = 3.50
    #scoring_function = sc.energymin
    max_score = -5 
    #prune_population = True
    scoring_args = [target]
    initial = 50
    file_name = 'terphenyl2.smi'

    (scores, population, generation,fitness, maxscore) = GA([population_size, file_name, generations,
                                           mating_pool_size, mutation_rate, scoring_args, max_score])
                                           
    print('done')

    
    Sl = 0
    for i,j in zip(population, scores):
        Sl = Sl+1
        print(Sl, " Molecule: ", i, "\nScore:",j, "\n") 
        
