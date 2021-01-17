# GA using SMILES finding terphenyl-based molecules maximizing LUMO

# Imports
from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

import numpy as np
import random

import Multipoint as co
import string_mutate as mu
import Sc_fn as sc
import ring_num as nm

# Read file containing inital population list:
def read_file(file_name):
  smiles_list = []
  with open(file_name,'r') as file:
    for smiles in file:
      smiles = smiles.replace('\n', '')
      smiles_list.append(smiles)

  return smiles_list

# Make initial population:
def make_initial_population(population_size,file_name):
  smiles_list = read_file(file_name)
  population = []
  for _ in range(population_size):
    smiles = random.choice(smiles_list)
    population.append(smiles)
    
  return population

# Function to find the degree of conjugation of molecule:
def conjugation(population):
    pc_unsat_list = []
    for molecule in population:
        db_count = 0
        for i in molecule:
            if i == '=' or i == '#':
                db_count+=1
            mol = Chem.MolFromSmiles(molecule)
            bond_count = mol.GetNumBonds()
        pc_unsat_list.append(db_count/bond_count)
    return pc_unsat_list     

# Fitness function: Needs updating! 
def calculate_normalized_fitness(scores, population, wt):
    flat_scores = np.array(scores)
    flat_scores = flat_scores.flatten()
    sum_lumo = sum(np.absolute(flat_scores))
    normalized_scores = [lumo/sum_lumo for lumo in flat_scores]
    pc_unsat_list = conjugation(population)
    normalized_fitness = []
    corrected_scores = []
    for i in range(len(population)):
        fit = (wt*normalized_scores[i]) - ((1-wt)*pc_unsat_list[i])
        corrected_scores.append(fit) 
        
# Scores need normalization because sigma p can't be > 1 in making mating pool
    OldMin = min(corrected_scores)
    OldMax = max(corrected_scores)
    fitness = [(((OldValue - OldMin) * (10 - 1)) / (OldMax - OldMin)) for OldValue in corrected_scores]
    sum_fitness = sum(fitness)
    normalized_fitness = [num/sum_fitness for num in fitness]    
    return normalized_fitness

# Makes mating pool where candidates with better fitness have better odds of appearing:
def make_mating_pool(population, fitness, mating_pool_size):
  mating_pool = []
  for i in range(mating_pool_size):
  	mating_pool.append(np.random.choice(population, p=fitness)) # Sigma p = 1

  return mating_pool
 
# Makes new molecules through crossover and mutation operators:
def reproduce(mating_pool,population_size,mutation_rate):
  new_population = []
  while len(new_population) < population_size:
    parent_A = random.choice(mating_pool)
    parent_B = random.choice(mating_pool)
    new_child = co.crossover(parent_A,parent_B)
    new_child = nm.rnumber(new_child)
    if new_child != None:
	    new_child = mu.mutate(new_child,mutation_rate)
	    if new_child != None:
	    	new_population.append(new_child)

  return new_population

# Compares ith generation with (i-1)th generation and returns the best 20 molecules:
def sanitize(population,scores,population_size):
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
    
    # Prints each generation with respective scores in result.out file
    f = open("Result.out", "a")
    print(new_population, file=f)
    print(new_scores,file=f)   
    f.close()
    
    return new_population, new_scores

# Genetic Algorithm:
def GA(args):
  population_size, file_name, generations, mating_pool_size, mutation_rate, wt = args
 
  population = make_initial_population(population_size,file_name)
  
  # Print initial generation:
  f = open("Result.out", "a")
  print(population, file=f)
  f.close()
  
  scores = sc.calculate_scores(population, initial)
  fitness = calculate_normalized_fitness(scores, population, wt)

  for generation in range(generations):
    mating_pool = make_mating_pool(population,fitness,mating_pool_size)
    new_population = reproduce(mating_pool,population_size,mutation_rate)
    
    # Print each generation
    f = open("Result.out", "a")
    print(new_population, file=f)
    f.close()
    
    new_scores = sc.calculate_scores(new_population, generation)
    population, scores = sanitize(population+new_population, scores+new_scores, population_size)  
    fitness = calculate_normalized_fitness(scores, population, wt)
    
  return (scores, population, generation+1, fitness)


if __name__ == "__main__":
    # Enter GA parameters here:
    population_size = 10 
    mating_pool_size = 10
    generations = 10
    mutation_rate = 0.5
    wt = 0.9
    initial = 50 # Just an identifier for initial population 
    file_name = 'terphenyl2.smi'
    
    # Part of Jenson code: to test if the molecule is legitimate 
    co.average_size = 39.15 
    co.size_stdev = 3.50 
    
    (scores, population, generation, fitness) = GA([population_size, file_name, generations,
                                           mating_pool_size, mutation_rate, wt])
                                           
    print('done')

    
    Sl = 0
    for i,j in zip(population, scores):
        Sl = Sl+1
        print(Sl, " Molecule: ", i, "\nScore:",j, "\n") 
        
