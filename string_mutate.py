# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 18:42:06 2019

@author: Pran Regu

MUTATION OPERATOR
"""

from rdkit import Chem
from rdkit.Chem import AllChem

import random
import numpy as np

import Multipoint as co

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

# List of possible mutations
def get_symbols():
    symbols = ['C(C#N)', 'C(C=O)', 'C(C(=O)C)', 'C(O)', 'C(C4=CC=CC=C4)', 'C(C(=O)O)'] 
               
    return symbols

# Mutation operator; Causes a mutation if random number less than specified mutation rate in string_GA.oy
def mutate(child,mutation_rate):
    if random.random() > mutation_rate:
        return child
    symbols = get_symbols()
    child = co.string2list(child)
   
    for i in range(50):
        random_number = random.random()
        mutated_gene = random.randint(0, len(child) - 1)
        random_symbol_number = random.randint(0, len(symbols)-1)
        new_child = list(child)
        random_number = random.random()
        
        # Finds a new position for mutation if the chosen position is a bond
        if (new_child[mutated_gene] == '='):
            m=1 #nothing
        
        else:
            new_child[mutated_gene] = symbols[random_symbol_number]
            new_child = co.list2string(new_child)
            if co.string_OK(new_child):
                 return new_child

    return co.list2string(child)
           
# Test case    
if __name__ == "__main__":
    co.average_size = 39.15
    co.size_stdev = 3.50
    mutation_rate = 1.0
    co.string_type = 'smiles'
    string = 'C1=CC=CC=C1C2=CC=C(C3=CC=CC=C3)C=C2'
    child = mutate(string,mutation_rate)
    print(child)
