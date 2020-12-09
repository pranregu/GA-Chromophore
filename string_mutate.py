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

def get_symbols():
    symbols = ['C(N(C)C)', 'C(N)', 'C(N(C))', 'C(O)', 'C(OC)', 'C(C)', 'C(F)', 
               'C(Cl)', 'C(SC)', 'C(Br)', 'C(C=O)', 'C(C(F)(F)F)', 'C(C#N)', 
               'C([N+](=O)[O-])', 'C(C3=CC=C4C=CC=CC4=C3)' ]
    return symbols

def mutate(child,mutation_rate):
    if random.random() > mutation_rate:
        return child
    symbols = get_symbols()
    child = co.string2list(child)
   
    for i in range(50):
        
        mutated_gene = random.randint(0, len(child) - 1)
        random_symbol_number = random.randint(0, len(symbols)-1)
        new_child = list(child)
        random_number = random.random()
        
        if new_child[mutated_gene] == 'C':
            if new_child[mutated_gene+1] == 'C':
                new_child[mutated_gene] = symbols[random_symbol_number]
                new_child = co.list2string(new_child)
                if co.string_OK(new_child):
                    return new_child
            elif new_child[mutated_gene+1] == '=' and new_child[mutated_gene-1] == 'C':     
                new_child[mutated_gene] = symbols[random_symbol_number]
                new_child = co.list2string(new_child)
                if co.string_OK(new_child):
                    return new_child
            elif new_child[mutated_gene+1] == '(':
                a = len(child)
                count = 0
                
                for j in range(mutated_gene, a):
                    if child[j] == '(':
                        count = count+1
                    elif child[j] == ')': 
                        count = count-1
                    if child[j] == ')' and count == 0:  
                        break
                #j = j+1
                new_child = co.list2string(new_child)
                str1 = new_child[mutated_gene : j+1]
                #print(str1)
                str2 =symbols[random_symbol_number]
                #print(str2)
                if len(str1)<13:
                    new_child = new_child.replace(str1, str2, 1)
                if co.string_OK(new_child):
                    return new_child 

    return co.list2string(child)
                    
if __name__ == "__main__":
    co.average_size = 39.15
    co.size_stdev = 3.50
    mutation_rate = 1.0
    co.string_type = 'smiles'
    string = 'C1=CC=CC=C1C2=CC=C(C3=CC=CC=C3)C=C2'
    child = mutate(string,mutation_rate)
    print(child)
