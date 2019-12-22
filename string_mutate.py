# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 18:42:06 2019

@author: Pran Regu
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

import random
import numpy as np

import Multipoint as co

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

def get_symbols():
    symbols = ['C(C#N)', 'C(C=O)', 'C(C(=O)C)', 'C(O)', 'C(C4=CC=CC=C4)', 'C(C(=O)O)', 
               'C(C4=CC=C(C5=CC=CC=C5)C=C4)', 'C(C4=C(C5=CC=CC=C5)C=CC=C4)', 'C(C4=CC(C5=CC=CC=C5)=CC=C4)']
    return symbols

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
        
        if (new_child[mutated_gene] == '='):
            m=1 #nothing
        
        else:
            new_child[mutated_gene] = symbols[random_symbol_number]
            new_child = co.list2string(new_child)
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
    Draw.MolToMPL(Chem.MolFromSmiles(child), size = (100, 100))
