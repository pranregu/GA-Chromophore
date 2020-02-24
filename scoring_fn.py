# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 15:29:18 2020

@author: Pran Regu
"""

from ase.io import read
from ase.calculators.emt import EMT

from pysmiles import read_smiles

import Multipoint as co

def energymin(mol):
    m = read_smiles(mol, explicit_hydrogen=True)
    #print(Chem.MolToMolBlock(m3), file = open('molfile.mol', 'w'))
    print(m, file = open('molfile.mol', 'w'))
    atoms = read('molfile.mol')
    atoms.set_calculator(EMT())
    score = atoms.get_potential_energy()
    return(score)

def calculate_scores(population,function):
  scores = []
  for gene in population:
    score = function(gene)
    scores.append(score)
    
  return scores

    
if __name__ == "__main__":
    co.average_size = 39.15
    co.size_stdev = 3.50
    string = 'C1=CC=CC=C1C2=CC=CC=C2'
    score = energymin(string)
    print(score)

