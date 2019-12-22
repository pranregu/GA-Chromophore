# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 23:47:06 2019

@author: Pran Regu
"""

from ase.io import read
from rdkit import Chem
from ase.calculators.emt import EMT

import Multipoint as co
import scoring_functions as sc

def energymin(mol, args):
    m = Chem.MolFromSmiles(mol)
    m3 = Chem.AddHs(m)
    print(Chem.MolToMolBlock(m3), file = open('molfile.mol', 'w'))
    atoms = read('molfile.mol')
    atoms.set_calculator(EMT())
    score = atoms.get_potential_energy()
    return(score)

def calculate_scores(population,function,scoring_args):
  scores = []
  for gene in population:
    score = function(gene,scoring_args)
    scores.append(score)
    
  return scores

    
if __name__ == "__main__":
    co.average_size = 39.15
    co.size_stdev = 3.50
    Compound = 'C1=CC=CC=C1C=C'
    target = Chem.MolFromSmiles(Compound)
    co.string_type = 'smiles'
    string = 'C1=CC=CC=C1C2=CC=CC=C2'
    score = energymin(string,[target])
    print(score)
