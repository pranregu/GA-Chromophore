# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 19:06:18 2020

@author: Pran Regu
"""

from ase.io import read
from rdkit import Chem

import xyz_gen as xg
import os
import time

def extract_LUMO(file):#change filename here
        with open(file,"r") as output_file:
            LUMO=[]
            count = 0
            for line in output_file:
                if count == 1:
                    LUMO.append(line.split()[0])
                    break
                if '-- Virtual --' in line:
                    count = 1
            LUMO = LUMO[0:6]
            for i in range(len(LUMO)):
                LUMO[i] = float(LUMO[i])
        return(LUMO)

def calculate_scores(population, args):
  scores = []
  #smiles to xyz
  for i in range(len(population)):
      name = "mol"+str(i)     
      #print(name)
      mol = Chem.MolFromSmiles(population[i])
      mol = Chem.AddHs(mol)
      xg.write_input_files(mol, name) 
      
  #xyz to qchem.in
  os.system('python3 gen.py')
 
  #qchem.in to qchem.run
  os.system('python3 makerun.py')
  
  #Submit jobs
  os.system('python3 submit_runs.py')

  #Wait
  time.sleep(720) 
  
  #LUMO from qchem.out
  for j in range(len(population)):  
      file = str(j+1)+".out"
      score = extract_LUMO(file)
      scores.append(score)
   
  return scores
