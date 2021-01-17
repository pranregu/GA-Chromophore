# Crossover function

from rdkit import Chem
from rdkit.Chem import AllChem

import random
import numpy as np

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

# Checks if the molecule is legitimate:
def string_OK(string):
  mol = string2mol(string)
  if not mol:
      return False
  try:
    Chem.SanitizeMol(mol)
    test_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    if test_mol == None:
      return None
    target_size = size_stdev*np.random.randn() + average_size #parameters set in GA_mol
    if mol.GetNumAtoms() > 5 and mol.GetNumAtoms() < target_size:
      return True
    else:
      return False
  except:
    return False

def list2string(list):
    string = ''.join(list)
    return string

def string2mol(string):
    try:
        mol = Chem.MolFromSmiles(string)
        return mol
    except:
        return None

# Selects a random splicing point:
def cut_point(parent):
    m = random.randint(0, len(parent) - 1)
    return m

# Multipoint crossover operator:
def crossover(parent_a,parent_b):
    
    for _ in range(50):
       
        str1 = ''
        str2 = ''
        
        #Splicing parent a
        while len(str1)<13:
            index_a = cut_point(parent_a)
            a = len(parent_a)
            count = 0
            i = index_a
            # Splicing parent_a
            if parent_a[index_a] != '(':
                for i in range(index_a, a):
                    if parent_a[i] == '(':
                        #index_a=1
                        break
            
            index_a = i 
        
            for i in range(index_a, a):
                if parent_a[i] == '(':
                    count = count+1
                elif parent_a[i] == ')': 
                    count = count-1
                if parent_a[i] == ')' and count == 0:  
                    break
                i = i+1
                
            str1 = parent_a[index_a : i+1] 
        
        #Splicing parent b        
        while len(str2)<13:
            index_b = cut_point(parent_b)
            b = len(parent_b)
            count = 0
            j = index_b
            # Splicing parent_b
            if parent_b[index_b] != '(':
                for j in range(index_b, b):
                    if parent_b[j] == '(':
                        #index_b=1
                        break
            
            index_b = j 
        
            for j in range(index_b, b):
                if parent_b[j] == '(':
                    count = count+1
                elif parent_b[j] == ')': 
                    count = count-1
                if parent_b[j] == ')' and count == 0:  
                    break
                j = j+1
                
            str2 = parent_b[index_b : j+1] 
 
        # Crossover - replacing str from parent a with str from b 
        child_string = parent_a.replace(str1, str2, 1)
       
        # Check if child_string is a valid molecule
        if string_OK(child_string):
            return (child_string)
   
    return None    
 

# Test case
if __name__ == "__main__":
    average_size = 39.15
    size_stdev = 3.50
    parent_a = 'C1=CC=C(C2=CC=CC=C2C2=CC(C4=C(C5=CC=CC=C5)C=CC=C4)=C(C3=CC=CC=C3)C=C2)C=C1'
    parent_b = 'C1=CC=C(C2=CC=CC=C2C2=CC=CC=C2C2=CC=CC=C2)C=C(C(=O)O)1'
    child = crossover(parent_a,parent_b)
    print(child)
