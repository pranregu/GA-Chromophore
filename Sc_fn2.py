# Scoring function (Evaluation function):
from rdkit import Chem

from ase.calculators.qchem import QChem
from ase.optimize import LBFGS

from ase.build import molecule
from ase import Atoms
from ase.io import read,write

import xyz_gen as xg

# Opens generated Qchem output file and extracts LUMO:
def extract_LUMO(file): 
        number = 0   
        with open(file,"r") as output_file:
            for line in output_file:
                if '-- Virtual --' in line:
                    number+=1
        #print(number)            
        with open(file,"r") as output_file:
            LUMO=[]
            count = 0
            for line in output_file:
                if count == number:
                    LUMO.append(line.split()[0])
                    break
                if '-- Virtual --' in line:
                    count+=1
            LUMO = LUMO[0:6]
            for i in range(len(LUMO)):
                LUMO[i] = float(LUMO[i])
        return(LUMO)

# Function to calculate the scores of the molecules in the population:
def calculate_scores(population, generation):
  scores = []
  #smiles to xyz
  for i in range(len(population)):
      name = str(generation)+"mol"+str(i)     
      mol = Chem.MolFromSmiles(population[i])
      mol = Chem.AddHs(mol)
      xg.write_input_files(mol, name) 

      # QChem Optimization
      comp = read(name+'A.xyz')
      calc = QChem(jobtype='opt',
                   label=name,
                   method='PBE',
                   basis='3-21G',
                   scf_convergence='8',
                   scf_max_cycles='200')
      comp.set_calculator(calc)
      opt = LBFGS(comp)
      opt.run()

      # LUMO from qchem.out
      file = name+".out"
      score = extract_LUMO(file)
      scores.append(score)

  #Printing scores:
  f = open("Result.out", "a")
  print(generation, file=f)
  print(scores, file=f)   
  f.close()  

  return scores
