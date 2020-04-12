# Phenyl_GA

This code is a genetic algorithm to generate terphenyl derivatives which uses a screening function
based on LUMO values of the molecules.

Prerequisites:
1. QChem
2. Python packages:
a. ase
b. rdkit
c. pymatgen

Installing:
Qchem can be downloaded from http://www.q-chem.com/. IMPORTANT: Qchem requires a license.
All of the above mentioned packages have conda installtions. 

Modules:
1. Multipoint: Crossover operator
   Takes 2 molecules as SMILES and produces a new molecule with features of both the parent 
   molecules.
2. string_mutate: Mutation operator
   Takes 1 molecule as SMILES and replaces one atom/functional group with another.
3. xyz_gen:
   Takes a SMILES representation and creates a xyz file
4. gen:
   Creates a Qchem input file
5. qchem_run:
   Makes slurm files. (file.run)
6. submit_runs:
   Submits slurm files.
7. New_scoring_fn: 
   Assigns scores to molecules based on their LUMO energies.
8. string_GA:
   Generates a population of molecules that evolve in each generation(iteration) with increasing
   values of LUMO energies.
9. terphenyl.smi
   File containing terphenyl SMILES, which are used to make the initial population in
   string_GA.
   
Running:
Run the file string_GA.
   
