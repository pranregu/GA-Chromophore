# Phenyl_GA
This code is a genetic algorithm to generate terphenyl derivatives which uses a screening function based on LUMO values of the molecules.

Prerequisites:

QChem
Python packages: a. ase b. rdkit c. pymatgen
Installing: Qchem can be downloaded from http://www.q-chem.com/. IMPORTANT: Qchem requires a license. /n
All of the above mentioned packages have conda installtions.

Modules:

Multipoint: Crossover operator Takes 2 molecules as SMILES and produces a new molecule with features of both the parent molecules. /n
string_mutate: Mutation operator Takes 1 molecule as SMILES and replaces one atom/functional group with another. /n
xyz_gen: Takes a SMILES representation and creates a xyz file /n
gen: Creates a Qchem input file /n
qchem_run: Makes slurm files. (file.run) /n
submit_runs: Submits slurm files. /n
New_scoring_fn: Assigns scores to molecules based on their LUMO energies. /n
string_GA: Generates a population of molecules that evolve in each generation(iteration) with increasing values of LUMO energies. /n
terphenyl.smi: File containing terphenyl SMILES, which are used to make the initial population in string_GA. /n
/n
Running: Run the file string_GA.
