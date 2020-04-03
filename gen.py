# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 15:50:44 2020

@author: Pran Regu
"""

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.inputs import QCInput
import glob

optRem ={
    'basis':'6-311G**',
    'job_type':'opt',
    'exchange':'B3LYP',
    'dft_d':'d3_bj',
    'scf_convergence':'8',
    'sym_ignore':'true',
    'geom_opt_max_cycles':'600',
} 

Energy=[['BASIS  =  6-311G**'],
            ['JOB_TYPE  =  SP'],
            ['MEM_STATIC  =  4000'],
            ['MEM_TOTAL  =  5000'],
            ['METHOD  =  B3LYP'],
            ['SCF_CONVERGENCE  =  8'],
            ['SCF_MAX_CYCLES  =  600'],
            ['SYMMETRY_IGNORE  =  1'],
            ['SYMMETRY_INTEGRAL  =  0'],
            ['UNRESTRICTED = 1']]

i = 0

for xyz in glob.glob("*.xyz"):
    #build molecule optimization files    
    #i = int(xyz.split('.')[0])
    i+=1    
    i_string = str(i)
    molM = Molecule.from_file(xyz)
    molM.set_charge_and_spin(charge = 0,spin_multiplicity = 1)
    NewInput = QCInput(molecule = molM, rem = optRem) 
    file = i_string +'.inp'
    NewInput.write_file(file)
    #print(file)
    with open(file, 'a') as outfile:
        outfile.write('\n\n')
        outfile.write('@@@\n')
        outfile.write('$molecule\n')
        outfile.write('  read\n')
        outfile.write('$end\n\n')
        outfile.write('$rem\n')
        for j in range(0,len(Energy)):
            outfile.write('  ')
            outfile.write(Energy[j][0])
            outfile.write('\n')
        outfile.write('$end\n\n')
    #print(file)
