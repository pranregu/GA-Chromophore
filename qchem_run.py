import os

os.system('python3 gen.py')
#execfile('gen.py')
 
#qchem.in to qchem.run
os.system('python3 makerun.py')
#execfile('makerun.py')s
  
#Submit jobs
os.system('python3 submit_runs.py')
#execfile('submit_runs.py')