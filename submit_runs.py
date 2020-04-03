import subprocess
import os
import glob

for run in glob.glob('*.run'):
    os.system('sbatch '+ run)
