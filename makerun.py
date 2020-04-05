import glob

for inp in glob.glob('*.inp'):
    name = inp[0:-4]
    f = open(name + ".run", "w")
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH --ntasks-per-node=1\n')
    f.write('#SBATCH --mem-per-cpu=2GB\n')
    f.write('#SBATCH --time=01:00:00\n')
    f.write('#SBATCH --output=' + name + '.out\n')
    f.write('#SBATCH --partition=sharada\n\n')

    f.write('export QC=/usr/usc/qchem/default\n')
    f.write('export QCAUX=$QC/qcaux\n')
    f.write('export QCPLATFORM=LINUX_Ix86_64\n')
    f.write('export QCRSH=ssh\n')
    f.write('export PATH=$QC/bin:$PATH\n')
    f.write('export QCSCRATCH=$TMPDIR\n\n')

    f.write('ulimit -c 0\n')
    f.write('cd ${SLURM_SUBMIT_DIR}\n')
    f.write('hpc-login2\n')
    f.write('source /usr/usc/qchem/default/qcenv.sh\n')
    f.write('qchem -nt 1 ' + inp + '\n')
    f.write('cp -R "$TMPDIR" "$SLURM_SUBMIT_DIR" ')
    f.close()
