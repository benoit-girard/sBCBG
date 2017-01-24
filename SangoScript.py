#!/usr/bin/python
# -*- coding: utf-8 -*-

import shlex
import subprocess

print "*** /!\ For simulations on Sango, check that interactive=False in LGNeurons.py ***"
print "*** and that the code to be executed has been made executable: chmod +x toto.py ***"

header = '#!/bin/bash \n\n'

IDstring = 'popo'

slurmOptions = ['#SBATCH --time=00:10:00 \n',
                '#SBATCH --partition=compute \n',
                '#SBATCH --mem-per-cpu=1G \n',
                '#SBATCH --ntasks=1 \n',
                '#SBATCH --cpus-per-task=2 \n',
                '#SBATCH --job-name=sBCBG_'+IDstring+'\n',
                '#SBATCH --input=none\n',
                '#SBATCH --output="runLogs/'+IDstring+'.out" \n',
                '#SBATCH --error="runLogs/'+IDstring+'.err" \n',
                '#SBATCH --mail-user=benoit.girard@isir.upmc.fr \n',
                '#SBATCH --mail-type=BEGIN,END,FAIL \n',
                ]

moduleLoad = ['module load nest/2.10 \n']

# build a parameter string strParams

# test param string:
strParams = '2644. 53. 8. 25. 14. 3000. 100. 9. 4.37 1.3 1.38 1.3 1. 13. 11. 100. 1. 1. 30. 70. 50. 1. 2. 25. 9. 15. 25. 9. 25. 9. 8. 2644. 25. 2644. 8. 23. 9.'

# number of neurons:
#strParams = '2644. 53. 8. 25. 14. 3000. 100. 9. '

# synaptic weight gains:
# strParams += ''

# write the script file accordingly
script = open(IDstring+'.slurm','w')
script.writelines(header)
script.writelines(slurmOptions)
script.writelines(moduleLoad)
script.writelines('./testFullBG.py '+strParams)
script.close()

print './testFullBG.py '+strParams

# execute the script file
p= []
command = 'sbatch '+IDstring+'.slurm'
p.append(subprocess.Popen(shlex.split(command)))

'''
command = 'srun '+options+'"python testFullBG.py'+params+'" &'
print command
formattedComm = shlex.split(command)

print "Launching BG simulation"
p.append(subprocess.Popen(formattedComm))
'''
