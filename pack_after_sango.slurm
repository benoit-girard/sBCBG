#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --time=04:00:00 
#SBATCH --partition=compute 
#SBATCH --mem-per-cpu=500M 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --job-name=pack_$2_$1
#SBATCH --input=none
#SBATCH --mail-user=
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --dependency=afterany:$1

if [ -d "$2" ]; then 
 echo "gathering params_score.csv files..." && time find $2 -name 'params_score.csv' > $2.files && echo "done." && echo "gathering base experiment files..." && time find $2 -maxdepth 1 -type f >> $2.files && echo "done." && echo "packing all files into a single archive..." && time cpio -oa -O $2.cpio < $2.files && echo "done." && echo "removing original experiment files..." && time find $2 -maxdepth 1 -name '???' -exec rm -r {} +
else
  echo "could not locate $2"
fi 

echo "done."
EOT
