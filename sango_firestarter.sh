#!/bin/bash

# This script initialize the run directory according to:
# - the name of parent directories, which together form the task number of the job
# - the possible values over which each parameter is screened, with one value per line in files named 'XYZ.txt' (where XYZ is the parameter name; for example, a file named 'IeMSN.txt' and containing 3 lines define 3 different values for IeMSN)
# The last part of the script is to be filled in by run.py, and should contain the copy commands as well as the run command (such as 'cp ../LGneurons.py /' and 'python testFullBG.py')

#########
# 1. retrieve xp number from directories
#########

dir=$(pwd)
xpnumber="${dir: -11:3}${dir: -7:3}${dir: -3}"
xpnumber=`expr $xpnumber + 0`
echo "XPNUMBER IS: $xpnumber"

xpbase="${dir: 0:-12}"
#echo "XPBASE IS: $xpbase"

xpname="${xpbase##*/}"
#echo "XPNAME IS: $xpname"

##########
# 2. retrieve number of parameters and number of variations somehow
##########

declare -A parameters
for p in `ls -1 $xpbase/*.txt | sort`
do
  basep=${p: 0:-4}
  namep="${basep##*/}"
  linep=`cat $p | wc -l`
  #echo "$namep has $linep params"
  parameters["$namep"]="$linep"
done

##########
# 3. compute offset of individual parameters
##########

declare -A offset
remxp=$xpnumber
for p in "${!parameters[@]}"
do
  offset["$p"]=$(($remxp%${parameters[$p]}))
  #echo " offset for $p is ${offset[$p]}"
  remxp=$(($remxp-${offset[$p]}))
  #echo " remainder is $remxp"
  remxp=$(($remxp/${parameters[$p]}))
  #echo " next value: $remxp"
done

##########
# 4. build the adapted modelParams.py using the skeleton and the offsets
##########

# Work in a local temporary directory
tempdir=`mktemp -d --tmpdir jean-lienard.XXXXXXXXXXXXXXXXXX`
cd $tempdir

cp $xpbase/baseModelParams.py $(pwd)/modelParams.py
echo "" >> $(pwd)/modelParams.py

for p in "${!parameters[@]}"
do
  pval=`head -n $((1+${offset["$p"]})) $xpbase/$p.txt | tail -n 1`
  #echo "selected val for $p: $pval"
  echo "params.update({'$p': $pval})" >> $(pwd)/modelParams.py
done

##########
# 5. copy other files from main directory and start the script
##########

mkdir $(pwd)/log # the log directory

#### SHOULD BE AUTO-GENERATED FROM THIS POINT FORWARD BY `run.py`

