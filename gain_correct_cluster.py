"""THESE ARE INSTRUCTIONS"""

#fix the defect file
sed '/4664 2/s//4663 4/' defects.txt >defects-4row.txt

file 1:20170617_gaincor.sh

#!/bin/bash

counter=0

while IFS='' read -r line || [[-n "$line" ]]; do

        echo "Counter is at $counter"

        echo "Starting $line ${line:0:14}_gaincorr.mrcs"

        /programs/x86_64-linux/imod/4.9.2/bin/clip norm -m 1 -D defects-4row.txt $line SuperRef_correctinfohere.dm4 ${line:0:14}_gaincorr.mrcs

        echo "Finished with $line"

        counter=$((counter+1))

done < "$1"



file 2: qsub_gaincor.sh

#!/bin/bash

#PBS -S /bin/bash

# Inherit all current environment variables

#PBS -V

# Job name

#PBS -N 0614_gaincor

# Keep Output and Error in real time

#PBS -k eo

# Queue name

#PBS -q br

# Specify the number of nodes and thread (ppn) for your job.

#PBS -l nodes=1:ppn=40

# Do not exceed this amount of memory per processor

##PBS -l pmem=4gb

# Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS

#PBS -l walltime=1000:00:00

#PBS -l cput=1000:00:00

###################################################################################################

# Source

#source /programs/sbgrid.shrc

# Switch to the working directory;

#cd $PBS_O_WORKDIR

./20170617_gaincor.sh input_tif.txt



make an input file

ls *.tif >input_tif.txt



submit to cluster

qsub -k oe qsub_gaincor.sh
