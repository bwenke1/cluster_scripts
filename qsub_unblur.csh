#!/bin/tcsh
#setenv MALLOC_CHECK_ 0


#PBS -S /bin/tcsh
### Inherit all current environment variables
#PBS -V
### Job name
#PBS -N XXXnameXXX
### Keep Output and Error
#PBS -k eo
### Queue name
#PBS -q XXXqueueXXX

# Send email to user
#PBS -M bwenke@caltech.edu

### Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=1:ppn=20

### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=10000:00:00
#################################
### Switch to the working directory;
# Environment
source /programs/labcshrc
cd $PBS_O_WORKDIR

mpirun -n 5 --prefix /programs/x86_64-linux/system/openmpi/1.8.4 --hostfile $PBS_NODEFILE unblur --j 4 --memory_per_thread 4 

echo "done"
echo mpirun --prefix /programs/x86_64-linux/system/openmpi/1.8.4 --hostfile $PBS_NODEFILE XXXcommandXXX

