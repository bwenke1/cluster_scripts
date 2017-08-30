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
### Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=1:ppn=XXXmpinodesXXX
### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=200:00:00
#################################
### Switch to the working directory;
# Environment
source /programs/labcshrc
setenv NCPUS XXXmpinodesXXX
cd $PBS_O_WORKDIR

XXXcommandXXX

echo "done"
echo mpirun --prefix /programs/x86_64-linux/system/openmpi/1.8.4 --npernode XXXthreadsXXX --hostfile $PBS_NODEFILE XXXcommandXXX

