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
#PBS -l nodes=XXXmpinodesXXX:ppn=XXXthreadsXXX
#setenv OMP_NUM_THREADS 1

### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=10000:00:00
#PBS -l cput=10000:00:00

#################################
### Switch to the working directory;
# Environment
source /programs/labcshrc
cd $PBS_O_WORKDIR

#mpirun --prefix /programs/x86_64-linux/system/openmpi/1.8.4 --mca btl_tcp_if_include 172.20.100.0/24 XXXcommandXXX
mpirun --mca btl_tcp_if_include 172.20.100.0/24 XXXcommandXXX

#echo "done"
#echo mpirun --prefix /programs/x86_64-linux/system/openmpi/1.8.4 --npernode XXXthreadsXXX --hostfile $PBS_NODEFILE XXXcommandXXX

