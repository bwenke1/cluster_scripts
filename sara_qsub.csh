#!/bin/tcsh      
#PBS -S /bin/tcsh
# Inherit all current environment variables
#PBS -V
# Job name
#PBS -N XXXnameXXX
# Keep output and error
#PBS -k eo
# Queue name
#PBS -q XXXqueueXXX
# Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=XXXnodesXXX:ppn=XXXdedicatedXXX:gpus=8
# Do not exceed this amount of memory per processor
##PBS -l pmem=4gb
# Tell PBS the anticipated run time for your job, where walltime=HH:MM:SS
#PBS -l walltime=1000:00:00
#PBS -l cput=1000:00:00

#######################################################################################
# Source the SBGrid script
source /programs/sbgrid.cshrc
which mpiexec

# Switch to working directory
cd $PBS_O_WORKDIR

# Run
mpirun -n XXXmpinodesXXX XXXcommandXXX
