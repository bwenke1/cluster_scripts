###################################################################################################
#!/bin/tcsh
#PBS -S /bin/tcsh
# Inherit all current environment variables
#PBS -V
# Job name
#PBS -N XXXnameXXX
# Keep Output and Error
#PBS -k eo
# Queue name
#PBS -q XXXqueueXXX
# Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=XXXnodesXXX:ppn=XXXdedicatedXXX
# Do not exceed this amount of memory per processor
#PBS -l pmem=3gb
# Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=10000:00:00
#PBS -l cput=10000:00:00
###################################################################################################
# Source the SBGrid script
source /programs/sbgrid.shrc
# Switch to the working directory;
cd $PBS_O_WORKDIR
# Run
mpirun --mca btl_tcp_if_include 172.20.100.0/24 XXXcommandXXX
