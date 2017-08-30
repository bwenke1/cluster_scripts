#!/bin/tcsh
#setenv MALLOC_CHECK_ 0
#PBS -e $PBS_O_WORKDIR/XXXerrfileXXX
#PBS -o $PBS_O_WORKDIR/XXXoutfileXXX
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
### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=200:00:00
#################################
### Switch to the working directory;
# Environment
source /programs/labcshrc
@ np = XXXmpinodesXXX * XXXthreadsXXX
cd $PBS_O_WORKDIR
#echo $PBS_NODEFILE
set nodes = `cat $PBS_NODEFILE | sort -u`
#echo $nodes
set i = 1
set hostlist = $nodes[$i]
@ i++
while ( $i <= $#nodes )
set hostlist = `echo $hostlist,$nodes[$i]`
@ i++
end

#-## Run:
mpirun --prefix /programs/x86_64-linux/system/openmpi/1.8.4 -bynode -H $hostlist -np $np XXXcommandXXX
#mpirun --prefix /programs/x86_64-linux/system/openmpi/1.8.4 --map-by node $hostlist -np $np XXXcommandXXX
echo "done"
