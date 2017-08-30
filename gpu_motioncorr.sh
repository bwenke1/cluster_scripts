#!/bin/tcsh      
#PBS -S /bin/tcsh
# Inherit all current environment variables
#PBS -V
# Job name
#PBS -N MotionCorr/job065
# Keep output and error
#PBS -k eo
# Queue name
#PBS -q gpu
# Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=1:ppn=4:gpus=8
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
/programs/x86_64-linux/system/openmpi/1.8.4/bin/mpirun -n 4 'which relion_run_motioncorr_mpi' --i ~/Desktop/data/polara/20170522_BENH/Import/job056/movies.star --o ~/Desktop/data/polara/20170522_BENH/MotionCorr/job065/ --save_movies  --first_frame_sum 1 --last_frame_sum 30 --bin_factor 2 --motioncorr_exe /programs/x86_64-linux/motioncor2/20170130/MotionCor2-01-30-2017 --bfactor 150 --use_motioncor2 --angpix 0.7 --patch_x 10 --patch_y 10 --gpu "0:1:2:3" --dose_weighting --voltage 300 --dose_per_frame 1.65 --preexposure 0
