#!/bin/bash 
#SBATCH -J 8x8_darwin
#SBATCH -A ACF-UTK0105
#SBATCH --exclude=clr0812,clr0813,ilp1119,ilp1120
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/home/ecarr/logfiles/Gmp2z1.%j.out
#SBATCH --error=/nfs/home/ecarr/logfiles/Gmp2z1.%j.err
#SBATCH --qos=campus
#SBATCH --partition=campus
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo ---Modules 
#module load netcdf/4.4.1.1
#module load  gcc/6.3.0
#module load  hdf5/1.10.1
#module load  lapack/3.7.0
module purge
module load intel-compilers/2021.2.0
module load intel-mpi/2021.2.0
module load PE-intel
module list
NCDIR=/lustre/isaac/proj/UTK0105/usr
echo $NCDIR
export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
echo $LD_LIBRARY_PATH
shopt -s expand_aliases
alias python='/usr/bin/python2'
python2 --version
python  --version
VALpath=/lustre/isaac/proj/UTK0105/Darwin/GITLAB/gud-gud/verification/GUD_closure_DD
USERpath=/lustre/isaac/scratch/ecarr/runs
RUN_ID=$SLURM_JOB_ID
TAG_ID='Gm_p2z2'
DIR_NAME=run_$TAG_ID\_$RUN_ID
USER_DIR=$USERpath/$DIR_NAME
RUN_DIR=$VALpath/$DIR_NAME
RESULTS_DIR=mm_$RUN_ID
GUDBpath=/lustre/isaac/proj/UTK0105/Darwin/GITLAB/gud-gud
#
#Setup run 
echo RUN_ID    : $RUN_ID
echo TAG_ID    : $TAG_ID
echo VALpath   : $VALpath
echo USERpath  : $USERpath
echo RUN_DIR   : $RUN_DIR
echo USER_DIR   : $USER_DIR
cd $USERpath
if [ -d "$USER_DIR" ]; then
    # Will enter here if $DIRECTORY exists, even if it contains spaces
    echo "User run dir exists  failed"
    exit 1
fi

mkdir  $USER_DIR
lfs setstripe $USER_DIR -S 32m -i -1 -c 1
cd $VALpath 
if [ -d "$RUN_DIR" ]; then
    # Will enter here if $DIRECTORY exists, even if it contains spaces
    echo "RUN area ln dir exists  failed"
    exit 1
fi
cd $VALpath 
ln -s $USER_DIR ./$DIR_NAME
cd   $RUN_DIR
#ln -s $USER_DIR $RUN_DIR

#build

cp -r $VALpath/ver_$TAG_ID/* $RUN_DIR
# chg mortality modifier
# relative links to wokr you must be in the runarea symbolic link
cd   $RUN_DIR
pwd
cd ./code
# copy correct input for cluster
#cd $RUN_DIR
#cd ./input
#cp ./isaac/* .
#
cd $RUN_DIR
#
mkdir ./build
cd ./build

if 
    $GUDBpath/tools/genmake2 \
      -rootdir $GUDBpath \
      -mods ../code \
      -mpi  \
      -optfile $GUDBpath/tools/build_options/linux_amd64_ifort+impi_stampede2_skx_isaac
then
    echo "genmake2 succeeded"
else
    echo "genmake2 failed"
    exit 1
fi

if 
  make depend
then
    echo "make depend succeeded"
else
    echo "make depend failed"
    exit 1
fi

if 
  make 
then
    echo "make  succeeded"
else
    echo "make  failed"
    exit 1
fi
#
cd ..
ln -s ./input/* .
cp ./build/mitgcmuv .
#rm -rf ./build  &
#

#run model
#mpirun -np 16 ./mitgcmuv


#run model
unset I_MPI_PMI_LIBRARY 
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0   # the option -ppn only works if you set this before
cd $RUN_DIR
pwd
if
    /usr/bin/time --output=outtime_$RUN_ID.log -p sh -c ' /sw/isaac/compilers/intel/oneAPI_2021.2.0/mpi/latest/bin/mpirun -np 64  ./mitgcmuv 2>&1 | tee output.log'
    
then
    echo "RUN  succeeded"
else
    echo "RUN  failed"
    exit 1
fi

#Submit Output processing
#qsub ./scripts/NC_$TAG_ID.sh  -N NC$PBS_JOBID -v Pass_DIR=$RUN_DIR,Pass_ID=$RUN_ID,Pass_RN=$RESULTS_DIR
sbatch --job-name=NC_npzd --export=RUN_DIR=$USER_DIR,MODEL_RUNID=$RUN_ID  ./scripts/NC_createTracer.sbatch
