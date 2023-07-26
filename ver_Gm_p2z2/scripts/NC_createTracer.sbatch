#!/bin/bash -x
#SBATCH -J NCproc
#SBATCH -A ACF-UTK0105
#SBATCH --exclude=ilp1119,ilp1120
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --time=48:00:00
#SBATCH --output=/nfs/home/ecarr/logfiles/NCproc.%j.out
#SBATCH --error=/nfs/home/ecarr/logfiles/NCproc.%j.err
#SBATCH --qos=condo
#SBATCH --partition=condo-dtalmy

module purge
module load zlib/1.2.11-gcc
#load anaconda
module load anaconda3/
source $ANACONDA3_SH
conda activate ma38

RUN_ID=$SLURM_JOB_ID

#setup dir info
if [ -z "$RUN_DIR" ] 
then
  echo "Please set \$RUN_DIR"
  echo "Input failed"
  exit 1
else
	echo "\$RUN_DIR   : $RUN_DIR"
fi
#RUN_DIR=/lustre/isaac/proj/UTK0105/ProjAnal/T41_mort/Model_runs/run_t41mort_5389956.apollo-acf
RESpath=/NC_trace
RESanal=/NC_analysis
RESscripts=/scripts
echo Result processing 

#echo RUN_DIR   : $RUN_DIR
#echo RESpath   : $RESpath
echo mod. Result Dir   : $RUN_DIR$RESpath
echo Analysis Dir      : $RUN_DIR$RESanal
echo script Dir        : $RUN_DIR$RESscripts
echo Model Run ID      : $MODEL_RUNID
#setup Output processing
alias python=python3
cd $RUN_DIR
mkdir $RUN_DIR$RESpath
mkdir $RUN_DIR$RESanal
#
# gud files
process_grid=true
if [ "$process_grid" = true ]; then
    echo Processing Grid
    $RUN_DIR$RESscripts/gluemncbig -o grid.nc  $RUN_DIR/ecco_gud*/grid*.nc
    nccopy -k netCDF-4 -d 2 -s  grid.nc $RUN_DIR$RESpath/grid.v4c.nc
    rm grid.nc
fi
#
# gud files
# this has not been updated
process_gud=true
if [ "$process_gud" = true ]; then
    echo Processing GUD files
    ModelFilePre=3d
    #ModelFilePre=gudnp2z2d
    ModelFilePreOut=3d
    declare -a gudArray=("TRAC01" "TRAC02" "TRAC03" "TRAC04" "TRAC05" \
                        "TRAC06" "TRAC07" "TRAC08" "TRAC09" "TRAC10" \
                        "TRAC11" "TRAC12" "TRAC13" "TRAC14" "TRAC15" \
                        "TRAC16" "TRAC17" "TRAC18" "TRAC19" "TRAC20" \
                        "TRAC21" "TRAC22" "TRAC23" "TRAC24" "TRAC25" \
                        "TRAC26" "TRAC27" "TRAC28" "TRAC29" "TRAC30" \
                        "TRAC31" "TRAC32" "TRAC33" "TRAC34" "TRAC35" \
                       # "TRAC36" "TRAC37" "TRAC38" "TRAC39" "TRAC40" \
                       # "TRAC41" "TRAC42" "TRAC43" "TRAC44" "TRAC45" \
                       # "TRAC46" "TRAC47" "TRAC48" "TRAC49" "TRAC50" \
                       # "TRAC51" "TRAC52" "TRAC53" "TRAC54" "TRAC55" \
                       # "TRAC56" "TRAC57" "TRAC58" "TRAC59" "TRAC60" \
                        "PAR01" "PP" "BP")
    
    
    
    
    for Tracer_Name in ${gudArray[@]}; do
        
        echo $Tracer_Name
        $RUN_DIR$RESscripts/gluemncbig -o $RUN_DIR/$ModelFilePre.$Tracer_Name.nc -v $Tracer_Name $RUN_DIR/mnc_test*/$ModelFilePre.0000000*.nc
        nccopy -k netCDF-4 -d 2 -s  $ModelFilePre.$Tracer_Name.nc $RUN_DIR$RESpath/$ModelFilePreOut.$Tracer_Name.v4c.nc
        rm $ModelFilePre.$Tracer_Name.nc
        # add calendar 360day to T for appropriate cftime interpretation
        ncatted -a calendar,T,c,c,'360_day' $RUN_DIR$RESpath/$ModelFilePreOut.$Tracer_Name.v4c.nc
        # I put Rename proj usr bin
        python3 $RUN_DIR$RESscripts/MU_renameNCtrace.py -fn $ModelFilePreOut.$Tracer_Name.v4c.nc -D $RUN_DIR$RESpath
	rm $ModelFilePreOut.$Tracer_Name.v4c.nc 
    done
else
	echo "\$process_gud not run : $process_gud"    
fi
#
# ptracers files
process_ptracers=false
if [ "$process_ptracers" = true ]; then
    echo Processing ptracers files
    ModelFilePre=ptracers
    declare -a ptArray=("DIC" "NH4" "NO2" "NO3" "PO4"\
                        "SiO2" "FeT" "DOC" "DON" "DOP"\
                        "DOFe" "POC" "PON" "POP" "POSi"\
                        "POFe" "PIC" "ALK" "O2" "c01"\
                        "c02" "c03" "c04" "c05" "c06"\
                        "c07" "c08" "c09" "c10" "c11"\
                        "c12" "c13" "c14" "c15" "c16"\
                        "c17" "c18" "c19" "c20" "c21"\
                        "c22" "c23" "c24" "c25" "c26"\
                        "c27" "c28" "c29" "c30" "c31"\
                        "c32" "c33" "c34" "c35" "c36"\
                        "c37" "c38" "c39" "c40" "c41")
    for Tracer_Name in ${ptArray[@]}; do
        
        echo $Tracer_Name
        $RUN_DIR$RESscripts/gluemncbig -o $RUN_DIR/$ModelFilePre.$Tracer_Name.nc -v $Tracer_Name $RUN_DIR/ecco_gud*/$ModelFilePre.0000000*.nc
        nccopy -k netCDF-4 -d 2 -s  $ModelFilePre.$Tracer_Name.nc $RUN_DIR$RESpath/$ModelFilePre.$Tracer_Name.v4c.nc
        rm $ModelFilePre.$Tracer_Name.nc
		# add calendar 360day to T for appropriate cftime interpretation
        ncatted -a calendar,T,c,c,'360_day' $RUN_DIR$RESpath/$ModelFilePre.$Tracer_Name.v4c.nc
        # I put Rename proj in script dir
        #python3 $RUN_DIR$RESscripts/MU_renameNCtrace.py -fn $ModelFilePre.$Tracer_Name.v4c.nc -D $RUN_DIR$RESpath
    done
else
	echo "\$process_ptracers not run : $process_ptracers"    
fi
#
# remember 360 day calendar
#
# output formated env variables
echo 'RUN_DIR='"$RUN_DIR"
echo 'MODEL_RUNID='"$MODEL_RUNID"
echo 'RESpath='"$RESpath"
echo 'RESanal='"$RESanal"
echo 'RESscripts='"$RESscripts"
####
# file to source in shell to allow copy and past of most of these commands
printf '%s\n' 'RUN_DIR='"$RUN_DIR" >>  $RUN_DIR$RESscripts'/envvar.txt'
printf '%s\n' 'MODEL_RUNID='"$MODEL_RUNID" >> $RUN_DIR$RESscripts'/envvar.txt'
printf '%s\n' 'RESpath='"$RESpath" >> $RUN_DIR$RESscripts'/envvar.txt'
printf '%s\n' 'RESanal='"$RESanal" >> $RUN_DIR$RESscripts'/envvar.txt'
printf '%s\n' 'RESscripts='"$RESscripts" >> $RUN_DIR$RESscripts'/envvar.txt'
#
cd $RUN_DIR
######
##  simple line graphs of Tracer_Name
#echo run line graph
#python3 $RUN_DIR$RESscripts/DA_multitracesumplot.py \
#	-NCD $RUN_DIR$RESpath \
#	-T $MODEL_RUNID \
#	-D $RUN_DIR$RESanal \
#	-L [c01,c02,c03]
#######################################
## Output based on the 3 last years
##
#cd $RUN_DIR
#echo run DA_pt_p2z1_3yr
#python3 $RUN_DIR$RESscripts/DA_pt_p2z1_3yr.py \
#	-NCD $RUN_DIR$RESpath \
#	-T $MODEL_RUNID \
#	-D $RUN_DIR$RESanal
###########
## Graph like Ward et al
##
#echo run DA_stackCarbon
#python3 $RUN_DIR$RESscripts/DA_stackCarbon.py \
#	-NCD $RUN_DIR$RESpath \
#	-T $MODEL_RUNID \
#	-D $RUN_DIR$RESanal
#############
## globe map and hist of slopes
## each cell over a 3yr summed depth
#python3 $RUN_DIR$RESscripts/DA_RegGlbSlope.py \
#	-NCD $RUN_DIR$RESpath \
#	-T $MODEL_RUNID \
#	-D $RUN_DIR$RESanal
##
############################################
##
## Single year stuff using GRPD_yr
##
#echo run MU_ag_grp_yr.py
#python3 $RUN_DIR$RESscripts/MU_ag_grp_yr.py \
#	-MD $RUN_DIR
##
#echo Ext_ZP_Data
#python3 $RUN_DIR$RESscripts/MU_extpt_data.py \
#	-NC $RUN_DIR$RESpath/GRPD_yr.nc  \
#	-L $RUN_DIR/scripts/Ir_Hat_HACKcombo_unix.csv \
#	-T $MODEL_RUNID \
#	-PFO $RUN_DIR$RESpath/$MODEL_RUNID\_Ext_ZP_Data.csv 
## 
#echo Gen ZP multi cruise slope analysis
#python3 $RUN_DIR$RESscripts/DA_ZPmulticruise.py \
#	-DF $RUN_DIR$RESpath/$MODEL_RUNID\_Ext_ZP_Data.csv   \
#	-T $MODEL_RUNID \
#	-PFO $RUN_DIR$RESpath/PtVB_$MODEL_RUNID\_Ext_VB_Data.csv 
#
#
#
#	
#
echo process finished
