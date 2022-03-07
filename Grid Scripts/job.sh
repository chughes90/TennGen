#!/bin/bash 

##create calculation sub-directories######################################################################
#"To have completely unique seeds, need to multiply by a number which is greater than the number of jobs."
#"Protecting in case multiple jobs start at once."
export outdir=${PBS_O_WORKDIR}/$PBS_JOBID
#export outdir=/lustre/haven/user/chughe26/scratch/jetbin6/$PBS_JOBID
mkdir $outdir
echo $PBS_TASKNUM
JOBID=`echo $PBS_JOBID | cut -d . -f 1`
NEWTASKNUM=$((PBS_TASKNUM * 20000))
JOB_NUMBER=$((NEWTASKNUM + $JOBID))
echo $JOB_NUMBER 
pwd
mkdir $outdir/$JOB_NUMBER

##copy input codes#############################################################################
cp /nics/b/home/chughe26/TennGen/BkgrLoad.h $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/TennGen/BkgrLoad.cxx $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/TennGen/maindriver_PreGen.C $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/TennGen/TennGen_PreGen_Macro.C $outdir/$JOB_NUMBER/ 

##copy input files (HI ONLY !!!)###############################################################
echo $1
echo $PBS_TASKNUM
cp /lustre/haven/user/chughe26/scratch/BKGD_EVENTS/0_5_percent/v1_v5/FINAL_FILES/$1/$PBS_TASKNUM/*.root $outdir/$JOB_NUMBER/

##load modules#################################################################################
## Do NOT change the order of the following files
cd $outdir/$JOB_NUMBER
export PATH=$PATH:/lustre/haven/proj/UTK0019/python-virtual-environments/env/bin
export ALIBUILD_WORK_DIR=/lustre/haven/proj/UTK0019/aliroot5/alice/sw
eval "`alienv shell-helper`"
#alienv load AliPhysics/latest
alienv --no-refresh load AliPhysics/latest


##initiate calculation#########################################################################
root -b -l -q  "maindriver_PreGen.C(100000000, $JOB_NUMBER , 0, 0, 1, kFALSE , kTRUE)"


# num_bkgd events / seed / Harmonic Flag centrality bin / Runtime Stats Flag / Run on Grid Flag

################################################################################################
#rm -r $outdir/$JOB_NUMBER


################################################################################################
#rm -r $outdir/$JOB_NUMBER

