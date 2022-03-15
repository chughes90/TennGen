#!/bin/bash
#for (( i=1949417; i<=1959470; i++ ))
#for (( i=1896303; i<=1896304; i++ ))
#change to your directory
for (( i=3647039; i<=3649556; i++ )) 
do
   JOB_NUM=$((i))
   DIR_STR=$"/lustre/haven/user/tmengel/TennGen/running/TennGen/output/$JOB_NUM"".apollo-acf"
   echo $DIR_STR
   FILE_STR=$"$JOB_NUM"".root"
if [ -d $DIR_STR ]
then
hadd $FILE_STR $DIR_STR/*/TennGen_PreGen*.root	
fi
done
hadd merged_final.root /lustre/haven/user/tmengel/TennGen/running/TennGen/output/3*.root
