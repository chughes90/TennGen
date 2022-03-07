#!/bin/bash
for (( i=18118; i<=18144; i++ ))
do 
   echo $i 
   ./qsubmit.sh $i
done

#valid numbers 2474 - 2500 OR 7540 - 23704
