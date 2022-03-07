#!/bin/bash

echo $1
qsub -F "$1" submit.qsub
