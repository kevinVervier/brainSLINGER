#!/bin/bash

#GSE directives
#$ -N Tis-SLINGER
#$ -S /bin/bash
#$ -cwd
#$ -q JM
#$ -pe smp 8

#or change the number of CPUs for large data...

#load R module
module load R

Rscript train_model.R $alpha $chunk

