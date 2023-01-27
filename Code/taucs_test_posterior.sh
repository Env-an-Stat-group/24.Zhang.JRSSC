#!/bin/csh
#$ -M jzhang19@nd.edu
#$ -m ae
#$ -pe smp 1
#$ -q long
#$ -N degree1_post_NS
#$ -t 1-21

module purge
module load R/4.2.0/gcc/8.5.0 gdal geos udunits
Rscript taucs_test_posterior.R
