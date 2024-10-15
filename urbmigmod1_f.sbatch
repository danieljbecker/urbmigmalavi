#!/bin/bash
#
#SBATCH --partition=longjobs
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=1024
#SBATCH --output=urbmigmod1_f_output_%J.txt
#SBATCH --error=urbmigmod1_f_error_%J.txt
#SBATCH --time=144:00:00
#SBATCH --job-name=urbmigmod1_f
#SBATCH --mail-user=tbverrett@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/verrett/R
#
#################################################
module load R/4.2.1-foss-2022a
Rscript OSCER_urbmigmalavi_05 pcms_noLeuco.R