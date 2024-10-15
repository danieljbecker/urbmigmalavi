#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=urbmigmod1_f9_output_%J.txt
#SBATCH --error=urbmigmod1_f9_error_%J.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=urbmigmod1_f9
#SBATCH --mail-user=tbverrett@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/verrett/R
#
#################################################
module load R/4.2.1-foss-2022a
Rscript OSCER_urbmigmalavi_05pcms_noLeuco_F9.R