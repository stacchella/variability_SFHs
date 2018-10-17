#!/bin/bash
### Name of the job
### Requested number of cores
#SBATCH -n 1
### Requested number of nodes
#SBATCH -N 1
### Requested computing time in minutes
#SBATCH -t 10080
### Partition or queue name
#SBATCH -p conroy,itc_cluster,hernquist
### memory per cpu, in MB
#SBATCH --mem-per-cpu=4000
### Job name
#SBATCH -J 'VARSFH_standard'
### output and error logs
#SBATCH -o VARSFH_standard_%a.out
#SBATCH -e VARSFH_standard_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
srun -n 1 python /n/conroyfs1/stacchella/variability_SFHs/scripts/run_make_SFH.py \
--idx_key="${SLURM_ARRAY_TASK_ID}" \
--filename_SFH="SFH_" \
--sfh_res=0.001 \
--redshift_start=0.3 \
--redshift_end=0.0 \
--redshift_lum_start=0.2 \
--redshift_lum_end=0.0 \
--number_galaxies=10 \
--scatter_MS_0=1.0 \
--list_of_slopes 0.0 0.5  \
--list_of_breaks 10.0 30.0 \
--aliasTbin=0.1 \
--logzsol=0.0 \
--dust2=0.2 \





