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
#SBATCH -J 'VARSFHc_standard'
### output and error logs
#SBATCH -o VARSFHc_standard_%a.out
#SBATCH -e VARSFHc_standard_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
srun -n 1 python /n/conroyfs1/stacchella/variability_SFHs/scripts/run_make_SFH_conv.py \
--idx_key="${SLURM_ARRAY_TASK_ID}" \
--filename_SFH="SFHc_" \
--sfh_res=0.001 \
--redshift_start=0.2 \
--redshift_end=0.0 \
--number_galaxies=1000 \
--scatter_MS_0=1.0 \
--list_of_slopes 0.0 0.5 1.0 1.2 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.8 3.0 3.5 4.0 \
--list_of_breaks 10.0 30.0 50.0 70.0 80.0 90.0 100.0 110.0 120.0 140.0 160.0 200.0 250.0 300.0 400.0 600.0 1000.0 \
--aliasTbin=1.0 \
--sub_res_time=2 \


