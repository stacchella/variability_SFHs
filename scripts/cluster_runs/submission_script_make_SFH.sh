#!/bin/bash
### Name of the job
### Requested number of cores
#SBATCH -n 1
### Requested number of nodes
#SBATCH -N 1
### Requested computing time in minutes
#SBATCH -t 10080
### Partition or queue name
#SBATCH -p conroy-intel,conroy,itc_cluster,general,hernquist,serial_requeue
### memory per cpu, in MB
#SBATCH --mem-per-cpu=4000
### Job name
#SBATCH -J 'SFH_z4_Z_fid_02_withoutT'
### output and error logs
#SBATCH -o SFH_z4_Z_fid_02_%a.out
#SBATCH -e SFH_z4_Z_fid_02_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
slopes = [0.0, 1.0, 1.5, 2.0, 2.5, 3.0]
breaks = [10.0, 30.0, 100.0, 300.0, 1000.0]
srun -n 1 python /n/conroyfs1/stacchella/variability_SFHs/scripts/run_make_SFH.py \
--idx_key="${SLURM_ARRAY_TASK_ID}" \
--filename_SFH="SFH_" \
--sfh_res=0.001 \
--redshift_start=0.5 \
--redshift_end=0.0 \
--number_of_galaxies=10 \
--scatter_MS_0=1.0 \
--list_of_slopes=slopes\
--list_of_breaks=breaks \
--aliasTbin=0.1 \







