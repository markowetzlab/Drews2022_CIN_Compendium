#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=MC_Simulation
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ruben.drews@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH --output=/home/%u/data/phd/cnsigs2/cnsigs2_revisions/slurm.%A.%x-%a.out
#SBATCH --error=/home/%u/data/phd/cnsigs2/cnsigs2_revisions/slurm.%A.%x-%a.err

###############################################################
### Execute script                                         ####
###############################################################

## Call this script on the headnode with the following command:
## sbatch -N 1 --mem 30GB /Users/drews01/phd/prjcts/cnsigs2/cnsigs2_revisions/scripts/Sig_activity_with_noise.sh

echo -e "JobID: $SLURM_JOB_ID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

BASE="/Users/drews01/data/phd/cnsigs2/cnsigs2_revisions"
UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)

Rscript --vanilla ${BASE}/scripts/Sig_activity_with_noise.R $UUID



  
