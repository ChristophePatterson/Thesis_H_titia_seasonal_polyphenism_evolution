#!/bin/bash

#SBATCH -c 64 
#SBATCH --mem=15G            # memory required, in units of k,M or G, up to 250G.
#SBATCH --gres=tmp:15G       # $TMPDIR space required on each compute node, up to 400G.
#SBATCH -t 40:00:00         # time limit in format dd-hh:mm:ss
#SBATCH --output=/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC/Hetaerina_titia_ddRAD_titia_dg/SNAPP/slurm-%x.%j.out

# Commands to execute start here
module purge
module load gcc/native
module load gdal
module load geos/3.10.1
module load r/4.2.1
module load bioinformatics
# module load beast
# May need to run the code
# /home/tmjj24/apps/beast/bin/packagemanager -add SNAPP

# Output directory
output_dir=(/nobackup/tmjj24/ddRAD/Demultiplexed_seq_processing/SNP_libraries_SDC_3/Hetaerina_titia_ddRAD_titia_dg/SNAPP/Hydro_5_max_cov)

# Input file prefix
phy_file=(Hetaerina_titia_ddRAD_titia_dg-SNAPP-hydro_5_max_cov)

cd $output_dir

# Run ruby script for creating xml SNAPP config file
# Uses https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md

ruby /home/tmjj24/scripts/job_scripts/Master-demulitiplex-scripts/ddRAD_Durham_and_Sheffield_CUAJ/SNAPP/snapp_prep.rb -p $phy_file.phy -t $phy_file.txt \
-c $phy_file.con.txt -l 1000000 -x $phy_file.xml -o $phy_file

## ALERT I HAVE CHANGED THE LOCATION OF BEAST FROM HOME TO NOBACKUP SINCE RUNNING THIS SCRIPT AS SUCH PREPARE FOR ERRORS

/nobackup/tmjj24/apps/beast/bin/beast -threads 64 -overwrite $phy_file.xml > $phy_file.screen.log

/nobackup/tmjj24/apps/beast/bin/treeannotator -burnin 10 $phy_file.trees $phy_file.trees.Anon

cd ~