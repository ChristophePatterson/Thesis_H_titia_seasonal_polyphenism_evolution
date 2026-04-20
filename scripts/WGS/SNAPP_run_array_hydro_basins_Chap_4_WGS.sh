#!/bin/bash
# Christophe Patterson
# 10/04/25
# for running on the UoN HPC Ada

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=84
#SBATCH --mem=80g
#SBATCH --time=167:00:00
#SBATCH --array=1
#SBATCH --job-name=SNAPP_WGS
#SBATCH --output=/gpfs01/home/mbzcp2/slurm_outputs/slurm-%x-%j.out

# Commands to execute start here
# OLD HAMILTON MODULES
###### module purge
###### module load gcc/native
###### module load gdal
###### module load geos/3.10.1
###### module load r/4.2.1
###### module load bioinformatics
# NEW ADA modules
module load R-uoneasy/4.2.1-foss-2022a
module load ruby-uoneasy/3.4.2-GCCcore-13.3.0

# module load beast
# May need to run the code
# /home/tmjj24/apps/beast/bin/packagemanager -add SNAPP

base_dir="/gpfs01/home/$USER/code/Github/Thesis_H_titia_seasonal_polyphenism_evolution"
cd $base_dir

hydro="5"
# Output directory
output_dir=(${base_dir}/data/SNAPP/Hydro_${hydro}_WGS_max_cov)
input_dir=(${base_dir}/data/SNPs/SNAPP)
mkdir -p $output_dir

# Input file prefix
phy_file=(titia.mysnps-SNAPP-hydro_${hydro}_max_cov)

# Run ruby script for creating xml SNAPP config file
# Uses https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md

ruby ${base_dir}/scripts/snapp_prep.rb -p $input_dir/$phy_file.phy -t $input_dir/$phy_file.txt \
-s $input_dir/${phy_file}_start_tree.newick -c $input_dir/$phy_file.con.txt -l 1000000 -x $output_dir/$phy_file.xml -o $output_dir/$phy_file
## ALERT I HAVE CHANGED THE LOCATION OF BEAST FROM HOME TO NOBACKUP SINCE RUNNING THIS SCRIPT AS SUCH PREPARE FOR ERRORS

~/apps/beast/bin/beast -threads $SLURM_CPUS_PER_TASK -overwrite $output_dir/$phy_file.xml > $output_dir/$phy_file.screen.log

## Rscript /gpfs01/home/mbzcp2/code/Github/Thesis_H_titia_seasonal_polyphenism_evolution/scripts/SNAPP/SNAPP_tree_QC.R $output_dir/$phy_file

~/apps/beast/bin/treeannotator -burnin 20 $output_dir/$phy_file.trees $output_dir/$phy_file._BI20.Anon.nex
~/apps/beast/bin/treeannotator -burnin 20 -height median $output_dir/$phy_file.trees $output_dir/$phy_file._BI20_HGTmed.Anon.nex

cd ~