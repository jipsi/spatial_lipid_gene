#!/bin/bash
#SBATCH --job-name=sd22.5.1_sd22.8.4
#SBATCH --mail-type=END                   
#SBATCH --mail-user=shoumit.dey@york.ac.uk          
#SBATCH --ntasks=1                             
#SBATCH --cpus-per-task=1                      
#SBATCH --mem=40gb                           
#SBATCH --time=04:00:00                        
#SBATCH --output=sd22.5.3.nmf.log               
#SBATCH --account=hyms-htleish-2019              
#SBATCH --partition=gpu                       
#SBATCH --gres=gpu:1                           

module load system/CUDA/10.0.130
module load lang/Miniconda3/4.9.2

source activate cell2loc_env2

command -v python 

python nmf_compartments.py 

source deactivate 
