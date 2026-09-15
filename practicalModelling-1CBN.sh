#!/bin/bash
#
#SBATCH --job-name=NeutralNetworkProtein1CBN
#SBATCH --output=NeutralNetworkProtein1CBN.out
#SBATCH --error=NeutralNetworkProtein1CBN.err
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT
#SBATCH --cpus-per-task=8
#SBATCH --mem=3000M
#SBATCH --time=00-24:00:00

# Exit the slurm script if a command fails
set -e

# Run the assembly
python /lisc/home/user/vester/PracticalModellingProject/NeutralNetworkProtein.py ./1CBN_Results/ ./1CBN.pdb \
    --threads "$SLURM_CPUS_PER_TASK" \
    --memory "$(($SLURM_MEM_PER_NODE / 1000))" \
    --tmp-dir "$TMPDIR" \
    --test