#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --mail-user=jcorn427@uw.edu
#SBATCH --mail-type=END
#SBATCH -J Age_02_QC

cd /vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/trimgalore_results/

PATH_FASTQ="/vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/trimgalore_results/"
samples=("$PATH_FASTQ"*.fq.gz)

mkdir trimgalore_fastqc_results

for sample in  ${samples[*]}
do
	srun -c 16 /vol01/ngs_tools/FastQC/fastqc "$sample" --noextract -t 16 -o ./trimgalore_fastqc_results/
	wait
done