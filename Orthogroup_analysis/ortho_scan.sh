#PBS -l cput=100:00:00

#PBS -l walltime=100:00:00

#PBS -l nodes=1:ppn=16:oracle9

#PBS -m abe

source /*/miniconda3/etc/profile.d/conda.sh

cd ortho

conda activate ortho39

orthofinder -f proteomes/
