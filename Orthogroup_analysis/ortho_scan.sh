#PBS -l cput=100:00:00

#PBS -l walltime=100:00:00

#PBS -l nodes=1:ppn=16:oracle9

#PBS -M calum.stewart.2@glasgow.ac.uk

#PBS -m abe

source /export/home2/cs424d/*/miniconda3/etc/profile.d/conda.sh

cd ortho

conda activate ortho39

orthofinder -f proteomes/
