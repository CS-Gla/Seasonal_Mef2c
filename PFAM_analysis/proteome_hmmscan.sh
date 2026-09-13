#PBS -l cput=100:00:00

#PBS -l walltime=100:00:00

#PBS -l nodes=1:ppn=16


#PBS -M calum.stewart.2@glasgow.ac.uk

#PBS -m abe
source /export/home2/cs424d/*/miniconda3/etc/profile.d/conda.sh

cd errata

cd data

conda activate hmmer

mkdir -p pfam_scan

for faa in prots/*.fasta; do
    base=$(basename "$faa" .fasta)
    
    hmmscan --cpu 16 \
        --domtblout pfam_scan/${base}.domtblout \
        Pfam-A.hmm \
        "$faa" \
        > pfam_scan/${base}.hmmscan.txt
done
