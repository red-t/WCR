#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=20G
#SBATCH -c 36
#SBATCH --array=42-200%5
#SBATCH --partition=12hours
#SBATCH --output=/data/tusers.ds/zhongrenhu/WCR/logs/blastn-log-%A-%a


# Print a little info for debugging
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


# Parameters initiation
BLAST_IMAGE=/data/tusers/zhongrenhu/Software/singularity_images/blastn_lynnlangit.sif
OUTPUTDIR=/data/tusers/zhongrenhu/WCR/draft_assembly/blast


# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""


# Copy inputs to temp dir
echo "Copying blastn database to $TMPDIR"
cp /data/tusers/zhongrenhu/WCR/draft_assembly/blast/GCF_917563875.1_PGI_DIABVI_V3a_genomic.fasta* .

echo "Copying query inputs to $TMPDIR"
if [ $SLURM_ARRAY_TASK_ID -lt 10 ];then
    QUERY=wcr_v0.3.0.part-00$SLURM_ARRAY_TASK_ID.fasta
elif [ $SLURM_ARRAY_TASK_ID -lt 100 ];then
    QUERY=wcr_v0.3.0.part-0$SLURM_ARRAY_TASK_ID.fasta
else
    QUERY=wcr_v0.3.0.part-$SLURM_ARRAY_TASK_ID.fasta
fi
cp /data/tusers/zhongrenhu/WCR/draft_assembly/blast/$QUERY .
echo "Done."
echo ""


# Process the data
echo "Running BLASTN for $QUERY" && date
singularity exec -B ${PWD}:/home/${USER} ${BLAST_IMAGE} blastn -num_threads 24 -query ${QUERY} -db GCF_917563875.1_PGI_DIABVI_V3a_genomic.fasta.blastdb -perc_identity 95 -evalue 1e-30 -word_size 50 -out wcr_v0.3.0.$SLURM_ARRAY_TASK_ID.blast.out -outfmt 7
echo "Done."
echo ""
date


# Copy files to output dir
echo "Copying results to destination..."
cp wcr_v0.3.0.$SLURM_ARRAY_TASK_ID.blast.out $OUTPUTDIR/
echo "Done."
echo ""


# Clean up
echo "Cleaning up..."
cd $HOME
echo "Deleting temp dir: " $TMPDIR
rm -rd $TMPDIR
echo "/tmp contents:"
ls -lh /tmp
echo ""
echo "Script complete."
date