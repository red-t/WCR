#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=20G
#SBATCH -c 30
#SBATCH --array=2-200%5
#SBATCH --partition=12hours
#SBATCH --output=/data/tusers.ds/zhongrenhu/WCR/logs/lachesis-log-%A-%a


# Print a little info for debugging
echo "HOSTNAME: " $(hostname)
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
date
echo ""


# Create temporary directory for work
echo "Creating temporary working dir..."
TMPDIR=$(mktemp -d -t zhongrenhu-tmp-XXXXXXXX)
cd $TMPDIR
echo "Changing wd: " $(pwd)
echo ""


# Workdir initialization
mkdir -p draft_assembly
mkdir -p round$SLURM_ARRAY_TASK_ID/cached_data
mkdir -p HiC
OUTPUTDIR=/data/tusers/zhongrenhu/WCR/lachesis_scaffold


# Copy inputs to temp dir
echo "Copying cached GLM file to $TMPDIR/round$SLURM_ARRAY_TASK_ID/cached_data ..."
cp /data/tusers/zhongrenhu/WCR/lachesis_scaffold/round1/cached_data/all.GLM round$SLURM_ARRAY_TASK_ID/cached_data

echo "Copying cached TrueMapping file to $TMPDIR/round$SLURM_ARRAY_TASK_ID/cached_data ..."
cp /data/tusers/zhongrenhu/WCR/lachesis_scaffold/round1/cached_data/TrueMapping.assembly.txt round$SLURM_ARRAY_TASK_ID/cached_data

echo "Copying cached RE sites information to $TMPDIR/draft_assembly ..."
cp /data/tusers/zhongrenhu/WCR/draft_assembly/*AAGCTT* ./draft_assembly

echo "Copying cached contig names to $TMPDIR/draft_assembly ..."
cp /data/tusers/zhongrenhu/WCR/draft_assembly/*names ./draft_assembly

echo "Copying cached draft-assembly contigs & reference genome to $TMPDIR/draft_assembly ..."
cp /data/tusers/zhongrenhu/WCR/draft_assembly/*fasta ./draft_assembly

echo "Copying Hi-C alignments to $TMPDIR/HiC ..."
cp /data/tusers/zhongrenhu/WCR/HiC/HiC-WCR-ov2-10cycle.REduced.paired_only.bam ./HiC
echo "Done."
echo ""


# INI file initialization
INI=$(touch round$SLURM_ARRAY_TASK_ID.ini)
# not finished yet


# Process the data
echo "Running LACHESIS" && date
singularity exec -B $PWD:/home/$USER /data/tusers/zhongrenhu/Software/singularity_images/lachesis_huzr14830.sif Lachesis $INI
echo "Done."
echo ""
date


# Copy files to output dir
echo "Copying results to destination..."
mv ./*jpg round$SLURM_ARRAY_TASK_ID
mv ./*txt round$SLURM_ARRAY_TASK_ID
mv $INI round$SLURM_ARRAY_TASK_ID
rm -r round$SLURM_ARRAY_TASK_ID/cached_data
cp -a round$SLURM_ARRAY_TASK_ID $OUTPUTDIR/
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