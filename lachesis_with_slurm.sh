#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=40G
#SBATCH -c 30
#SBATCH --array=1001-2000%30
#SBATCH --partition=4hours
#SBATCH --output=/data/tusers.ds/zhongrenhu/WCR/logs/lachesis-log-%A-%a

MINWAIT=120
MAXWAIT=240
sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))

# Random integer function
function rand_int() {
    SEED=$RANDOM
    RND_NUM=`echo $SEED $1 $2|awk '{srand($1);printf "%d",rand()*10000%($3-$2)+$2}'`
    echo $RND_NUM
}
# rnd=$(rand_int 1 50)

# Random float function
function rand_float() {
    NUM=$(rand_int $1 $2)
    RND_NUM=`awk -v NUM=$NUM -v SEED=$RANDOM 'BEGIN { srand(SEED);printf("%.3f\n", NUM+rand()) }'`
    echo $RND_NUM
}


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

echo "Copying cached RE site counts information to $TMPDIR/draft_assembly ..."
cp /data/tusers/zhongrenhu/WCR/draft_assembly/*counts_AAGCTT* ./draft_assembly

echo "Copying cached contig names to $TMPDIR/draft_assembly ..."
cp /data/tusers/zhongrenhu/WCR/draft_assembly/*names ./draft_assembly

#echo "Copying cached draft-assembly contigs & reference genome to $TMPDIR/draft_assembly ..."
#cp /data/tusers/zhongrenhu/WCR/draft_assembly/*fasta ./draft_assembly

echo "Copying Hi-C alignments to $TMPDIR/HiC ..."
cp /data/tusers/zhongrenhu/WCR/HiC/HiC-WCR-ov2-10cycle.REduced.paired_only.bam ./HiC
echo "Done."
echo ""


# INI file initialization
touch round$SLURM_ARRAY_TASK_ID.ini
INI=round$SLURM_ARRAY_TASK_ID.ini
CLUSTER_MIN_RE_SITES=$(rand_int 8 20)
CLUSTER_MAX_LINK_DENSITY=$(rand_float 2 3)
CLUSTER_NONINFORMATIVE_RATIO=$(rand_float 3 5)
ORDER_MIN_N_RES_IN_TRUNK=$(rand_int 30 100)
ORDER_MIN_N_RES_IN_SHREDS=$(rand_int 8 20)
REPORT_QUALITY_FILTER=2

echo -e "SPECIES = wcr\n" >> $INI
echo -e "OUTPUT_DIR = round$SLURM_ARRAY_TASK_ID\n" >> $INI
echo -e "DRAFT_ASSEMBLY_FASTA = draft_assembly/wcr_v0.3.0.fasta\n" >> $INI
echo -e "SAM_DIR = HiC\n" >> $INI
echo -e "SAM_FILES = HiC-WCR-ov2-10cycle.REduced.paired_only.bam\n" >> $INI
echo -e "RE_SITE_SEQ = AAGCTT\n" >> $INI
echo -e "USE_REFERENCE = 1\n" >> $INI
echo -e "SIM_BIN_SIZE = 0\n" >> $INI
echo -e "REF_ASSEMBLY_FASTA = draft_assembly/GCF_917563875.1_PGI_DIABVI_V3a_genomic.fasta\n" >> $INI
echo -e "BLAST_FILE_HEAD = draft_assembly/wcr_v0.3.0\n" >> $INI
echo -e "DO_CLUSTERING = 1\n" >> $INI
echo -e "DO_ORDERING   = 1\n" >> $INI
echo -e "DO_REPORTING  = 1\n" >> $INI
echo -e "OVERWRITE_GLM = 0\n" >> $INI
echo -e "OVERWRITE_CLMS = 0\n" >> $INI
echo -e "CLUSTER_N = 10\n" >> $INI
echo -e "CLUSTER_CONTIGS_WITH_CENS = -1\n" >> $INI
echo -e "CLUSTER_MIN_RE_SITES = $CLUSTER_MIN_RE_SITES\n" >> $INI
echo -e "CLUSTER_MAX_LINK_DENSITY = $CLUSTER_MAX_LINK_DENSITY\n" >> $INI
echo -e "CLUSTER_NONINFORMATIVE_RATIO = $CLUSTER_NONINFORMATIVE_RATIO\n" >> $INI
echo -e "CLUSTER_DRAW_HEATMAP = 1\n" >> $INI
echo -e "CLUSTER_DRAW_DOTPLOT = 1\n" >> $INI
echo -e "ORDER_MIN_N_RES_IN_TRUNK = $ORDER_MIN_N_RES_IN_TRUNK\n" >> $INI
echo -e "ORDER_MIN_N_RES_IN_SHREDS = $ORDER_MIN_N_RES_IN_SHREDS\n" >> $INI
echo -e "ORDER_DRAW_DOTPLOTS = 1\n" >> $INI
echo -e "REPORT_EXCLUDED_GROUPS = -1\n" >> $INI
echo -e "REPORT_QUALITY_FILTER = $REPORT_QUALITY_FILTER\n" >> $INI
echo -e "REPORT_DRAW_HEATMAP = 1\n" >> $INI


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