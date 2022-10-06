#!/bin/bash

# usage: run_blastn.sh -N 200 -d /data/tusers/zhongrenhu/WCR/draft_assembly -i /data/tusers/zhongrenhu/Software/singularity_images/blastn_lynnlangit.sif -r /data/tusers.ds/zhongrenhu/WCR/draft_assembly/GCF_917563875.1_PGI_DIABVI_V3a_genomic.fasta -q /data/tusers.ds/zhongrenhu/WCR/draft_assembly/wcr_v0.3.0.fasta

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-N <N_chunk>\tsplit the input file into N_chunk pieces."
    echo -e "\t-d <work_dir>\tworking directory."
    echo -e "\t-i <image>\tpath of blast image(singularity)."
    echo -e "\t-r <ref_fa>\tpath of reference fasta."
    echo -e "\t-q <query_fa>\tpath of query fasta."
    echo -e "\t-h \tShow this information"
}
######## Help Information ########


######## Getting parameters ########
while getopts ":N:d:i:r:q:h" OPTION; do
    case $OPTION in
        N)  N_CHUNK=$OPTARG;;
        d)  WORK_DIR=$OPTARG;;
        i)  BLAST_IMAGE=$OPTARG;;
        r)  REF_FA=$OPTARG;;
        q)  QUERY_FA=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done

REF_NAME=`basename ${REF_FA}`
QUERY_NAME=`basename ${QUERY_FA}`
PREFIX=`${QUERY_NAME%.*}`
SUFFIX=`${QUERY_NAME##*.}`
######## Getting parameters ########


### checking ###
echo -e "N_CHUNK:\t${N_CHUNK}"
echo -e "WORK_DIR:\t${WORK_DIR}"
echo -e "BLAST_IMAGE:\t${BLAST_IMAGE}"
echo -e "REF_FA:\t${REF_FA}"
echo -e "QUERY_FA:\t${QUERY_FA}"

echo -e "ENTER WORK DIRECTORY:  ${WORK_DIR}"
cd ${WORK_DIR}

echo -e "COPY INPUT FILE INTO WORK DIRECTORY"
[ ! -f ${REF_NAME} ] && cp ${REF_FA} .
[ ! -f ${QUERY_NAME} ] && cp ${QUERY_FA} .
### checking ###


### Main ###

# build blastn database
[ ! -f ${REF_NAME}.blastdb.nhr ] && singularity exec -B ${PWD}:/home/${USER} ${BLAST_IMAGE} makeblastdb -in ${REF_NAME} -dbtype=nucl -out ${REF_NAME}.blastdb -logfile log_blastn

# split input file
fasta-splitter.pl --n-parts ${N_CHUNK} ${REF_NAME}

# run blast
for CHUNK in {1..${N_CHUNK}};do
    [ ! -f ${PREFIX}.${CHUNK}.blast.out ] && singularity exec -B ${PWD}:/home/${USER} ${BLAST_IMAGE} blastn -query ${PREFIX}.part-${CHUNK}.${SUFFIX} -db ${REF_NAME}.blastdb -perc_identity 95 -evalue 1e-30 -word_size 50 -out ${PREFIX}.${CHUNK}.blast.out -outfmt 7 &
done

### Main ###