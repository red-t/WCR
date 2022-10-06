#!/bin/bash

# usage: run_blastn.sh -d /data/tusers/zhongrenhu/WCR/draft_assembly -i /data/tusers/zhongrenhu/Software/singularity_images/blastn_lynnlangit.sif -r /data/tusers.ds/zhongrenhu/WCR/draft_assembly/GCF_917563875.1_PGI_DIABVI_V3a_genomic.fasta -q /data/tusers.ds/zhongrenhu/WCR/draft_assembly/wcr_v0.3.0.fasta

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-d <work_dir>\tworking directory."
    echo -e "\t-i <image>\tpath of blast image(singularity)."
    echo -e "\t-r <ref_fa>\tpath of reference fasta."
    echo -e "\t-q <query_fa>\tpath of query fasta."
    echo -e "\t-h \tShow this information"
}
######## Help Information ########


######## Getting parameters ########
while getopts ":d:i:r:q:h" OPTION; do
    case $OPTION in
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
######## Getting parameters ########


### checking ###
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
[ ! -f ${REF_NAME}.blastdb.nhr ] && singularity exec -B ${PWD}:/home/${USER} ${BLAST_IMAGE} makeblastdb -in ${REF_NAME} -dbtype=nucl -out ${REF_NAME}.blastdb -logfile log_blastn
[ ! -f ${QUERY_NAME}.blast.out ] && singularity exec -B ${PWD}:/home/${USER} ${BLAST_IMAGE} blastn -query ${QUERY_NAME} -db ${REF_NAME}.blastdb -perc_identity 95 -evalue 1e-30 -word_size 50 -out ${QUERY_NAME}.blast.out -outfmt 7 &
### Main ###