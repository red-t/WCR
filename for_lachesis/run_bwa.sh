#!/bin/bash

######## Help Information ########
function help_info(){
    echo `basename $0`
    echo -e "\t-1 <fq/fq.gz>\tFirst FASTQ."
    echo -e "\t-2 <fq/fq.gz>\tSecond FASTQ."
    echo -e "\t-I <idx prefix>\tPrefix of BWA Index."
    echo -e "\t-h \tShow this information"
}


######## Getting parameters ########
while getopts ":1:2:I:h" OPTION; do
    case $OPTION in
        1)  FIRST=$OPTARG;;
        2)  SECOND=$OPTARG;;
        I)  INDEX=$OPTARG;;
        h)  help_info && exit 1;;
        *)  help_info && exit 1;;
    esac
done


DIR=`dirname ${FIRST%.*1*f*q}`
PREFIX=`basename ${FIRST%.*1*f*q}`

### checking ###
echo -e "which BWA:\t`which bwa`"
echo -e "FIRST:\t$FIRST"
echo -e "SECOND:\t$SECOND"
echo -e "DIR:\t$DIR"
echo -e "Align command:\tbwa -PS5 -t 32 $INDEX $FIRST $SECOND -o $DIR/$PREFIX.sam"
### checking ###

bwa mem -PS5 -t 32 $INDEX $FIRST $SECOND -o $DIR/$PREFIX.sam