#!/bin/bash

#sub script for gunzip_and_cutadapt_mirna_SL_list_in.sh
#7/30/2014 mod to do simultaneous 5' adaptor trimming on mod v3 adaptors that are all the same length.

GZ=$1
slid=$2

export PATH=/share/samdata/sorenar/miRNA/pipeline/src/cutadapt-1.7.1/bin/:$PATH
result_path="/share/samdata/sorenar/miRNA/pipeline/results/"$slid"/"


THREE_PRIME_AD_SEQ="ACGGGCTAATATTTATCGGTGGAGCATCACGATCTCGTAT"
FIVE_PRIME_AD1_SEQ="^CAGTCG"
FIVE_PRIME_AD2_SEQ="^TGACTC"
FIVE_PRIME_AD3_SEQ="^GCTAGA"
FIVE_PRIME_AD4_SEQ="^ATCGAT"


SL_TRIM_FILE=$result_path$slid"_trim.fastq"
NO_3AD_FILE=$result_path$slid"_NO3AD.fastq"
NO_5AD_FILE=$result_path$slid"_NO5AD.fastq"
TOO_SHORT_FILE=$result_path$slid"_SHORT_FAIL.fastq"

echo $SL_TRIM_FILE
echo $NO_3AD_FILE
echo $NO_5AD_FILE
echo $TOO_SHORT_FILE

cutadapt -a $THREE_PRIME_AD_SEQ -e 0.25 --match-read-wildcards --untrimmed-output=$NO_3AD_FILE $GZ | cutadapt -e 0.34 --match-read-wildcards --no-indels -m 15 -O 6 -n 1 -g $FIVE_PRIME_AD1_SEQ -g $FIVE_PRIME_AD2_SEQ -g $FIVE_PRIME_AD3_SEQ -g $FIVE_PRIME_AD4_SEQ --untrimmed-output=$NO_5AD_FILE --too-short-output=$TOO_SHORT_FILE - > $SL_TRIM_FILE
