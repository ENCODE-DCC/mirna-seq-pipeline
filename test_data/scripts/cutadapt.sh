#!/bin/bash

#sub script for gunzip_and_cutadapt_mirna_SL_list_in.sh
#7/30/2014 mod to do simultaneous 5' adaptor trimming on mod v3 adaptors that are all the same length.

slid=$1
setid=$2

export PATH=/share/samdata/sorenar/miRNA/pipeline/src/cutadapt-1.7.1/bin/:$PATH
result_path="/share/samdata/sorenar/miRNA/pipeline/results/"$slid"/"
file_path="/share/samdata/sorenar/miRNA/pipeline/data/"

if [ $setid -eq 1 ]; then
	AD1="^CAGTCG";AD2="^TGACTC";AD3="^GCTAGA";AD4="^ATCGAT"
elif [ $setid -eq 2 ]; then
	AD1="^TCGCAG";AD2="^CTCTGA";AD3="^AGAGCT";AD4="^GATATC"
elif [ $setid -eq 3 ]; then
	AD1="^CACGTG";AD2="^TGTACC";AD3="^GCGTAA";AD4="^ATACGT"
elif [ $setid -eq 4 ]; then
        AD1="^CAGCGT";AD2="^TGCTAC";AD3="^GCAGTA";AD4="^ATTACG"
elif [ $setid -eq 5 ]; then
        AD1="^AGCGTC";AD2="^GATCCT";AD3="^CTGAAG";AD4="^TCATGA"
elif [ $setid -eq 6 ]; then
        AD1="^GCATCG";AD2="^ATGCTC";AD3="^TGCAGA";AD4="^CATGAT"
elif [ $setid -eq 7 ]; then
        AD1="^GTCGCA";AD2="^ACTCTG";AD3="^TAGAGC";AD4="^CGATAT"
fi

THREE_PRIME_AD_SEQ="ACGGGCTAATATTTATCGGTGGAGCATCACGATCTCGTAT"
FIVE_PRIME_AD1_SEQ=$AD1
FIVE_PRIME_AD2_SEQ=$AD2
FIVE_PRIME_AD3_SEQ=$AD3
FIVE_PRIME_AD4_SEQ=$AD4


SL_TRIM_FILE=$result_path$slid"_trim.fastq"
NO_3AD_FILE=$result_path$slid"_NO3AD.fastq"
NO_5AD_FILE=$result_path$slid"_NO5AD.fastq"
TOO_SHORT_FILE=$result_path$slid"_SHORT_FAIL.fastq"

cutadapt -a $THREE_PRIME_AD_SEQ -e 0.25 --match-read-wildcards --untrimmed-output=$NO_3AD_FILE $file_path$slid".fastq" | cutadapt -e 0.34 --match-read-wildcards --no-indels -m 15 -O 6 -n 1 -g $FIVE_PRIME_AD1_SEQ -g $FIVE_PRIME_AD2_SEQ -g $FIVE_PRIME_AD3_SEQ -g $FIVE_PRIME_AD4_SEQ --untrimmed-output=$NO_5AD_FILE --too-short-output=$TOO_SHORT_FILE - > $SL_TRIM_FILE
