#!/bin/sh
# $1 test case
./run.sh aa_align
rm -r  stage/ work/
./run.sh aa_align_fasta
rm -r  stage/ work/
./run.sh aa_align_fg
rm -r  stage/ work/
./run.sh aa_feature_group
rm -r  stage/ work/
./run.sh convert_input
rm -r  stage/ work/
./run.sh mixed_all_protein
rm -r  stage/ work/
./run.sh nt_all
rm -r  stage/ work/
./run.sh protein_input
rm -r  stage/ work/
./run.sh aa_all
rm -r  stage/ work/
./run.sh aa_hemoglobin
rm -r  stage/ work/
./run.sh convert_input_align
rm -r  stage/ work/
./run.sh mixed_align_fasta_fg
rm -r  stage/ work/
./run.sh mixed_fasta_fg
rm -r  stage/ work/
./run.sh nt_feature_group
rm -r  stage/ work/
./run.sh genome_group_covid19
rm -r stage/ work/