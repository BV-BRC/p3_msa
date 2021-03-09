#!/bin/sh
# $1 test case
perl App-MSA.pl https://p3.theseed.org/services/app_service/ ~/p3_msa/p3_msa/app_specs/MSA.json ~/p3_msa/p3_msa/tests/$1.json 1> ~/p3_msa/tests/$1.out 2> ~/p3_msa/tests/$1.err
