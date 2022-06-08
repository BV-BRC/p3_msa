#!/bin/sh
# $1 test case
MODULE=~/dev_container/modules/p3_msa
~/dev_container/bin/App-MSA https://p3.theseed.org/services/app_service/ $MODULE/app_specs/MSA.json $MODULE/tests/$1.json 1> $MODULE/tests/out/$1.out 2> $MODULE/tests/out/$1.err
