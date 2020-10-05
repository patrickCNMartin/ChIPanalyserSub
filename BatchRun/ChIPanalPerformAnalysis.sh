#!/bin/bash

#
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -M pm16057@essex.ac.uk

data=$'DataInputNULLgeometric.txt'

arg1=$(sed "$1q;d" $data | awk '{print $1}')
arg2=$(sed "$1q;d" $data | awk '{print $2}')
arg3=$(sed "$1q;d" $data | awk '{print $3}')
arg4=$(sed "$1q;d" $data | awk '{print $4}')
arg5=$(sed "$1q;d" $data | awk '{print $5}')
arg6=$(sed "$1q;d" $data | awk '{print $6}')
arg7=$(sed "$1q;d" $data | awk '{print $7}')
arg8=$(sed "$1q;d" $data | awk '{print $8}')
arg9=$(sed "$1q;d" $data | awk '{print $9}')
arg10=$(sed "$1q;d" $data | awk '{print $10}')
arg11=$(sed "$1q;d" $data | awk '{print $11}')
arg12=$(sed "$1q;d" $data | awk '{print $12}')
arg13=$(sed "$1q;d" $data | awk '{print $13}')
arg14=$(sed "$1q;d" $data | awk '{print $14}')
arg15=$(sed "$1q;d" $data | awk '{print $15}')


./ChIPanalPerformAnalysis.R $arg1 $arg2 $arg3 $arg4 $arg5 $arg6 $arg7 $arg8 $arg9 $arg10 $arg11 $arg12 $arg13 $arg14 $arg15
