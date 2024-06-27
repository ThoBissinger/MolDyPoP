#!/bin/bash
################################################################################
## NAMEVEC_FROM_SAMPFILE.SH
## Extracts names of all variables stored in a sampfile. Can be used to
## copy-paste variable names to matlab
## Input parameters
## 1         sampfile      Path to sampfile
##
## Call as
## >> namevec_from_sampfile.sh "SAMPFILE"
################################################################################
sampfile=$1
cat $sampfile | cut -d '=' -f 1 | while read line
do
	echo "\"$line\""
done
