#!/bin/bash
################################################################################
## INFILE_CHANGE_OPTION.SH
## Changes an option in the infile
## Input parameters
## 1         infile        Path to infile that should be changed
## 2         optname       Name of the option to be changed
## 3         optval        Value of the option to be changed
##
## Call as
## >> infile_change_option.sh "INFILE" "OPTNAME" "OPTVAL"
################################################################################
infile=$1
optname=$2
optval=$3
if [[ ! $3 ]]; then
  echo 'Input for infile_change_option:
$1: infile   (path to infile)
$2: optname  (name of option to be changed)
$3: optval   (new value for the option)'
else
  ## If input is correct, execute this sed command
  sed -i "/^$optname =/ c $optname = $optval" $infile
fi
