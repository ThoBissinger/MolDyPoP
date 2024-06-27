#!/bin/bash

basedir="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/FinalPaperPlots"
pdfdir=$basedir
tifdir=$basedir/tiff_files
newdir=$basedir/RenamedPlots
NAMES=(mxy_M_Evolution
mxy_Magnetization
mxy_AbsSusceptibility
mxy_Mhistogram_T_.17
mxy_gr_sqrtN_256
mxy_DiffCoeff_sqrtN_256
mxy_snap_sqrtN_256_T_0.010
mxy_snap_sqrtN_256_T_0.090
mxy_SpinSCF_sqrtN_256
mxy_SpinSCF_FS_3in1
mxy_chi_Multiplot
mxy_Eta_Compare
LepriRuffo_ShortTime_Multiplot_T_0.170
ACF_Multiplot
Compare_Corrs_sqrtN_128_q_0.255
mxy_Smperp_T_0.140_q_0.679
mxy_omegagamma_Compare_Fits_sqrtN_128_T_0.170
mxy_omegagamma_SpinwaveGammaexp)

N_N=${#NAMES[@]}

NewNAMES=(Fig_1
Fig_2
Fig_3
Fig_4
Fig_5_a
Fig_5_b
Fig_5_c
Fig_5_d
Fig_6_a
Fig_6_b
Fig_7
Fig_8
Fig_9
Fig_10
Fig_11
Fig_12
Fig_13
Fig_14)



for i in `seq 0 $((N_N-1))`; do
  oldname=${NAMES[$i]}
  newname=${NewNAMES[$i]}
  cp "$pdfdir/$oldname".pdf "$newdir/$newname."pdf
  cp "$tifdir/$oldname".tif "$newdir/$newname."tif
done

for i in `seq 0 $((N_N-1))`; do
  oldname=${NAMES[$i]}
  newname=${NewNAMES[$i]}
  #echo "$oldname $(grep $oldname $basedir/manuscript_submission.tex)"
  sed -i "s/\.\/img\/$oldname/$newname/" $basedir/manuscript_submission.tex
done
