#!/bin/bash
# Creates and copies all figures that are relevant into a folder
basedir=/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots
figfold=$basedir/FinalPaperPlots
# rm $figfold/*
plotfold=$basedir/plots

########################################################################################################
## Figure 1: Magnetization Evolution
########################################################################################################
filepath=Magnetization_Evolution/EvolutionOnSphere/mxy_M_Evolution
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 2: Magnetization vs temperature
########################################################################################################
# Plot_Magnetization.m
filepath=static/mxy_Magnetization
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 3: Susceptibility
########################################################################################################
#filepath=static/mxy_AbsSusceptibility_oldversion
filepath=static/mxy_AbsSusceptibility
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 4: Magnetization histogram
########################################################################################################
filepath=Mhistogram/mxy_Mhistogram_T_.17
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 5: structural properties
########################################################################################################
#filepath=structure/mxy_strutcutre_sqrtN_256
#cp "$plotfold/$filepath".* $figfold
filepath=/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2022_ThesisPlots/gr
cp "$filepath"/*.pdf $figfold
cp "$filepath"/*.png $figfold

########################################################################################################
## Figure 5 cont: G6
########################################################################################################
# /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder/Plot_G6_Paper.m
filepath=/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/BondOrder/G6/mxy_G6_sqrtN_128
cp "$filepath".* $figfold


########################################################################################################
## Figure 6: spin correlation function
########################################################################################################
filepath=Spin_SCF/mxy_SpinSCF_sqrtN_256
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 7a: spin correlation function FS
########################################################################################################
filepath=Spin_SCF/mxy_SpinSCF_FS_3in1
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 7b: eta analysis
########################################################################################################
filepath=Spin_SCF/mxy_Eta_Compare
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 8: chimq Multiplot
########################################################################################################
filepath=chi_plots/Multiplot/mxy_chi_Multiplot
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 9: Lepri-Ruffo short-time behavior
########################################################################################################
filepath=LepriRuffo/LepriRuffo_ShortTime_Multiplot_T_0.170
cp "$plotfold/$filepath".* $figfold
#filepath=LepriRuffo/mxy_LepriRuffo_ShortTime_T_0.170
#cp "$plotfold/$filepath".* $figfold
#filepath=LepriRuffo/fmxy_LepriRuffo_ShortTime_T_0.170
#cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 10: Lepri-Ruffo long-time multiplot
########################################################################################################
filepath=LepriRuffo/ACF_Multiplot
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 11: Multiplot of Spin TCF properties
########################################################################################################
# Plot_TCF_FFT_Multiplot.m
filepath=TCF/Multiplot/Compare_Corrs_sqrtN_128_q_0.255
#filepath=TCF/Multiplot/Compare_Corrs_sqrtN_128_q_0.170
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 12: Nelson-Fisher law for high frequencies
########################################################################################################
# NF_FFT/NF_FFT_PictureCreate.m
#filepath=TCF/NF_HighFreq/mxy_NFHighFreq_T_0.17_q_0.34_spline_factor_0_smoothened
filepath="/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/NF_FFT/plots/mxy_Smperp_T_0.140_q_0.679"
#TCF/NF_HighFreq/mxy_NFHighFreq_T_0.17_q_0.34_spline_factor_0_smoothened
cp "$filepath".* $figfold

########################################################################################################
## Figure 13: Finite-size scaling of peak height (Probably discarded)
########################################################################################################
filepath=TCF/FFT_FS/mperp_mxy_TCF_FFT_FS_T_0.170_nq_3
cp "$plotfold/$filepath".* $figfold
filepath=TCF/FFT_FS/mperp_fmxy_TCF_FFT_FS_T_0.170_nq_3
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 14: Spin wave propagation speed
########################################################################################################
filepath=omega_gamma_plots/mxy_c_spinwavespeed_omega_1
cp "$plotfold/$filepath".* $figfold


########################################################################################################
## Figure 15: Omega & gamma plots
########################################################################################################
FILEPATHS=("omega_gamma_plots/mxy_omega_1_vs_T"
  "omega_gamma_plots/mxy_gamma_vs_T"
  "omega_gamma_plots/mxy_omega_1_vs_q_log"
  "omega_gamma_plots/mxy_gamma_vs_q_log")
for filepath in ${FILEPATHS[@]}; do
  cp "$plotfold/$filepath".* $figfold
done

########################################################################################################
## Figure 15 Replacement: Omega & gamma plots: Fit compare
########################################################################################################
filepath=omega_gamma_plots/mxy_omegagamma_Compare_Fits_sqrtN_128_T_0.170
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 15 Replacement: Omega & gamma plots: Size compare
########################################################################################################
filepath=omega_gamma_plots/mxy_omegagamma_Compare_Size_T_0.170
cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 16 Replacement: Omega & gamma plots: Spin wave and gamma
########################################################################################################
#filepath=omega_gamma_plots/mxy_omegagamma_SpinwaveGammaexp
filepath=/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/SymmAntisymmFit/mxy_omegagamma_SpinwaveGammaexp
cp "$filepath".* $figfold

#cp "$plotfold/$filepath".* $figfold

########################################################################################################
## Figure 17 Replacement: c^2 by Upsilon
########################################################################################################

filepath=/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/plots/omega_gamma_plots/mxy_omegagamma_c_by_Upsilon
cp "$filepath".* $figfold


cd $figfold
# Creating tiff files
ls *.pdf | while read line
do
  name=${line/.pdf/}
  echo $line
  pdftoppm -tiff -r 300 $name.pdf tiff_files/$name
  mv "tiff_files/$name-1.tif" tiff_files/$name.tif
done

cd -
