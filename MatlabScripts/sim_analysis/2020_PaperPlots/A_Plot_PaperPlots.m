% Runs all plotting functions relevant for the paper

% Figure 1 Magnetization Evolution on spehre
clear; close all; 
run Plot_Magnetization_Evolution_EvoOnSphere

% Figure 2 Magnetization vs temperature
clear; close all; 
run Plot_Magnetization

% Figure 3 Susceptibility
clear; close all; 
run Plot_Susceptibility.m
% run CasiulisFig3

% Figure 4 Magnetization histogram
clear; close all; 
run Plot_MHistograms

% Figure 5 structural properties
clear; close all; 
run Plot_Structure

% Figure 6 spin correlation function
clear; close all; 
run Plot_SpinSCF

% Figure 7 a spin correlation function FS
clear; close all; 
run Plot_SpinSCF_FSCompare_SinglePlot

% Figure 7 b eta analysis
clear; close all; 
run Plot_EtaXi

% Figure 8 chimq Multiplot
clear; close all; 
run Plot_chimq_Multiplot

% % Figure 9 Lepri-Ruffo short-time behavior
% clear; close all; 
% run Plot_LepriRuffo.m

% Figure 9 & 10 Lepri-Ruffo short-time and long-time multiplot
clear; close all; 
run Plot_LepriRuffo_Multiplot.m

% Figure 11 Multiplot of Spin TCF properties
clear; close all; 
run Plot_TCF_FFT_Multiplot.m

% Figure 12 Nelson-Fisher law for high frequencies
clear; close all; 
run Plot_NF_HighFreq.m

% Figure 13 Finite-size scaling of peak height (Probably discarded)
clear; close all; 
run Plot_TCF_FFT_FS.m

% Figure 14 Spin wave propagation speed
clear; close all; 
run Plot_omega_spinwavedispersion.m

% Figure 15 Omega & gamma plots
clear; close all; 
run Plot_omega_gamma_vs_q.m
clear; close all; 
run Plot_omega_gamma_vs_T.m
clear; close all; 
run Plot_omega_gamma_FitCompare.m

% Finally, run a bash script to collect the resulting plots
system('./gather_paper_plots.sh');


