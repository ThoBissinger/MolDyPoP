clear all
close all
T_values = 1.30 : 0.02: 1.60;
sqrtN = 128;
% pathbase = '/data/scc/thobi/200213_XYModel/sqrtN_64/T_1.34/run_1/output';
path_base = sprintf('/data/scc/thobi/200213_XYModel/nh/sqrtN_%d/', sqrtN);

filename='sampling_output_samp.m';
for i = 1 : length(T_values)
    T_cur = T_values(i);
    path_suffix = sprintf('T_%0.2f/run_1/output', T_cur);
    path_cur = strcat(path_base, path_suffix);
    
    create_data_package(path_cur,filename);
    
    N(i,:) = N_cur;
    Tmax(i,:) = Tmax_cur;
    averaging_times(i,:) = averaging_times_cur;
    rbin(i,:) = rbin_cur;
    qbin(i,:) = qbin_cur;
    H(i,:) = H_cur;
    H_2(i,:) = H_2_cur;
    W(i,:) = W_cur;
    W_2(i,:) = W_2_cur;
    M(i,:,:) = M_cur;
    M_2(i,:) = M_2_cur;
    M_4(i,:) = M_4_cur;
    absM(i,:) = absM_cur;
    Theta(i,:) = Theta_cur;
    temperature(i,:) = temperature_cur;
    abs_vortices(i,:) = abs_vortices_cur;
    signed_vortices(i,:) = signed_vortices_cur;
    mxq(i,:) = mxq_cur;
    myq(i,:) = myq_cur;
%     eq(i,:) = eq_cur;
    wq(i,:) = wq_cur;
    chimxq(i,:) = chimxq_cur;
    chimyq(i,:) = chimyq_cur;
    chieq(i,:) = chieq_cur;
    chiwq(i,:) = chiwq_cur;
    SCF_Spin(i,:) = SCF_Spin_cur;
    SCF_anglediff(i,:) = SCF_anglediff_cur;
    TCF_times(i,:) = TCF_times_cur;
    ACF_Spin(i,:) = ACF_Spin_cur;
    ACF_anglediff(i,:) = ACF_anglediff_cur;
    gxx(i,:) = gxx_cur;
    gxy(i,:) = gxy_cur;
    gxw(i,:) = gxw_cur;
    gxe(i,:) = gxe_cur;
    gyy(i,:) = gyy_cur;
    gyw(i,:) = gyw_cur;
    gye(i,:) = gye_cur;
    gww(i,:) = gww_cur;
    gwe(i,:) = gwe_cur;
    gee(i,:) = gee_cur;
    H_av(i) = H_av_cur;
    H_2_av(i) = H_2_av_cur;
    H_var(i) = H_var_cur;
    W_av(i) = W_av_cur;
    W_2_av(i) = W_2_av_cur;
    W_var(i) = W_var_cur;
    M_av(i,:) = M_av_cur;
    M_2_av(i) = M_2_av_cur;
    M_4_av(i) = M_4_av_cur;
    M_var(i) = M_var_cur;
    Binder_cum(i) = Binder_cum_cur;
    absM_av(i) = absM_av_cur;
    Theta_av(i) = Theta_av_cur;
    temperature_av(i) = temperature_av_cur;
    abs_vortices_av(i) = abs_vortices_av_cur;
    signed_vortices_av(i) = signed_vortices_av_cur;
    mxq_av(i,:) = mxq_av_cur;
    myq_av(i,:) = myq_av_cur;
    eq_av(i,:) = eq_av_cur;
    wq_av(i,:) = wq_av_cur;
    chimxq_av(i,:) = chimxq_av_cur;
    chimyq_av(i,:) = chimyq_av_cur;
    chieq_av(i,:) = chieq_av_cur;
    chiwq_av(i,:) = chiwq_av_cur;
    SCF_Spin_av(i,:) = SCF_Spin_av_cur;
    SCF_anglediff_av(i,:) = SCF_anglediff_av_cur;

end

indexmat = create_indexmat(length(averaging_times_cur),length(qbin_cur));
oldFold = cd(path_base);
save 'summarized_data'
cd(oldFold);