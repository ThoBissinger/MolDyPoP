%% 0 Initialization
clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/SymmAntisymmFit
run ../initialization_script
saveswitch=1;
basedir=sprintf('%s/SymmAntisymmFit',fig_base);
i_model = 1;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    modelname="MXY";

    sqrtN_vals = [16 32 64 128 256];
    dir="/data/scc/thobi/211201_LongerTime/mxy_3.00";
    dir_more_q="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";

    T_vals = [.03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25];
    T_dirs = {"T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25"};
%     T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
%     T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"};
    
    sampfilename="samp_Dynamics";
    
    L_vals=[9.25,18.5,37,74,148];
%     q_select=[1,3,6,9,12];
    q_select=[1,3,6,9];
    q_thresh = 9;

elseif (i_model == 2)
    
    L_vals=[16,32,64,128,256];
    
    
elseif (i_model == 3)
    L_vals=[9.25,18.5,37,74,148];
    
elseif (i_model == 4)
    L_vals=[16, 32, 64, 128];
    
end
N_N = numel(sqrtN_vals);
N_T = numel(T_vals);
N_q = numel(q_select);
markersym='o';

fitfunc_symm="exp(-gamma*x/2)*(cos(omega_1*x) + gamma/2/omega_1*sin(omega_1*x))";
fitfunc_asymm="exp(-gamma*x/2)*cos(omega_1*x)";

%% 1 Data assembly
eta_vals = zeros(1,N_T);
for i_T = 1:N_T
    T_dir=T_dirs{i_T};
    absM_vec=zeros(1,N_N);
    for i_N = 1:N_N
        sqrtN = sqrtN_vals(i_N);
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_more_q,sqrtN,T_dir,sampfilename);
        load(curfile,"absM_av");
        absM_vec(i_N) = absM_av;
    end
    eta_fitob = fit_eta_Magnetization_FS(absM_vec(1:N_N),L_vals(1:N_N));
    eta_vals(i_T) = eta_fitob.eta;
end


n_period = 10;
weightexp = 1;

gamma_vals=zeros(N_N,N_T,N_q);
omega_1_vals=zeros(N_N,N_T,N_q);
omega_0_vals_om0fit=zeros(N_N,N_T,N_q);
gamma_vals_om0fit=zeros(N_N,N_T,N_q);

gamma_vals_symm=zeros(N_N,N_T,N_q);
omega_1_vals_symm=zeros(N_N,N_T,N_q);
gamma_vals_asymm=zeros(N_N,N_T,N_q);
omega_1_vals_asymm=zeros(N_N,N_T,N_q);


max_FFT_vals=zeros(N_N,N_T,N_q);
gamma_FFT_vals=zeros(N_N,N_T,N_q);
omega_2_FFT_vals=zeros(N_N,N_T,N_q);
% omega_0_vals=zeros(N_T,N_q);
q_vals_collect=zeros(N_N,N_T,N_q);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    fprintf("%d ",sqrtN);
    for i_T = 1:N_T
        T = T_vals(i_T);
        T_dir = T_dirs{i_T};
        fprintf("%.3f ",T);
    
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir,sqrtN,T_dir,sampfilename);
        if (~ isfile(curfile))
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_more_q,sqrtN,T_dir,sampfilename);
        end
        load(curfile,"averaging_times","gmperpmperp","qbin");
        t = averaging_times;
        n_t = numel(t);
        n_q = numel(qbin);
        for i_q = 1:N_q
            if q_select(i_q) > n_q - 1
                curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_more_q,sqrtN,T_dir,sampfilename);
                load(curfile,"averaging_times","gmperpmperp","qbin");
                t = averaging_times;
                n_t = numel(t);
                n_q = numel(qbin);
            end
            q_vals_collect(i_N,i_T,i_q) = qbin(q_select(i_q));
    
            cf=real(gmperpmperp(q_select(i_q):n_q:end));
            cf=cf/cf(1);
    
            c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_1');
            gamma_vals(i_N,i_T,i_q) = c(1);
            omega_1_vals(i_N,i_T,i_q) = c(2);

            startvec=abs([c(1),c(2)]);

            c = fit_DampedOscillator_RealSpace(t,cf,n_period,weightexp,'omega_0');
            gamma_vals_om0fit(i_N,i_T,i_q) = c(1);
            omega_0_vals_om0fit(i_N,i_T,i_q) = c(2);

            fitob = fit(t(:),cf(:)/cf(1),fitfunc_symm,...
                'StartPoint',startvec,...
                'Lower',[0,0],'Upper',10*startvec);
            gamma_vals_symm(i_N,i_T,i_q) = fitob.gamma;
            omega_1_vals_symm(i_N,i_T,i_q) = fitob.omega_1;

            fitob = fit(t(:),cf(:)/cf(1),fitfunc_asymm,...
                'StartPoint',startvec,...
                'Lower',[0,0],'Upper',10*startvec);
            gamma_vals_asymm(i_N,i_T,i_q) = fitob.gamma;
            omega_1_vals_asymm(i_N,i_T,i_q) = fitob.omega_1;
    
            [ft_vals,om_vals]=FT_correlation(t, cf, 0);
            ft_vals=real(ft_vals);
            [ft_max,i_ft_max]=max(ft_vals);
            omega_2_FFT_vals(i_N,i_T,i_q)=abs(om_vals(i_ft_max));
            max_FFT_vals(i_N,i_T,i_q)=ft_max;
    %         omega_0_vals(i_N,i_T,i_q) = sqrt(c(2)^2 - gamma_cur^2/4);
        end
        fprintf("\n");
    end
end
save(sprintf('%s/%s_dataset',basedir,curmodel));


return
omega_0_vals = sqrt(omega_1_vals.^2 + gamma_vals.^2/4);
gamma_FFT_vals_a = 2./max_FFT_vals;
gamma_FFT_vals_b = omega_0_vals.^2 .* max_FFT_vals / 2;

bad_fit_indices = find(omega_1_vals.^2 < gamma_vals.^2 / 4);
omega_0_vals_adjusted = omega_0_vals;
omega_0_vals_adjusted(bad_fit_indices) = 0;
omega_1_vals_adjusted = omega_1_vals;
omega_1_vals_adjusted(bad_fit_indices) = 0;
omega_2_vals_adjusted = omega_2_FFT_vals;
omega_2_vals_adjusted(bad_fit_indices) = 0;
omega_0_vals_om0fit_adjusted = omega_0_vals_om0fit;
omega_0_vals_om0fit_adjusted(bad_fit_indices) = 0;

gamma_a_vals=zeros(N_N,N_T);
gamma_sigma_vals=zeros(N_N,N_T);
omega_0_a_vals=zeros(N_N,N_T);
omega_0_sigma_vals=zeros(N_N,N_T);
omega_1_a_vals=zeros(N_N,N_T);
omega_1_sigma_vals=zeros(N_N,N_T);
c_spinwavespeed_omega_0_vals=zeros(N_N,N_T);
c_spinwavespeed_omega_0_om0fit_vals=zeros(N_N,N_T);
c_spinwavespeed_omega_1_vals=zeros(N_N,N_T);
c_spinwavespeed_omega_2_vals=zeros(N_N,N_T);
fitfunc_coeffs=@(a,sigma,x) a*x.^sigma;
fitfunc_spinwave=@(c,x) c*x;
% ind_select=1:numel(q_vals_collect);

for i_N = 1:N_N
    for i_T = 1:N_T
        q_vec=q_vals_collect(i_N,i_T,:);
        ind_select=1:numel(q_vec);
        vec_cur=gamma_vals(i_N,i_T,ind_select);
        init_sigma=2;
        init_a=vec_cur(end)/q_vec(end)^init_sigma;
        fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
            'StartPoint',[init_a,init_sigma]);
        gamma_a_vals(i_N,i_T) = fitob.a;
        gamma_sigma_vals(i_N,i_T) = fitob.sigma;
    
        vec_cur=omega_1_vals(i_N,i_T,ind_select);
        init_sigma=1;
        init_a=vec_cur(end)/q_vec(end)^init_sigma;
        fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
            'StartPoint',[init_a,init_sigma]);
        omega_1_a_vals(i_N,i_T) = fitob.a;
        omega_1_sigma_vals(i_N,i_T) = fitob.sigma;
    
        vec_cur=omega_0_vals_om0fit_adjusted(i_N,i_T,ind_select);
        init_c = .2;
        fitob=fit(q_vec(:),vec_cur(:),fittype(fitfunc_spinwave),...
            'StartPoint',[init_c]);
        c_spinwavespeed_omega_0_om0fit_vals(i_N,i_T) = fitob.c;
    
        vec_cur=omega_0_vals_adjusted(i_N,i_T,ind_select);
        init_c = .2;
        fitob=fit(q_vec(:),vec_cur(:),fittype(fitfunc_spinwave),...
            'StartPoint',[init_c]);
        c_spinwavespeed_omega_0_vals(i_N,i_T) = fitob.c;

        vec_cur=omega_1_vals_adjusted(i_N,i_T,ind_select);
        init_c = .2;
        fitob=fit(q_vec(:),vec_cur(:),fittype(fitfunc_spinwave),...
            'StartPoint',[init_c]);
        c_spinwavespeed_omega_1_vals(i_N,i_T) = fitob.c;
    
        vec_cur=omega_2_vals_adjusted(i_N,i_T,ind_select);
        init_c = .2;
        fitob=fit(q_vec(:),vec_cur(:),fittype(fitfunc_spinwave),...
            'StartPoint',[init_c]);
        c_spinwavespeed_omega_2_vals(i_N,i_T) = fitob.c;

        vec_cur=omega_0_vals(i_N,i_T,ind_select);
        init_sigma=1;
        init_a=vec_cur(end)/q_vec(end)^init_sigma;
        fitob = fit(q_vec(:),vec_cur(:),fittype(fitfunc_coeffs),...
            'StartPoint',[init_a,init_sigma]);
        omega_0_a_vals(i_N,i_T) = fitob.a;
        omega_0_sigma_vals(i_N,i_T) = fitob.sigma;
    end
end

save(sprintf('%s/%s_dataset',basedir,curmodel));








