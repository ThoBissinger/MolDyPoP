%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/correlation_times',fig_base);

%% 0 Preparing data
model="mxy";modeltag="\textbf{MXY}";modeldir="mxy_3.00"; T = .11; T_dir = "T_.11"; simdir='210715_LinearTimeSampling'; simdir='220201_ReducedSmapleStepDeltat';
% model="fmxy";modeltag="\textbf{FMXY}";modeldir="fmxy"; T = .17; T_dir = "T_.17"; simdir='210715_LinearTimeSampling';
% model="xy";modeltag="\textbf{XY}";modeldir="xy_s"; T = .89; T_dir = "T_.89"; simdir='211201_LongerTime';

% T_dirs=["T_.18","T_.185", "T_.19","T_.195","T_.20","T_.205","T_.21","T_.22","T_.23","T_.24","T_.25","T_.27","T_.29","T_.31","T_.33","T_.35","T_.37","T_.40","T_.43","T_.46","T_.49","T_.52"];
% T_vals=[.18,.185, .19,.195,.20,.205,.21,.22,.23,.24,.25,.27,.29,.31,.33,.35,.37,.40,.43,.46,.49,.52];

T_dirs=["T_.18","T_.185", "T_.19","T_.195","T_.20","T_.205","T_.21","T_.22","T_.23","T_.24"];
T_vals=[.18,.185, .19,.195,.20,.205,.21,.22,.23,.24];

% sqrtN_vals = [16 32 64 128];
% L_vals = [9.25 18.5 37 74];
sqrtN_vals=64;
L_vals=37;

% last_q_switch = 1; % If 0, all q values can be used, skips the last entry if 1
n_T=numel(T_vals);
n_N=numel(sqrtN_vals);
t_cell=cell(n_N,n_T);
r_cell=cell(n_N,n_T);
q_cell=cell(n_N,n_T);
tcf_cell=cell(n_N,n_T);
scf_cell=cell(n_N,n_T);
scf_q_cell=cell(n_N,n_T);
absM_vals=zeros(n_N,n_T);
for i_N = 1:n_N
    sqrtN=sqrtN_vals(i_N);
    for i_T = 1:n_T
        T = T_vals(i_T);
        T_dir = T_dirs(i_T);
        file=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics.mat',...
            simdir,modeldir,sqrtN,T_dir);
        S=load(file,'gmperpmperp','averaging_times','qbin',...
            'chimperpq_av','SCF_Spin_av','rbin',...
            'absM_av');
        t_cell{i_N,i_T}=S.averaging_times;
        r_cell{i_N,i_T}=S.rbin;
        q_cell{i_N,i_T}=S.qbin;
        tcf_cell{i_N,i_T}=real(S.gmperpmperp)/sqrtN^2;
        scf_cell{i_N,i_T}=S.SCF_Spin_av;
        scf_q_cell{i_N,i_T}=real(S.chimperpq_av)/sqrtN^2;
        absM_vals(i_N,i_T)=S.absM_av;
    end
end

% sqrtN=128; L=74;
runmax=500;

dt_sim = .01;
% 211201_LongerTime
% 210715_LinearTimeSampling
res_factor=3/4*log(runmax);
res_function = resolution_Gauss;
% res_function=@(t,tau) resolution_Laplace_pleateau(t,tau,tau);

% fitfunc="exp(-gamma*x/2)*(cos(omega*x)+gamma/omega/2*sin(omega*x))";

%% Fitting
fitfunc="a*exp(-x/b)";
fitfunc_offset="a*exp(-x/b)+c";
fitfunc_double="a*exp(-x/tau_1)+b*exp(-x/tau_2)";

tau_cell=cell(n_N,n_T);
xi_vals=zeros(n_N,n_T);
xi_a_vals=zeros(n_N,n_T);
xi_b_vals=zeros(n_N,n_T);
q_non0_cell=cell(n_N,n_T);
% tau_1_cell=cell(n_N,n_T);
% tau_2_cell=cell(n_N,n_T);
fitob_cell=cell(n_N,n_T);
fitob_xi_cell=cell(n_N,n_T);
for i_N = 1:n_N
    sqrtN=sqrtN_vals(i_N);
    for i_T = 1:n_T
        T = T_vals(i_T);
        T_dir = T_dirs(i_T);
        fprintf('Fit -- sqrtN = %d -- T = %.3f\n',sqrtN,T);
        q_vals=q_cell{i_N,i_T};
        t=t_cell{i_N,i_T};
        cf_full=tcf_cell{i_N,i_T};
        n_q=numel(q_vals);
        nonzero_ind=find(cf_full(1:n_q));
        q_vals_non0 = q_vals(nonzero_ind);
        n_q_non0 = numel(q_vals_non0);

        tau_vals=zeros(1,n_q_non0);
%         tau_1_vals=zeros(1,n_q - last_q_switch);
%         tau_2_vals=zeros(1,n_q - last_q_switch);
        fitob_sub_cell=cell(1,n_q_non0);
        for i_q_non0 = 1:n_q_non0
            i_q = nonzero_ind(i_q_non0);
            cf=cf_full(i_q:n_q:end);
            fitob=fit(t(:),cf(:)/cf(1),fitfunc,'StartPoint',[1,1]);
%             fitob_double=fit(t(:),cf(:)/cf(1),fitfunc_double,'StartPoint',[fitob.a,fitob.tau,1,1]);
            fitob_sub_cell{i_q_non0}=fitob;
            tau_vals(i_q_non0)=fitob.b;
%             tau_1_vals(i_q)=fitob_double.tau_1;
%             tau_2_vals(i_q)=fitob_double.tau_2;
        end
        tau_cell{i_N,i_T}=tau_vals;
        q_non0_cell{i_N,i_T}=q_vals_non0;
%         tau_1_cell{i_N,i_T}=tau_1_vals;
%         tau_2_cell{i_N,i_T}=tau_2_vals;
        fitob_cell{i_N,i_T}=fitob_sub_cell;
    
        r_vals=r_cell{i_N,i_T};
        scf=scf_cell{i_N,i_T};
%         fitob=fit(r_vals(:),scf(:),fitfunc_offset,'StartPoint',[1,1,absM_vals(i_N,i_T)^2]);
        fitob=fit(r_vals(:),scf(:),fitfunc,'StartPoint',[1,1]);
        xi_vals(i_N,i_T)=fitob.b;
        xi_a_vals(i_N,i_T)=fitob.a;
%         xi_c_vals(i_N,i_T)=fitob.c;
        fitob_xi_cell{i_N,i_T}=fitob;
    end
end
%% Plotting
i_N=1;
i_T=10;
z=1;
fprintf("%.3f\n",T_vals(i_T));
q_vals=q_cell{i_N,i_T};
n_q=numel(q_vals);
t=t_cell{i_N,i_T}; 
n_t=numel(t);
cf_full=tcf_cell{i_N,i_T};
eta=.25;
for i_q=1:10
    xi=xi_vals(i_N,i_T);
    q=q_vals(i_q);
    cf=cf_full(i_q:n_q:end);
    if(max(abs(cf))~=0)
        dispname = sprintf('$q=%.3f$',q);
        plot(t*q^(z),q^((z+2-eta))*cf,'DisplayName',dispname);
        hold on;
    end
end
xlim([0 200])
hLegend=legend('Location','northeast','Interpreter','latex');
hold off;
