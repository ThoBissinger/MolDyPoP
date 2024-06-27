%% Setting the stage
basedir="/data/scc/thobi/00_MolDyPop_Final/Example/Matlab/Plots"; % CHANGE TO SETUP
datadir="/data/scc/thobi/00_MolDyPop_Final/Example/integ_data/mxy";
matfilename="samp_integ.mat";

sqrtN_vals=[16,32];
L_vals=[9.25,18.5];
T_str = [".05" ".09" ".13" ".15" ".17" ".19" ".21" ".23" ".25" ];
% T_str = [".05" ".13" ".15" ".17" ".19" ".21" ".23" ".25" ];
T_vals=str2double(T_str);
N_N=numel(sqrtN_vals);
N_T=numel(T_vals);

%% Collecting the values
absM_vals=zeros(N_N,N_T);
qbin_cell=cell(N_N,N_T);
t_cell=cell(N_N,N_T);
TCF_cell=cell(N_N,N_T);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    for i_T = 1:N_T
        T = T_str(i_T);
        curmatfile=datadir + "/sqrtN_" + sqrtN + "/T_" + T + "/" + matfilename;
        S=load(curmatfile,'absM_av',"gmperpmperp","qbin","averaging_times");
        absM_vals(i_N,i_T) = S.absM_av;
        qbin_cell{i_N,i_T} = S.qbin;
        t_cell{i_N,i_T} = S.averaging_times;
        TCF_cell{i_N,i_T} = S.gmperpmperp;
    end
end

%% Creating the Figure
figure
c_map = lines(N_N);

TCF_cur=TCF_cell{1,2};
q_cur=qbin_cell{1,2};
t_cur=t_cell{1,2};
TCF_cur=real(TCF_cur(1,1:numel(q_cur):end));
plot(t_cur,TCF_cur/TCF_cur(1),...
    'DisplayName',"$N = (16)^2$",...
    'Color','red');
hold on;
TCF_cur=TCF_cell{2,2};
q_cur=qbin_cell{2,2};
t_cur=t_cell{2,2};
TCF_cur=real(TCF_cur(1,2:numel(q_cur):end));
plot(t_cur,TCF_cur/TCF_cur(1),...
    'DisplayName',"$N = (32)^2$",...
    'Color','blue');

legend show;
legend('interpreter','latex','FontSize',14);
xlabel('$t$',"interpreter","latex",'FontSize',14);
ylabel('$C_{m\perp}(q,t)/\chi_{m\perp}(q)$',"interpreter","latex",'FontSize',14);
set(gca,'FontSize',12)
xlim([0,1e2]);

exportgraphics(gcf,"TCF.pdf")