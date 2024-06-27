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
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    for i_T = 1:N_T
        T = T_str(i_T);
        curmatfile=datadir + "/sqrtN_" + sqrtN + "/T_" + T + "/" + matfilename;
        S=load(curmatfile,'absM_av');
        absM_vals(i_N,i_T) = S.absM_av;
    end
end

%% Creating the Figure
figure
c_map = lines(N_N);
for i_N = 1:N_N
    plot(T_vals,absM_vals(i_N,:),...
        'DisplayName',sprintf("$N = (%d)^2$",sqrtN_vals(i_N)),...
        'Color',c_map(i_N,:),...
        'Marker','o');
    hold on;
end
legend show;
legend('interpreter','latex','FontSize',14);
xlabel('$T$',"interpreter","latex",'FontSize',14);
ylabel('$|M|$',"interpreter","latex",'FontSize',14);
set(gca,'FontSize',12)

exportgraphics(gcf,"Magnetization.pdf")