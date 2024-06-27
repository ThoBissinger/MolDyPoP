%% 0 Initialization
clear
run initialization_script;
saveswitch=1;
basedir=sprintf('%s/plots/static/magnetization_fits',fig_base);


% model="mxy";modeltag="\textbf{MXY}";modeldir="mxy_3.00"; simdir='210715_LinearTimeSampling'; % simdir='220201_ReducedSmapleStepDeltat';
% model="fmxy";modeltag="\textbf{FMXY}";modeldir="fmxy"; T = .17; T_dir = "T_.17"; simdir='210715_LinearTimeSampling';
model="xy";modeltag="\textbf{XY}";modeldir="xy_s"; T = .89; T_dir = "T_.89"; simdir='201207_equilibration/xy_s/scale';
sqrtN_vals = [16 32 64 128 256];
L_vals = [9.25 18.5 37 74 148];
T_vals = [.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
T_dirs = ["T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40" "T_.43" "T_.46" "T_.49" "T_.52"];

if model == "xy"
    sqrtN_vals = [16 32 64 128 256];
    L_vals = sqrtN_vals;
    T_vals=[.10, .20, .30, .40, .50, .60, .70, .80, .85, .87, .89, .90, .91, .93, .95, .97, 1.00, 1.03, 1.06, 1.09, 1.10, 1.12, 1.15, 1.18, 1.20, 1.21, 1.24, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00]; 
   T_dirs=["T_.10", "T_.20", "T_.30", "T_.40", "T_.50", "T_.60", "T_.70", "T_.80", "T_.85", "T_.87", "T_.89", "T_.90", "T_.91", "T_.93", "T_.95", "T_.97", "T_1.00", "T_1.03", "T_1.06", "T_1.09", "T_1.10", "T_1.12", "T_1.15", "T_1.18", "T_1.20", "T_1.21", "T_1.24", "T_1.30", "T_1.40", "T_1.50", "T_1.60", "T_1.70", "T_1.80", "T_1.90", "T_2.00"]; 
end


n_T=numel(T_vals);
n_N=numel(sqrtN_vals);

alpha_vals=[sqrt(2),1.3585,1,.5];
alpha_names=["$\alpha = \sqrt{2}$", "$\alpha = 1.3585$", "$\alpha = 1$", "$\alpha = 1/2$"];
alpha_symbols=["^","v","d","s"];
n_alpha=numel(alpha_vals);

absM_vals=zeros(n_N,n_T);
for i_N = 1:n_N
    sqrtN=sqrtN_vals(i_N);
    for i_T = 1:n_T
        T = T_vals(i_T);
        T_dir = T_dirs(i_T);
        file=sprintf('/data/scc/thobi/%s/%s/sqrtN_%d/%s/samp_Dynamics.mat',...
            simdir,modeldir,sqrtN,T_dir);
        if model == "xy"
            file=sprintf('/data/scc/thobi/%s/sqrtN_%d/%s/samp_eq.mat',...
                simdir,sqrtN,T_dir);
            if sqrtN == 16
                file=sprintf('/data/scc/thobi/201207_equilibration/xy_s/anneal/sqrtN_16/%s/samp_eq.mat',...
                    T_dir);
            end
        end
        S=load(file,'absM_av');
        absM_vals(i_N,i_T)=S.absM_av;
    end
end

%% Fitting
clear T_C T_star T_KT crossover_M mag_fitcuns;
T_KT_cell=cell(1,n_alpha);
T_star_cell=cell(1,n_alpha);
T_C_cell=cell(1,n_alpha);
cross_M_cell=cell(1,n_alpha);
mag_fitfuncs_cell=cell(1,n_alpha);
for i_alpha=1:n_alpha
    alpha=alpha_vals(i_alpha);
    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_vals,sqrtN_vals,alpha);
    T_KT_cell{i_alpha}=T_KT;
    T_star_cell{i_alpha}=T_star;
    T_C_cell{i_alpha}=T_C;
    cross_M_cell{i_alpha}=crossover_M;
    mag_fitfuncs_cell{i_alpha}=mag_fitfuncs;
end

clear T_C T_star T_KT crossover_M mag_fitcuns;
n_points=300;
T_KT_longvec=zeros(1,n_points);
T_star_longvec=zeros(1,n_points);
T_C_longvec=zeros(1,n_points);
T_C_fit_1=zeros(1,n_points);
T_C_fit_2=zeros(1,n_points);
T_C_fit_3=zeros(1,n_points);
T_C_fit_spline=zeros(1,n_points);
aa = linspace(.1,3,300);
beta_exp=3*pi^2/128;
for i_alpha = 1:n_points
    alpha = aa(i_alpha);
    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_vals(end,:),sqrtN_vals(end),alpha);
    T_KT_longvec(i_alpha)=T_KT;
    T_star_longvec(i_alpha)=T_star;
    T_C_longvec(i_alpha)=T_C;

    fitfunc_TC="a*(b-x)^.2313";
    i_nearest=find(T_vals<T_star,1,'last');
    startpoint=[crossover_M/(T_C-T_star)^beta_exp,T_C];

    try
        ind=i_nearest:i_nearest+1;
        T_cur=T_vals(ind);
        M_cur=absM_vals(end,ind);
        fitob=fit(T_cur(:),M_cur(:),fitfunc_TC,'StartPoint',startpoint);
        T_C_fit_1(i_alpha)=fitob.b;
    catch err
        T_C_fit_1(i_alpha)=0;
    end

    try
        ind=i_nearest-1:i_nearest+2;
        T_cur=T_vals(ind);
        M_cur=absM_vals(end,ind);
        fitob=fit(T_cur(:),M_cur(:),fitfunc_TC,'StartPoint',startpoint);
        T_C_fit_2(i_alpha)=fitob.b;
    catch err
        T_C_fit_2(i_alpha)=0;
    end

    try
        ind=i_nearest-2:i_nearest+3;
        T_cur=T_vals(ind);
        M_cur=absM_vals(end,ind);
        fitob=fit(T_cur(:),M_cur(:),fitfunc_TC,'StartPoint',startpoint);
        T_C_fit_3(i_alpha)=fitob.b;
    catch err
        T_C_fit_3(i_alpha)=0;
    end
    
    M_spline=spline(T_vals,absM_vals(end,:));
    diffT=T_C-T_star;
    xx=linspace(T_star-2*diffT,T_star+.75*diffT);
    yy=fnval(M_spline,xx);
    fitob=fit(xx(:),yy(:),fitfunc_TC,'StartPoint',startpoint);
    T_C_fit_spline(i_alpha)=fitob.b;

end

%% Plotting: all in one
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=linspecer(3);
for i_alpha=[1,2]
    alpha=alpha_vals(i_alpha);
    semilogx(sqrtN_vals.^2,T_C_cell{i_alpha},...
        'Color',c_map(1,:),'Marker',alpha_symbols(i_alpha),...
        "HandleVisibility","off");
    hold on;
    semilogx(sqrtN_vals.^2,T_star_cell{i_alpha},...
        'Color',c_map(2,:),'Marker',alpha_symbols(i_alpha),...
        "HandleVisibility","off");
    semilogx(sqrtN_vals.^2,T_KT_cell{i_alpha},...
        'Color',c_map(3,:),'Marker',alpha_symbols(i_alpha),...
        "HandleVisibility","off");
end
semilogx(NaN,NaN,'Color',c_map(1,:),"DisplayName","$T_{C}$");
semilogx(NaN,NaN,'Color',c_map(2,:),"DisplayName","$T^*$");
semilogx(NaN,NaN,'Color',c_map(3,:),"DisplayName","$T_{BKT}$");

hlegend=legend('Location','northeast','Interpreter','latex');
hXLabel = xlabel('$N$','interpreter','latex',...
    'FontName', 'cmr12');
hYLabel = ylabel('$T$','interpreter','latex',...
    'FontName', 'cmr12');

figname=sprintf('%s/%s_AllInOne',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% Plotting: triptychon
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 pagewidth_cm columnwidth_cm]);
T_cells={T_C_cell,T_star_cell,T_KT_cell};
T_names=["$T_{C}$","$T^*$","$T_{BKT}$"];
for i_T=1:3
    subplot(1,3,i_T)
    c_map=linspecer(n_alpha);
    for i_alpha=1:n_alpha
        alpha=alpha_vals(i_alpha);
        T_cell=T_cells{i_T};
        T_trans_vals=T_cell{i_alpha};
        semilogx(sqrtN_vals.^2,T_trans_vals,...
            'Color',c_map(i_alpha,:),'Marker',alpha_symbols(i_alpha),...
            "DisplayName",alpha_names(i_alpha),'linewidth',1.0);
        hold on;
    end
    ylim([.1 .25]);
    if model == "xy"
        ylim([.6 1.25]);
    end
    if i_T == 1
        hlegend=legend('Location','southwest','Interpreter','latex');
    end
    htitle=title(T_names(i_T),'Interpreter','latex');
    hXLabel = xlabel('$N$','interpreter','latex',...
        'FontName', 'cmr12');
    hYLabel = ylabel('$T$','interpreter','latex',...
        'FontName', 'cmr12');
end

figname=sprintf('%s/%s_Triptychon',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% Plotting: alphacompare 
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=linspecer(3);
plot(aa,T_C_longvec,...
    'Color',c_map(1,:),'linewidth',2,...
    "DisplayName","$T_{C}(N_{\max})$");
hold on;
plot(aa,T_star_longvec,...
    'Color',c_map(2,:),'linewidth',2,...
    "DisplayName","$T^*(N_{\max})$");
plot(aa,T_KT_longvec,...
    'Color',c_map(3,:),'linewidth',2,...
    "DisplayName","$T_{BKT}(N_{\max})$");
hlegend=legend('Location','southeast','Interpreter','latex');
hXLabel = xlabel('$\alpha$','interpreter','latex',...
    'FontName', 'cmr12');
hYLabel = ylabel('$T$','interpreter','latex',...
    'FontName', 'cmr12');
ylim([.1 .2]);
if model == "xy"
    ylim([.6 1.00]);
end
set(gca,'XGrid','on','YGrid','on');

figname=sprintf('%s/%s_alphacompare',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end


%% Plotting: T_C_fitcompare
figure
set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
c_map=linspecer(5);
plot(aa,T_C_longvec,...
    'Color',c_map(1,:),'linewidth',2,...
    "DisplayName","$T_{C}$ standard fit");
hold on;
plot(aa,T_C_fit_1,...
    'Color',c_map(2,:),'linewidth',2,...
    "DisplayName","$T_{C}$ fit 1");
plot(aa,T_C_fit_2,...
    'Color',c_map(3,:),'linewidth',2,...
    "DisplayName","$T_{C}$ fit 2");
plot(aa,T_C_fit_3,...
    'Color',c_map(4,:),'linewidth',2,...
    "DisplayName","$T_{C}$ fit 3");
plot(aa,T_C_fit_spline,...
    'Color',c_map(5,:),'linewidth',2,...
    "DisplayName","$T_{C}$ fit (spline)");
hlegend=legend('Location','southeast','Interpreter','latex');
hXLabel = xlabel('$\alpha$','interpreter','latex',...
    'FontName', 'cmr12');
hYLabel = ylabel('$T$','interpreter','latex',...
    'FontName', 'cmr12');
ylim([.17 .2]);
if model == "xy"
    ylim([.9 1.00]);
end
set(gca,'XGrid','on','YGrid','on');

figname=sprintf('%s/%s_T_C_fitcompare',basedir,model);
if(saveswitch == 1)
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end
