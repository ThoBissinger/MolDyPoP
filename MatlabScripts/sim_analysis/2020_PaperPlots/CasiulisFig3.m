%% 0 Initialization
run initialization_script;
basedir=sprintf('%s/plots/static',fig_base);

saveswitch=1;
i_model = 1;
if (i_model == 1)
    curmodel="mxy";
    curtitle="MXY model";
    dir='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
    sqrtN_vals=[16 32 64 128 256];
    T_dirs = {"T_.01" "T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25" "T_.27" "T_.29" "T_.31" "T_.33" "T_.35" "T_.37" "T_.40"};
    T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40];

    collectfilename='samp_Dynamics';

%     data=load(mxydata_name);
    runmax=500;

    L_vals=[9.25,18.5,37,74,148];
elseif (i_model == 2)
    curmodel="xy";
    curtitle="SXY model";
    collectfilename='samp_Dynamics.mat';
    
    data=load(xydata_name);
    runmax=250;

    L_vals=[16,32,64,128,256];
    T_max = 3.5;
elseif (i_model == 3)
    curmodel="fmxy";
    curtitle="FMXY model"; 
    collectfilename='samp_Dynamics_collect.mat';
    
    data=load(fmxydata_name);
    runmax=250;
    L_vals=[9.25,18.5,37,74,148];
    T_max = .4;
elseif (i_model == 4)
    curmodel="xy_s";
    curtitle="SXY model";
    runmax=125;
    collectfilename='samp_eq_collect.mat';
    
    data=matfile(xysdata_name);
    FS_Mag_fitdata=matfile(xysfit_FSMag_name);

    L_vals=[16 32 64 128 256];
    T_max=2;
end
N_N = numel(sqrtN_vals);
N_T = numel(T_dirs);



%% 1 Data Assembly
absM_mean=zeros(N_N,N_T);
absM_err=zeros(N_N,N_T);
M_2_mean=zeros(N_N,N_T);
M_2_err=zeros(N_N,N_T);
M_4_mean=zeros(N_N,N_T);
M_4_err=zeros(N_N,N_T);
chi_mean=zeros(N_N,N_T);
chi_err=zeros(N_N,N_T);
for i_N = 1:N_N
    sqrtN = sqrtN_vals(i_N);
    N = sqrtN^2;
    fprintf('sqrtN = %d',sqrtN);
    for i_T = 1:N_T
        T = T_vals(i_T);
        fprintf(' %.3g',T);
        curfile=sprintf('%s/sqrtN_%d/%s/%s_collect.mat',dir,sqrtN,T_dirs{i_T},collectfilename);
        load(curfile,'absM_av_collect','M_2_av_collect','M_4_av_collect')
        absM_mean(i_N,i_T)=mean(absM_av_collect);
        absM_err(i_N,i_T)=std(absM_av_collect) / sqrt(runmax);
        M_2_mean(i_N,i_T)=mean(M_2_av_collect);
        M_2_err(i_N,i_T)=std(M_2_av_collect) / sqrt(runmax);
        M_4_mean(i_N,i_T)=mean(M_4_av_collect);
        M_4_err(i_N,i_T)=std(M_4_av_collect) / sqrt(runmax);
        chicur=N/T*(M_2_av_collect - absM_av_collect.^2);
        chi_mean(i_N,i_T)=mean(chicur);
        chi_err(i_N,i_T)=std(chicur) / sqrt(runmax);
    end
    fprintf(' \n');
end
%% Fit
eta_fitobs=cell(1,N_T);
eta_vals=0*T_vals;
for i_T = 1:N_T
    absM_vec=absM_mean(:,i_T);
    eta_fitobs{i_T}=fit_eta_Magnetization_FS(absM_vec,sqrtN_vals);
    eta_vals(i_T)=eta_fitobs{i_T}.eta;
end




%% 2 The plot
figure
hold on
c_map = linspecer(N_N);
N_vals=sqrtN_vals.^2;
[max_chiabsm,i_max_vec]=max((T_vals.*chi_mean)');
chi_err_vec=T_vals(i_max_vec).*chi_err(i_max_vec)*sqrt(runmax);
% i_max_vec=zeros(size(sqrtN_vals));
% chiabsm=cell(size(sqrtN_vals));
% chierr=cell(size(sqrtN_vals));
for i_N = 1 : N_N
    sqrtN=sqrtN_vals(i_N);
    N=sqrtN^2;
    chiabsm_vec=chi_mean(i_N,:);
%         chiabsm_vec=zeros(size(T_vals));
%         chiabsm_vec=N_vals(i_N)*(M_2_av(i_N,:) - absM_av(i_N,:).^2) ./T_vals;
%     var_vec = M_2_av(i_N,:) - absM_av(i_N,:).^2;
%     chiabsm_vec = N_vals(i_N) * var_vec ./ T_vals;
%     chiabsm{i_N} = chiabsm_vec;

%     varerr_vec=sqrt(1/runmax * (M_2_av(i_N,:) - (runmax - 3) / (runmax - 1) *var_vec.^2));
    % Apparently pg 483, Rao, 1973
%     chierr{i_N} = N_vals(i_N) * varerr_vec ./ T_vals;

%     [max_chiabsm(i_N),i_max_vec(i_N)] = max(T_vals.*chiabsm_vec);
%     fitfunc=mag_fitfuncs{i_N};

%         hchiabsM_line(i_N) = line(T_vals,chiabsm_vec);
    exponents=2*(1-eta_vals);
    hchiabsM_line(i_N) = line(T_vals,chiabsm_vec./sqrtN_vals(i_N).^exponents);
    set(hchiabsM_line(i_N), ...
        'LineStyle', '-', 'LineWidth', 1, ...
        'HandleVisibility','off',...
        'Color', c_map(i_N,:))

    
%         hchiabsM_dots(i_N) = line(T_vals,chiabsm_vec);
    dispname = sprintf('$N = (%d)^2$',sqrtN);
    hchiabsM_dots(i_N) = line(T_vals,chiabsm_vec./sqrtN_vals(i_N).^exponents);
    set(hchiabsM_dots(i_N), ...
        'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', 5, ...
        'HandleVisibility', 'off', ...
        'Color',c_map(i_N,:),...
        'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:));
    
    hlabel(i_N)=plot(NaN,NaN,...
        'LineStyle', '-', ...
        'LineWidth', 1,...
        'Marker', 'o', 'MarkerSize', 5, ...
        'DisplayName', dispname, ...
        'Color',c_map(i_N,:),...
        'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:));
    
end

hLegend = legend('Location', 'NorthWest','interpreter','latex');

hXLabel = xlabel('$T$','interpreter','latex');
hYLabel = ylabel('$\chi_{m}$','interpreter','latex');
h_axis = gca;
xlim([0, max(T_vals)]);

% 3a Font
set_fonts_default


% 3a Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', .5)
ax_full = gca;

% 3a Inset
ax_inset = axes('Position',[.62 .62 .26 .26]);



hchimax_line = errorbar(N_vals,max_chiabsm./N_vals,chi_err_vec./N_vals);
f_POW = fittype('a*(x)^(2-eta)');
init_exponent = .25;
init_factor = 1;
c0=[init_factor,init_exponent];
Lower_lim=[0, 0];
Upper_lim=[Inf, Inf];

%     c = fit(sqrtN_vals(1:end-1)', max_chiabsm(1:end-1)', f_POW,'StartPoint',c0,...
%         'Lower',Lower_lim, 'Upper', Upper_lim);
c = fit(sqrtN_vals(1:end)', max_chiabsm(1:end)', f_POW,'StartPoint',c0,...
    'Lower',Lower_lim, 'Upper', Upper_lim,...
    'Algorithm','Trust-Region', 'Robust', 'LAR');

f_POW_fixedeta = fittype('a*(x)^(2-.25)');
init_factor = 1;
c0=[init_factor];
Lower_lim=[0];
Upper_lim=[Inf];

c_fixedeta = fit(sqrtN_vals(:), max_chiabsm(:), f_POW_fixedeta,'StartPoint',c0,...
    'Lower',Lower_lim, 'Upper', Upper_lim);
NN=logspace(1,6);

hchimax_fit = line(NN,c_fixedeta.a * sqrt(NN).^(- .25));
%     hchimax_fit = line(NN,c.a * sqrt(NN).^(- c.eta));

c_map = linspecer(3);
set(hchimax_line, ...
    'LineStyle', 'none', 'LineWidth', 1, ...
    'Marker', '^', 'MarkerSize', 5, ...
    'Color',c_map(1,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(1,:))

set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'XTick', [1e1 1e2 1e3 1e4], 'YTick', 0:1e-2:5e-2, ...
    'LineWidth', .5, 'Xscale', 'log','Yscale','lin')

hXLabel = xlabel('$N$','interpreter','latex');
hYLabel = ylabel('$\chi_{\textrm{max}} / N$','interpreter','latex');

set_fonts_default


xlim([.5*N_vals(1),2*N_vals(end)])
ylim([0 .05])
%     ylim([.006,.025])

 % 3a Text
%      annotation('textbox','String',sprintf('$\\eta = %.2f$',c.eta),...
%      annotation('textbox','String',sprintf('$\\eta = %.2f$',.25),...
 annotation('textbox','String',sprintf('$\\eta = %.2f$',.25),...
     'Position',ax_inset.Position + [.07 .155 .0 .0],...
     'LineWidth',.5,...
     'HorizontalAlignment','left','Vert','bottom',...
     'interpreter','latex','FitBoxToText','on',...
     'BackgroundColor','White',...
     'FontSize', fontsize_annotation, 'FontName', 'cmr12')

 set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);

figname=sprintf('%s/%s_AbsSusceptibility',basedir,curmodel);
fprintf('Creating figure %s\n',figname)
if(saveswitch == 1)
   exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
   exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end





















%% Old implementation
run initialization_script;
basedir=sprintf('%s/plots/static',fig_base);

saveswitch=1;

i_model = 1;
% for i_model = [1,3,4]

    
    if (i_model == 1)
        curmodel="mxy";
        curtitle="MXY model";
        collectfilename='samp_Dynamics_collect.mat';

        data=load(mxydata_name);
        runmax=500;

        L_vals=[9.25,18.5,37,74,148];
        T_max = .4;
    elseif (i_model == 2)
        curmodel="xy";
        curtitle="SXY model";
        collectfilename='samp_Dynamics.mat';
        
        data=load(xydata_name);
        runmax=250;

        L_vals=[16,32,64,128,256];
        T_max = 3.5;
    elseif (i_model == 3)
        curmodel="fmxy";
        curtitle="FMXY model"; 
        collectfilename='samp_Dynamics_collect.mat';
        
        data=load(fmxydata_name);
        runmax=250;
        L_vals=[9.25,18.5,37,74,148];
        T_max = .4;
    elseif (i_model == 4)
        curmodel="xy_s";
        curtitle="SXY model";
        runmax=125;
        collectfilename='samp_eq_collect.mat';
        
        data=matfile(xysdata_name);
        FS_Mag_fitdata=matfile(xysfit_FSMag_name);

        L_vals=[16 32 64 128 256];
        T_max=2;
    end
    T_vals=data.('T_vals');
    sqrtN_vals=data.('sqrtN_vals');
    N_vals=sqrtN_vals.^2;
    absM_av=cell2mat(data.('absM_av'));
    M_2_av=cell2mat(data.('M_2_av'));
    M_4_av=cell2mat(data.('M_4_av'));
    
    T_dirs = data.('T_dirs');
    collectfile_dir = data.('dirs');
    sqrtN_dirs = data.('sqrtN_dirs');
    
    labels=["$N = (16)^{2}$", "$N = (32)^{2}$", "$N = (64)^{2}$", "$N = (128)^{2}$", "$N = (256)^{2}$"];
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    [T_C,T_star,T_KT,crossover_M,mag_fitfuncs] = fit_MagnetizationScaling(T_vals,absM_av,L_vals);
    %% 3a Plot a
%     figure(i_model*10)
    figure
    hold on
    N_N = numel(sqrtN_vals);
    c_map = linspecer(N_N);
    max_chiabsm=zeros(size(sqrtN_vals));
    i_max_vec=zeros(size(sqrtN_vals));
    chiabsm=cell(size(sqrtN_vals));
    chierr=cell(size(sqrtN_vals));
    for i_N = 1 : N_N
%         chiabsm_vec=zeros(size(T_vals));
%         chiabsm_vec=N_vals(i_N)*(M_2_av(i_N,:) - absM_av(i_N,:).^2) ./T_vals;
        var_vec = M_2_av(i_N,:) - absM_av(i_N,:).^2;
        chiabsm_vec = N_vals(i_N) * var_vec ./ T_vals;
        chiabsm{i_N} = chiabsm_vec;

        varerr_vec=sqrt(1/runmax * (M_2_av(i_N,:) - (runmax - 3) / (runmax - 1) *var_vec.^2));
        % Apparently pg 483, Rao, 1973
        chierr{i_N} = N_vals(i_N) * varerr_vec ./ T_vals;

        [max_chiabsm(i_N),i_max_vec(i_N)] = max(T_vals.*chiabsm_vec);
        fitfunc=mag_fitfuncs{i_N};

%         hchiabsM_line(i_N) = line(T_vals,chiabsm_vec);
        hchiabsM_line(i_N) = line(T_vals,chiabsm_vec./sqrtN_vals(i_N).^(0));
        set(hchiabsM_line(i_N), ...
            'LineStyle', '-', 'LineWidth', 1.5, ...
            'HandleVisibility','off',...
            'Color', c_map(i_N,:))

        
%         hchiabsM_dots(i_N) = line(T_vals,chiabsm_vec);
        hchiabsM_dots(i_N) = line(T_vals,chiabsm_vec/sqrtN_vals(i_N)^(0));
        set(hchiabsM_dots(i_N), ...
            'LineStyle', 'none', ...
            'Marker', 'o', 'MarkerSize', 5, ...
            'LineWidth', .5, ...
            'HandleVisibility','off',...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:),...
            'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))

        plot(NaN,NaN, ...
            'LineStyle', '-', ...
            'Marker', 'o', 'MarkerSize', 5, ...
            'LineWidth', .5, ...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:),...
            'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
    
    end
    

%     for i_N = 1 : N_sqrtN
%         hMag_dots(i_N) = line(T_vals,absM_av(i_N,:));
%         set(hMag_dots(i_N), ...
%             'LineStyle', 'none', ...
%             'Marker', 'o', 'MarkerSize', 4, ...
%             'DisplayName', labels(i_N), ...
%             'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
%         
%         hcrossover_line(i_N) = line(T_star(i_N),crossover_M(i_N));
%         set(hcrossover_line(i_N), ...
%             'LineStyle', 'none', ...
%             'Marker', 's', 'MarkerSize', 10, ...
%             'Color', c_map(i_N,:), 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
%     end


     % 3a Legend, axes etc
%     if (i_model == 1)
%         hLegend = legend([hchiabsM_line(1), hchiabsM_line(2), hchiabsM_line(3), hchiabsM_line(4)], ...
%             '$N = 2^{10}$', '$N = 2^{12}$', '$N = 2^{14}$', '$N = 2^{16}$');
% %             'Location', 'SouthWest','interpreter','latex','FontSize', 20);
%     else
%         hLegend = legend([hchiabsM_line(1), hchiabsM_line(2), hchiabsM_line(3), hchiabsM_line(4), hchiabsM_line(5)], ...
%             '$N = (16)^{2}$', '$N = (32)^{2}$', '$N = (64)^{2}$', '$N = (128)^{2}$', '$N = (256)^{2}$');
%             'Location', 'SouthWest','interpreter','latex','FontSize', 20);
        
%     end
    hLegend = legend('Location', 'NorthWest','interpreter','latex');
    % hlegend(2).LineStyle = '-';
    % lineEntry = findobj(hLegend.EntryContainer, 'Object',hchiabsM_line(2));
    % entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    % entryMarker.LineWidth = .5;
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$\chi_{m}$','interpreter','latex');
%     hYLabel = ylabel('$\chi_{m}/\sqrt{N}$','interpreter','latex');
    h_axis = gca;
    xlim([0, T_max]);

    % 3a Font
    set_fonts_default


    % 3a Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5)
    ax_full = gca;
    
    % 3a Inset
    ax_inset = axes('Position',[.62 .62 .26 .26]);
    
%     max_chiabsm=cell2mat(max_chiabsm);
    chi_var_vec=zeros(size(sqrtN_vals));
    chi_mean_vec=zeros(size(sqrtN_vals));
    chi_err_vec=zeros(size(sqrtN_vals));
    for i_N = 1:numel(sqrtN_vals)
        i_T = find(max_chiabsm(i_N)/sqrtN_vals(i_N)^2 == (M_2_av(i_N,:) - absM_av(i_N,:).^2));
        curfile=sprintf('%s/%s/%s/%s',collectfile_dir,sqrtN_dirs(i_N),T_dirs{i_T},collectfilename);
        curdata=load(curfile,'M_2_av_collect','absM_av_collect');
        M_2_av_collect=curdata.('M_2_av_collect');
        absM_av_collect=curdata.('absM_av_collect');
        chi_mean_vec(i_N)=mean(M_2_av_collect - absM_av_collect.^2);
        chi_var_vec(i_N)=sqrt(var(M_2_av_collect - absM_av_collect.^2));
        chi_err_vec(i_N) = chierr{i_N}(i_max_vec(i_N)) * T_vals(i_max_vec(i_N)) / 2;
    end
%     hchimax_line = line(N_vals,max_chiabsm./N_vals);
%     chi_err_vec=chi_var_vec/sqrt(runmax);
    hchimax_line = errorbar(N_vals,max_chiabsm./N_vals,chi_err_vec./N_vals);
    f_POW = fittype('a*(x)^(2-eta)');
    init_exponent = .25;
    init_factor = 1;
    c0=[init_factor,init_exponent];
    Lower_lim=[0, 0];
    Upper_lim=[Inf, Inf];
    
%     c = fit(sqrtN_vals(1:end-1)', max_chiabsm(1:end-1)', f_POW,'StartPoint',c0,...
%         'Lower',Lower_lim, 'Upper', Upper_lim);
    c = fit(sqrtN_vals(1:end)', max_chiabsm(1:end)', f_POW,'StartPoint',c0,...
        'Lower',Lower_lim, 'Upper', Upper_lim,...
        'Algorithm','Trust-Region', 'Robust', 'LAR');

    f_POW_fixedeta = fittype('a*(x)^(2-.25)');
    init_factor = 1;
    c0=[init_factor];
    Lower_lim=[0];
    Upper_lim=[Inf];
    
    c_fixedeta = fit(sqrtN_vals(:), max_chiabsm(:), f_POW_fixedeta,'StartPoint',c0,...
        'Lower',Lower_lim, 'Upper', Upper_lim);
    NN=logspace(1,6);
    
    hchimax_fit = line(NN,c_fixedeta.a * sqrt(NN).^(- .25));
%     hchimax_fit = line(NN,c.a * sqrt(NN).^(- c.eta));
    
    c_map = linspecer(3);
    set(hchimax_line, ...
        'LineStyle', 'none', 'LineWidth', 1, ...
        'Marker', '^', 'MarkerSize', 5, ...
        'Color',c_map(1,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(1,:))
    
    set(ax_inset, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'XTick', [1e1 1e2 1e3 1e4], 'YTick', 0:1e-2:5e-2, ...
        'LineWidth', .5, 'Xscale', 'log','Yscale','lin')
    
    hXLabel = xlabel('$N$','interpreter','latex');
%     hYLabel = ylabel('$\chi_{\textrm{max}} / N$','interpreter','latex');
    hYLabel = ylabel('$\sigma_{\textrm{max}}^2$','interpreter','latex');
    
    set_fonts_default

    
    xlim([.5*N_vals(1),2*N_vals(end)])
    ylim([0 .05])
%     ylim([.006,.025])
    
     % 3a Text
%      annotation('textbox','String',sprintf('$\\eta = %.2f$',c.eta),...
%      annotation('textbox','String',sprintf('$\\eta = %.2f$',.25),...
     annotation('textbox','String',sprintf('$\\eta = %.2f$',.25),...
         'LineWidth',.5,...
         'Position',ax_inset.Position + [.07 .155 .0 .0],...
         'HorizontalAlignment','left','Vert','bottom',...
         'interpreter','latex','FitBoxToText','on',...
         'BackgroundColor','White',...
         'FontSize', fontsize_annotation, 'FontName', 'cmr12')

     set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
    
    if(saveswitch == 1)
        figname=sprintf('%s/%s_AbsSusceptibility_oldversion',basedir,curmodel);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    





















   
    %% 2 Plot sigma/M
    figure
    hold on
    N_N = numel(sqrtN_vals);
    c_map = linspecer(N_N);
    for i_N = 1 : N_N
%         chiabsm_vec=zeros(size(T_vals));
%         chiabsm_vec=(M_2_av(i_N,:) - absM_av(i_N,:).^2);
%         chiabsm{i_N} = chiabsm_vec;
%         max_chiabsm(i_N) = max(chiabsm_vec);
%         fitfunc=mag_fitfuncs{i_N};

        hchiabsM_line(i_N) = line(T_vals,sqrt((M_2_av(i_N,:)./absM_av(i_N,:).^2 - 1)));
        set(hchiabsM_line(i_N), ...
            'LineStyle', '-', 'LineWidth', 2.5, ...
            'HandleVisibility','off',...
            'Color', c_map(i_N,:))

        hchiabsM_dots(i_N) = line(T_vals,sqrt((M_2_av(i_N,:)./absM_av(i_N,:).^2 - 1)));
        set(hchiabsM_dots(i_N), ...
            'LineStyle', 'none', ...
            'Marker', 'o', 'MarkerSize', 6, ...
            'DisplayName', labels(i_N), ...
            'Color',c_map(i_N,:),...
            'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
        
    end
%     for i_N = 1 : N_sqrtN
%         hMag_dots(i_N) = line(T_vals,absM_av(i_N,:));
%         set(hMag_dots(i_N), ...
%             'LineStyle', 'none', ...
%             'Marker', 'o', 'MarkerSize', 4, ...
%             'DisplayName', labels(i_N), ...
%             'Color',c_map(i_N,:),'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
%         
%         hcrossover_line(i_N) = line(T_star(i_N),crossover_M(i_N));
%         set(hcrossover_line(i_N), ...
%             'LineStyle', 'none', ...
%             'Marker', 's', 'MarkerSize', 10, ...
%             'Color', c_map(i_N,:), 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', c_map(i_N,:))
%     end


     % 2 Legend, axes etc
    hlegend = legend('Location', 'NorthWest','interpreter','latex');
    hXLabel = xlabel('$T$','interpreter','latex');
    hYLabel = ylabel('$\sigma/\langle m\rangle$','interpreter','latex');
    xlim([0, T_max]);

    % 2 Font
    set_fonts_default



    % 2 Adjust axes properties
    set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'YTick', 0:2:100, 'YGrid', 'on', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
        'LineWidth', .5)
    ax_full = gca;

    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
 
    
    figname=sprintf('%s/%s_SigmaOverM',basedir,curmodel);
    fprintf('Creating figure %s\n',figname)
    if(saveswitch == 1)
       exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',dpi_png);
       exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end
    
    
    
    
    
        
    
% end