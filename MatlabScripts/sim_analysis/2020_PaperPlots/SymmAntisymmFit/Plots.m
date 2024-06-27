clear
cd /home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/2020_PaperPlots/SymmAntisymmFit
run ../initialization_script
saveswitch=1;
basedir=sprintf('%s/SymmAntisymmFit',fig_base);

sqrtN_vals = [16 32 64 128];
sampfilename="samp_Dynamics";
    
    
% model="mxy";
model="mxy";
if model=="mxy"
    curmodel="mxy";
    sqrtN_vals = [16 32 64 128 256];
    L_vals=[9.25,18.5,37,74,148];
    T_vals = [.03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25];
    T_dirs = {"T_.03" "T_.05" "T_.07" "T_.09" "T_.11" "T_.13" "T_.14" "T_.15" "T_.155" "T_.16" "T_.165" "T_.17" "T_.175" "T_.18" "T_.185" "T_.19" "T_.195" "T_.20" "T_.205" "T_.21" "T_.22" "T_.23" "T_.24" "T_.25"};
    dir_lintime="/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00";
    dir_reduced_dt="/data/scc/thobi/220201_ReducedSmapleStepDeltat/mxy_3.00";
    dir_longer_t="/data/scc/thobi/211201_LongerTime/mxy_3.00";
elseif model=="xy_s"
    model="xy_s";
    curmodel="xy";
    sqrtN_vals = [16 32 64 128];
    L_vals = sqrtN_vals;
    T_dirs = {"T_.10", "T_.20", "T_.30", "T_.40", "T_.50", "T_.60", "T_.70", "T_.80", "T_.85", "T_.91", "T_.95", "T_1.00"};
    T_vals = [.10, .20, .30, .40, .50, .60, .70, .80, .85, .91, .95, 1.00];
    dir_reduced_dt="/data/scc/thobi/220201_ReducedSmapleStepDeltat/xy_s";
end

%load(sprintf('%s/%s_dataset',basedir,model));

fitfunc_symm="exp(-gamma*x/2)*(cos(omega_1*x) + gamma/2/omega_1*sin(omega_1*x))";
fitfunc_asymm="exp(-gamma*x/2)*cos(omega_1*x)";



%% Fixed N and T, varying q
i_N=4;
T_select=[.11, .14, .17, .18, .20];
% T_select=[.85 .91 .95 1.00];
% T_select=[.50 .60 .70 .80 .85 .91 .95 1.00];
% T=.14;
for T=T_select
    i_T=find(T_vals==T);
    sqrtN=sqrtN_vals(i_N);
    T=T_vals(i_T);
    fprintf('T = %.3f ',T);
    f_asymm=@(t,om,ga) exp(-ga*t/2).*cos(om*t);
    f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
    curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
    if ~isfile(curfile)
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
        if ~isfile(curfile)
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
        end
    end
    load(curfile,"averaging_times","gmperpmperp","qbin");
    t=averaging_times;
    
    q_integers=find((qbin/qbin(1)>round(qbin/qbin(1))-1e-3) .* (qbin/qbin(1)<round(qbin/qbin(1))+1e-3));
    q_nonzero=find(gmperpmperp(1:numel(qbin)) ~= 0);
    q_integers=intersect(q_integers,q_nonzero);
    
    n_q = numel(qbin);
    n_q_int = numel(q_integers);
    ga_s_vals=zeros(1,n_q_int);
    om_s_vals=zeros(1,n_q_int);
    ga_a_vals=zeros(1,n_q_int);
    om_a_vals=zeros(1,n_q_int);
    
    q_vals=qbin(q_integers);
    for i_q_int = 1:numel(q_vals)
        i_q = q_integers(i_q_int);
        q = q_vals(i_q_int);
        fprintf('q = %.3f, ',q);
        cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
        fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
        fitob_asym=fit(t(:),cf(:),fitfunc_asymm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
        ga_s_vals(i_q_int)=fitob_sym.gamma;
        om_s_vals(i_q_int)=fitob_sym.omega_1;
        ga_a_vals(i_q_int)=fitob_asym.gamma;
        om_a_vals(i_q_int)=fitob_asym.omega_1;
    end
    fprintf('\n');
    fitfunc_pow="a*x^b";
    fitob_ga_s=fit(q_vals(:),ga_s_vals(:),fitfunc_pow, ...
        'StartPoint',[ga_s_vals(1)/q_vals(1)^2,2], ...
        'Lower',[0,0]);
    fitob_ga_a=fit(q_vals(:),ga_a_vals(:),fitfunc_pow, ...
        'StartPoint',[ga_a_vals(1)/q_vals(1)^2,2], ...
        'Lower',[0,0]);
    
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    
    plot(q_vals,ga_s_vals./q_vals.^2,'-r','DisplayName','t-symmetric fit'); hold on;
    plot(q_vals,ga_a_vals./q_vals.^2,'-b','DisplayName','t-asymmetric fit');
    hlegend=legend('Location','best','Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_labels);
    hxlabel=xlabel('$q$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
    hylabel=ylabel('$\gamma / q^2$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
    htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$',sqrtN,T),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    subtitle_txt=sprintf('$\\sigma_{\\gamma}^{\\textrm{sym}}=%.3g$, $\\sigma_{\\gamma}^{\\textrm{asym}}=%.3g$',fitob_ga_s.b,fitob_ga_a.b);
    hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_labels);
    
    
    figname=sprintf('%s/figs/%s_ga_sqrtN_%d_T_%.3f',basedir,curmodel,sqrtN,T);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
end

%% Compare fit
i_N=4;
fig_indices=[1:4];
if model == "mxy"
    T_select=[.11, .14, .17, .18, .20];
elseif model == "xy_s"
    T_select=[.10 .20 .30 .40];%.50 .60 .70 .80 .85 .91 .95 1.00];
end
% T_select=[.85 .91 .95 1.00];
% T=.14;
i_q_int=4;
for T=T_select
    
    i_T=find(T_vals==T);
    sqrtN=sqrtN_vals(i_N);
    T=T_vals(i_T);
    f_asymm=@(t,om,ga) exp(-ga*t/2).*cos(om*t);
    f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
%     q_vals=reshape(q_vals_collect(i_N,i_T,:),[1,numel(q_vals_collect(i_N,i_T,:))]);
    % om_symm_vals=reshape(omega_1_vals_symm(i_N,i_T,:),[1,numel(omega_1_vals_symm(i_N,i_T,:))]);
    % om_asymm_vals=reshape(omega_1_vals_asymm(i_N,i_T,:),[1,numel(omega_1_vals_asymm(i_N,i_T,:))]);
    % ga_symm_vals=reshape(gamma_vals_symm(i_N,i_T,:),[1,numel(gamma_vals_symm(i_N,i_T,:))]);
    % ga_asymm_vals=reshape(gamma_vals_asymm(i_N,i_T,:),[1,numel(gamma_vals_asymm(i_N,i_T,:))]);
    
    curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
%     curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
    if ~isfile(curfile)
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
        if ~isfile(curfile)
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
        end
    end
    load(curfile,"averaging_times","gmperpmperp","qbin");
    
    q_integers=find((qbin/qbin(1)>round(qbin/qbin(1))-1e-3) .* (qbin/qbin(1)<round(qbin/qbin(1))+1e-3));
    i_q = q_integers(i_q_int);
    q=qbin(i_q);
    
    t=averaging_times;
    cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
    fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
    fitob_asym=fit(t(:),cf(:),fitfunc_asymm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
    ga_s=fitob_sym.gamma;
    om_s=fitob_sym.omega_1;
    ga_a=fitob_asym.gamma;
    om_a=fitob_asym.omega_1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIGURE: simple fit vs sim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(1,fig_indices)
        figure
        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        plot(t,cf,'-r','DisplayName','Simulation Data'); hold on;
        plot(t,f_symm(t,om_s,ga_s),'-b','DisplayName','t-symmetric fit');
        plot(t,f_asymm(t,om_a,ga_a),'-g','DisplayName','t-asymmetric fit');
        xmax=min([max(t),20/ga_s,2*pi*4/(om_s)]);
        xlim([0 xmax]);
        hlegend=legend('Location','Northeast','Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        hxlabel=xlabel('$t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        hylabel=ylabel('$C_{m\perp}(q,t)$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$, $q=%.3f$',sqrtN,T,q),'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_titles);
        subtitle_txt={sprintf('symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_s,ga_s),...
            sprintf('a-symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_a,ga_a)};
        hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        
        
        figname=sprintf('%s/figs/%s_cf_sqrtN_%d_T_%.3f_q_%.3f',basedir,curmodel,sqrtN,T,q);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIGURE: Peak Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(2,fig_indices)
        figure
        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
        ga=fitob_sym.gamma;
        om=fitob_sym.omega_1;
        peaksep=pi/om;
        interval_edges=peaksep/2:peaksep:max(t);
        
        t_peak=zeros(1,numel(interval_edges)-1);
        y_peak=t_peak;
        for i_interval = 1:numel(interval_edges)-1
            tmin=interval_edges(i_interval);
            tmax=interval_edges(i_interval+1);
            index_interval = find((t<tmax).*(t>tmin));
            [maxval,maxarg] = max(abs(cf(index_interval)));
            t_peak(i_interval) = t(index_interval(maxarg));
            y_peak(i_interval) = maxval;
            plot(t(index_interval),cf(index_interval)); hold on;
        end
        plot(t_peak,[y_peak;-y_peak],'k*')
        plot(interval_edges,0*interval_edges,'k.')
    
    %     hlegend=legend('Location','Northeast','Interpreter','latex',...
    %         'FontName','cmr14','fontsize',fontsize_labels);
        hxlabel=xlabel('$t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        hylabel=ylabel('$C_{m\perp}(q,t)$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$, $q=%.3f$',sqrtN,T,q),'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_titles);
        subtitle_txt={sprintf('symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_s,ga_s),...
            sprintf('a-symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_a,ga_a)};
        hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        
        
        figname=sprintf('%s/figs/%s_cf_peaks_sqrtN_%d_T_%.3f_q_%.3f',basedir,curmodel,sqrtN,T,q);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIGURE 3: Peak height compare sim to fit 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(3,fig_indices)
        figure
        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
        ga=fitob_sym.gamma;
        om=fitob_sym.omega_1;
        peaksep=pi/om;
        interval_edges=peaksep/2:peaksep:max(t);
        cf_fit=f_symm(t,om,ga);
        t_peak=zeros(1,numel(interval_edges)-1);
        y_peak=t_peak;
        tfit_peak=t_peak;
        yfit_peak=t_peak;
        for i_interval = 1:numel(interval_edges)-1
            tmin=interval_edges(i_interval);
            tmax=interval_edges(i_interval+1);
            index_interval = find((t<tmax).*(t>tmin));
            [maxval,maxarg] = max(abs(cf(index_interval)));
            t_peak(i_interval) = t(index_interval(maxarg));
            y_peak(i_interval) = maxval;
%             [maxval,maxarg] = max(abs(cf_fit(index_interval)));
%             tfit_peak(i_interval) = t(index_interval(maxarg));
%             yfit_peak(i_interval) = maxval;
            
    %         plot(t(index_interval),cf(index_interval)); hold on;
        end
        semilogy(t,exp(-ga*t/2),'b--','DisplayName','fit peaks');
        hold on;
        semilogy(t_peak,y_peak,'rx--','DisplayName','simulation data peaks');
        
        
        
    
        hlegend=legend('Location','SouthWest','Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        hxlabel=xlabel('$t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        hylabel=ylabel('$C_{m\perp}(q,t)$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$, $q=%.3f$',sqrtN,T,q),'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_titles);
        subtitle_txt={sprintf('symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_s,ga_s),...
            sprintf('a-symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_a,ga_a)};
        hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        ylim([max(.1*min(y_peak),exp(-ga*t(end)/2)),1]);
        
        
        figname=sprintf('%s/figs/%s_cf_peakcompare_sqrtN_%d_T_%.3f_q_%.3f',basedir,curmodel,sqrtN,T,q);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


        xlim([0 min(max(t),25*peaksep)]);
        set(gca,'YScale','lin')
        legend('Location','Best')
        figname=sprintf('%s/figs/%s_cf_peakcompare_zoomed_sqrtN_%d_T_%.3f_q_%.3f',basedir,curmodel,sqrtN,T,q);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIGURE 4: gamma(t)/t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(4,fig_indices)
        figure
        set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        
        y_gafree=f_symm(t,om_s,0);
        ind_nonzero = find(abs(y_gafree) > .95); % I just need the peaks
        expga_vals=cf(ind_nonzero)./y_gafree(ind_nonzero);
        ga_vals_cur=-2*log(abs(y_peak))./t_peak;
        plot(t_peak,ga_vals_cur,'.b','DisplayName','time-dependent damping');
    
        hlegend=legend('Location','Best','Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        hxlabel=xlabel('$t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        hylabel=ylabel('$\gamma(t)/t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
        htitle=title(sprintf('$N=(%d)^2$, $T=%.3f$, $q=%.3f$',sqrtN,T,q),'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_titles);
        subtitle_txt={sprintf('symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_s,ga_s),...
            sprintf('a-symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om_a,ga_a)};
        hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
            'FontName','cmr14','fontsize',fontsize_labels);
        ylim([min(0,min(ga_vals_cur)) 1.2*max(ga_vals_cur)]);
        
        figname=sprintf('%s/figs/%s_ga_TD_sqrtN_%d_T_%.3f_q_%.3f',basedir,curmodel,sqrtN,T,q);
        fprintf('Creating figure %s\n',figname)
        exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
        exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    end

end








%% Size compare peak heights
if model == "mxy"
    T_select=[.11, .14, .17, .18, .20];
    sqrtN_select = [32 64 128];
elseif model == "xy_s"
    T_select=[.10 .20 .30 .40 ];%.50 .60 .70 .80 .85 .91 .95 1.00];
    sqrtN_select = [16 32 64 128];
end
% T_select=[.85 .91 .95 1.00];
% T=.14;
i_q_int=1;

for T=T_select
    
    i_T=find(T_vals==T);
    T=T_vals(i_T);
    f_asymm=@(t,om,ga) exp(-ga*t/2).*cos(om*t);
    f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
    q_cell=cell(1,numel(sqrtN_select));
    t_peak_cell=cell(1,numel(sqrtN_select));
    y_peak_cell=cell(1,numel(sqrtN_select));
    ga_vals=zeros(1,numel(sqrtN_select));
    om_vals=zeros(1,numel(sqrtN_select));

    for i_N = 1:numel(sqrtN_select)
        sqrtN=sqrtN_select(i_N);
        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN,T_dirs{i_T},sampfilename);
    %     curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
        if ~isfile(curfile)
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN,T_dirs{i_T},sampfilename);
            if ~isfile(curfile)
                curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN,T_dirs{i_T},sampfilename);
            end
        end
        load(curfile,"averaging_times","gmperpmperp","qbin");
    
        q_integers=find((qbin/qbin(1)>round(qbin/qbin(1))-1e-3) .* (qbin/qbin(1)<round(qbin/qbin(1))+1e-3));
        i_q = q_integers(i_q_int*2^(i_N-1));
        q=qbin(i_q);
        q_vals{i_N}=q;
        
        t=averaging_times;
        cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
        fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
        ga=fitob_sym.gamma;
        om=fitob_sym.omega_1;
        ga_vals(i_N)=ga;
        om_vals(i_N)=om;

        peaksep=pi/om;
        interval_edges=peaksep/2:peaksep:max(t);
        t_peak=zeros(1,numel(interval_edges)-1);
        y_peak=t_peak;
        tfit_peak=t_peak;
        yfit_peak=t_peak;
        for i_interval = 1:numel(interval_edges)-1
            tmin=interval_edges(i_interval);
            tmax=interval_edges(i_interval+1);
            index_interval = find((t<tmax).*(t>tmin));
            [maxval,maxarg] = max(abs(cf(index_interval)));
            t_peak(i_interval) = t(index_interval(maxarg));
            y_peak(i_interval) = maxval;
        end
        t_peak_cell{i_N}=t_peak;
        y_peak_cell{i_N}=y_peak;
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIGURE: Peak height compare sim to fit for different system sizes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
        
    marksym=["s","o","x","*","d"];
    for i_N = 1:numel(sqrtN_select)
        sqrtN=sqrtN_select(i_N);
        
        ga=ga_vals(i_N);
        t_peak=t_peak_cell{i_N};
        y_peak=y_peak_cell{i_N};
        
        dispname=sprintf('$N=(%d)^2$',sqrtN);
        semilogy(t_peak,y_peak, ...
            'LineStyle','none',...
            'Marker',marksym(i_N),'DisplayName',dispname);
        hold on;
    end
    semilogy(t_peak,exp(-ga_vals(1)*t_peak/2),'k-','DisplayName','fit peaks');
    
        

    hlegend=legend('Location','SouthWest','Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_labels, ...
        'NumColumns',2);
    hxlabel=xlabel('$t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
    hylabel=ylabel('$C_{m\perp}(q,t)$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
    htitle=title(sprintf('$T=%.3f$, $q=%.3f$',T,q),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    subtitle_txt={sprintf('symm: $\\omega_1=%.2g$, $\\gamma=%.2g$',om,ga)};
    hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_labels);
    ylim([max(.1*min(y_peak),exp(-ga*t(end)/2)),1]);
    
    
    figname=sprintf('%s/figs/%s_peak_sizecompare_T_%.3f_q_%.3f',basedir,curmodel,T,q);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');


    xlim([0 min(max(t),25*peaksep)]);
    set(gca,'YScale','lin')
    legend('Location','Best')
    figname=sprintf('%s/figs/%s_peak_sizecompare_zoomed_T_%.3f_q_%.3f',basedir,curmodel,T,q);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIGURE: gamma value compare sim to fit for different system sizes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
    c_map_cur=lines(numel(sqrtN_select));
    marksym=["s","o","x","*","d"];
%     newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
%     colororder(newcolors)
    for i_N = 1:numel(sqrtN_select)
        sqrtN=sqrtN_select(i_N);
        
        ga=ga_vals(i_N);
        t_peak=t_peak_cell{i_N};
        y_peak=y_peak_cell{i_N};
        
        dispname=sprintf('$N=(%d)^2$',sqrtN);
        plot(t_peak,-2*log(y_peak)./t_peak, ...
            'LineStyle','-','LineWidth',.7,...
            'Marker',marksym(i_N),'DisplayName',dispname, ...
            'MarkerSize',5,...
            'Color',c_map_cur(i_N,:));
        hold on;
    end
%     newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
%     colororder(newcolors)
    yline(NaN,'DisplayName','fit result','Color','k');
    for i_N = 1:numel(sqrtN_select)
        yline(ga_vals(i_N),'HandleVisibility','off','Color',c_map_cur(i_N,:), ...
            'LineWidth',1.2);
    end
    
        

    hlegend=legend('Location','Best','Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_labels, ...
        'NumColumns',2);
    hxlabel=xlabel('$t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
    hylabel=ylabel('$\gamma(t)/t$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
    htitle=title(sprintf('$T=%.3f$, $q=%.3f$',T,q),'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_titles);
    subtitle_txt={sprintf('$N=(%d)^2$: $\\omega_1=%.2g$, $\\gamma=%.2g$',sqrtN_select(1),om_vals(1),ga_vals(1))};
    hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
        'FontName','cmr14','fontsize',fontsize_labels);
    
    
    figname=sprintf('%s/figs/%s_peak_ga_sizecompare_T_%.3f_q_%.3f',basedir,curmodel,T,q);
    fprintf('Creating figure %s\n',figname)
    exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
    exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');

end












%% Fixed N, sigma error bars
% for i_N=3:5
%     for fittype = ["lin","log"]
%         for qmode = ["qint","qfull"]
T_select=[.11, .14, .16, .17, .175, .18, .185, .20];
T_select=[.50 .60 .70 .80 .85 .91 .95 1.00];
            
for i_N=[3,4]
    for fittype = ["log","lin"]
        for qmode = ["qint"]
            n_q_int = numel(q_integers);
            fitob_ga_s_cell=cell(1,numel(T_select));
            fitob_ga_a_cell=cell(1,numel(T_select));
            sigma_s_vals=zeros(1,numel(T_select));
            sigma_a_vals=zeros(1,numel(T_select));
            sigma_err_s_vals=zeros(1,numel(T_select));
            sigma_err_a_vals=zeros(1,numel(T_select));
            sigma_s_cimin_vals=zeros(1,numel(T_select));
            sigma_a_cimin_vals=zeros(1,numel(T_select));
            sigma_s_cimax_vals=zeros(1,numel(T_select));
            sigma_a_cimax_vals=zeros(1,numel(T_select));
            
            for i=1:numel(T_select)
                T=T_select(i);
                i_T=find(T_vals==T);
                sqrtN=sqrtN_vals(i_N);
                T=T_vals(i_T);
                f_asymm=@(t,om,ga) exp(-ga*t/2).*cos(om*t);
                f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
                curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                q_off=1;
                if ~isfile(curfile)
                    curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                    q_off=1;
                    if ~isfile(curfile)
                        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                        q_off=0;
                    end
                end
                load(curfile,"averaging_times","gmperpmperp","qbin");
                t=averaging_times;
                
                if qmode == "qint"
                    q_integers=find((qbin(1:end-q_off)/qbin(1)>round(qbin(1:end-q_off)/qbin(1))-1e-3) ...
                        .* (qbin(1:end-q_off)/qbin(1)<round(qbin(1:end-q_off)/qbin(1))+1e-3));
                elseif qmode == "qfull"
                    q_integers=1:numel(qbin)-q_off;
                end
                q_integers=find(gmperpmperp(q_integers)~=0);
                
                n_q = numel(qbin);
                n_q_int = numel(q_integers);
                ga_s_vals=zeros(1,n_q_int);
                om_s_vals=zeros(1,n_q_int);
                ga_a_vals=zeros(1,n_q_int);
                om_a_vals=zeros(1,n_q_int);
                
                q_vals=qbin(q_integers);
                for i_q_int = 1:numel(q_vals)
                    i_q = q_integers(i_q_int);
                    cf=gmperpmperp(i_q:numel(qbin):end)/gmperpmperp(i_q);
                    fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
                    fitob_asym=fit(t(:),cf(:),fitfunc_asymm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
                    ga_s_vals(i_q_int)=fitob_sym.gamma;
                    om_s_vals(i_q_int)=fitob_sym.omega_1;
                    ga_a_vals(i_q_int)=fitob_asym.gamma;
                    om_a_vals(i_q_int)=fitob_asym.omega_1;
                end
                if fittype == "lin"
                    fitfunc_pow="a*x^b";
                    [fitob_ga_s,gof_s,out_s]=fit(q_vals(:),ga_s_vals(:),fitfunc_pow, ...
                        'StartPoint',[ga_s_vals(1)/q_vals(1)^2,2], ...
                        'Lower',[0,0]);
                    [fitob_ga_a,gof_a,out_a]=fit(q_vals(:),ga_a_vals(:),fitfunc_pow, ...
                        'StartPoint',[ga_a_vals(1)/q_vals(1)^2,2], ...
                        'Lower',[0,0]);
                elseif fittype == "log"
                    fitfunc_pow="b*log(x)+a";
                    [fitob_ga_s,gof_s,out_s]=fit(q_vals(:),log(ga_s_vals(:)),fitfunc_pow, ...
                        'StartPoint',[log(ga_s_vals(1)/q_vals(1)^2),2], ...
                        'Lower',[-Inf,0]);
                    [fitob_ga_a,gof_a,out_a]=fit(q_vals(:),log(ga_a_vals(:)),fitfunc_pow, ...
                        'StartPoint',[log(ga_a_vals(1)/q_vals(1)^2),2], ...
                        'Lower',[-Inf,0]);
                end
                fitob_ga_s_cell{i}=fitob_ga_s;
                fitob_ga_a_cell{i}=fitob_ga_a;
                sigma_s_vals(i)=fitob_ga_s.b;
                sigma_a_vals(i)=fitob_ga_a.b;
                
                ci=confint(fitob_ga_s);
                sigma_s_cimin_vals(i)=ci(1,2);
                sigma_s_cimax_vals(i)=ci(2,2);
                sigma_err_s_vals(i)=(ci(2,2)-ci(1,2))/3.92;
                ci=confint(fitob_ga_a);
                sigma_err_a_vals(i)=(ci(2,2)-ci(1,2))/3.92;
                sigma_a_cimin_vals(i)=ci(1,2);
                sigma_a_cimax_vals(i)=ci(2,2);
            
                
            end
            
            figure
            set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
            set(gca,'FontSize',fontsize_axis);
            errorbar(T_select,sigma_s_vals,sigma_err_s_vals,'-r','DisplayName','t-symmetric fit'); hold on;
            errorbar(T_select,sigma_a_vals,sigma_err_a_vals,'-b','DisplayName','t-asymmetric fit');
            plot(T_select,sigma_s_cimin_vals,'vr','HandleVisibility','off');%,'DisplayName','t-symmetric fit, ci min'); 
            plot(T_select,sigma_s_cimax_vals,'^r','HandleVisibility','off');%,'DisplayName','t-symmetric fit, ci max'); 
            hlegend=legend('Location','best','Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_labels);
            hxlabel=xlabel('$T$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
            hylabel=ylabel('$\sigma_{\gamma}$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
            htitle=title(sprintf('$N=(%d)^2$',sqrtN),'Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_titles);
            
            subtitle_txt=sprintf('Fittype "%s", q selection "%s"',fittype,qmode);
            hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_labels);
            
            
            figname=sprintf('%s/figs/%s_ga_err_%sfit_%s_sqrtN_%d',basedir,curmodel,fittype,qmode,sqrtN);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
end



%% Fixed N, sigma error bars, only one fit
if model == "mxy"
    T_select=[.11, .14, .16, .17, .175, .18, .185, .20];
elseif model == "xy_s"
    T_select=[.50 .60 .70 .80 .85 .91 .95 1.00];
end

% for i_N=3:5
%     for fittype = ["lin","log"]
%         for qmode = ["qint","qfull"]
for i_N=[3:4]
    for fittype = ["lin","log"]
        for qmode = ["qint"]
            n_q_int = numel(q_integers);
            fitob_ga_s_cell=cell(1,numel(T_select));
            sigma_s_vals=zeros(1,numel(T_select));
            sigma_err_s_vals=zeros(1,numel(T_select));
            sigma_s_cimin_vals=zeros(1,numel(T_select));
            sigma_s_cimax_vals=zeros(1,numel(T_select));
            
            for i=1:numel(T_select)
                T=T_select(i);
                i_T=find(T_vals==T);
                sqrtN=sqrtN_vals(i_N);
                T=T_vals(i_T);
                fprintf('T = %.3f ',T);
    
                f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
                curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                q_off=1;
                if ~isfile(curfile)
                    curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                    q_off=1;
                    if ~isfile(curfile)
                        curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                        q_off=0;
                    end
                end
                load(curfile,"averaging_times","gmperpmperp","qbin");
                t=averaging_times;
                
                if qmode == "qint"
                    q_integers=find((qbin(1:end-q_off)/qbin(1)>round(qbin(1:end-q_off)/qbin(1))-1e-3) ...
                        .* (qbin(1:end-q_off)/qbin(1)<round(qbin(1:end-q_off)/qbin(1))+1e-3));
                elseif qmode == "qfull"
                    q_integers=1:numel(qbin)-q_off;
                end
                q_integers=find(gmperpmperp(q_integers)~=0);
                
                n_q = numel(qbin);
                n_q_int = numel(q_integers);
                ga_s_vals=zeros(1,n_q_int);
                om_s_vals=zeros(1,n_q_int);
                
                q_vals=qbin(q_integers);
                for i_q_int = 1:numel(q_vals)
                    i_q = q_integers(i_q_int);
                    fprintf('q = %.3f, ',q_vals(i_q_int));
                    cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
                    fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
                    ga_s_vals(i_q_int)=fitob_sym.gamma;
                    om_s_vals(i_q_int)=fitob_sym.omega_1;
                end
                if fittype == "lin"
                    fitfunc_pow="a*x^b";
                    [fitob_ga_s,gof_s,out_s]=fit(q_vals(:),ga_s_vals(:),fitfunc_pow, ...
                        'StartPoint',[ga_s_vals(1)/q_vals(1)^2,2], ...
                        'Lower',[0,0]);
                elseif fittype == "log"
                    fitfunc_pow="b*log(x)+a";
                    [fitob_ga_s,gof_s,out_s]=fit(q_vals(:),log(ga_s_vals(:)),fitfunc_pow, ...
                        'StartPoint',[log(ga_s_vals(1)/q_vals(1)^2),2], ...
                        'Lower',[-Inf,0]);
                end
                fitob_ga_s_cell{i}=fitob_ga_s;
                sigma_s_vals(i)=fitob_ga_s.b;
                
                ci=confint(fitob_ga_s);
                sigma_s_cimin_vals(i)=ci(1,2);
                sigma_s_cimax_vals(i)=ci(2,2);
                sigma_err_s_vals(i)=(ci(2,2)-ci(1,2))/3.92;
            
                
            end
            
            figure
            set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
            set(gca,'FontSize',fontsize_axis);
            plot(T_select,sigma_s_cimin_vals,'xr','DisplayName','ci min'); hold on;
            plot(T_select,sigma_s_cimax_vals,'xg','DisplayName','ci max'); 
            errorbar(T_select,sigma_s_vals,sigma_err_s_vals,'-k','DisplayName','$\sigma_{\gamma}$ fit'); hold on;
            
            hlegend=legend('Location','best','Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_labels);
            hxlabel=xlabel('$T$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
            hylabel=ylabel('$\sigma_{\gamma}$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
            htitle=title(sprintf('$N=(%d)^2$',sqrtN),'Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_titles);
            
            subtitle_txt=sprintf('Fittype "%s", q selection "%s"',fittype,qmode);
            hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_labels);
            
            
            figname=sprintf('%s/figs/%s_ga_fullT_err_%sfit_%s_sqrtN_%d',basedir,curmodel,fittype,qmode,sqrtN);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
end




%% Fixed N, gamme-error included, sigma error bars, only one fit
T_select=[.11, .14, .16, .17, .175, .18, .185, .20];
% T_select=[.50 .60 .70 .80 .85 .91 .95 1.00];
% T_select=[.70];

n_fit=10;
fittypes = ["lin","log"];
n_fittypes = numel(fittypes);
% for i_N=3:5
%     for fittype = ["lin","log"]
%         for qmode = ["qint","qfull"]
for i_N=[5]
    for qmode = ["qint"]
%         n_q_int = numel(q_integers);
        fitob_sigma_cell=cell(n_fit,n_fittypes,numel(T_select));
        sigma_s_vals=zeros(n_fit,n_fittypes,numel(T_select));
        sigma_err_s_vals=zeros(n_fit,n_fittypes,numel(T_select));
        sigma_s_cimin_vals=zeros(n_fit,n_fittypes,numel(T_select));
        sigma_s_cimax_vals=zeros(n_fit,n_fittypes,numel(T_select));
        
        for i=1:numel(T_select)
            T=T_select(i);
            i_T=find(T_vals==T);
            sqrtN=sqrtN_vals(i_N);
            T=T_vals(i_T);
            fprintf('T = %.3f ',T);

            f_symm=@(t,om,ga) exp(-ga*t/2).*(cos(om*t)+ga/2/om*sin(om*t));
            curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_reduced_dt,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
            q_off=1;
            if ~isfile(curfile)
                curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_longer_t,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                q_off=1;
                if ~isfile(curfile)
                    curfile=sprintf('%s/sqrtN_%d/%s/%s.mat',dir_lintime,sqrtN_vals(i_N),T_dirs{i_T},sampfilename);
                    q_off=0;
                end
            end
            load(curfile,"averaging_times","gmperpmperp","qbin");
            t=averaging_times;
            
            if qmode == "qint"
                q_integers=find((qbin(1:end-q_off)/qbin(1)>round(qbin(1:end-q_off)/qbin(1))-1e-3) ...
                    .* (qbin(1:end-q_off)/qbin(1)<round(qbin(1:end-q_off)/qbin(1))+1e-3));
            elseif qmode == "qfull"
                q_integers=1:numel(qbin)-q_off;
            end
            q_integers=find(gmperpmperp(q_integers)~=0);
            
            n_q = numel(qbin);
            n_q_int = numel(q_integers);
            
            ga_s_vals=zeros(1,n_q_int);
            ga_s_ste_vals=zeros(1,n_q_int);
            ga_s_cimin_vals=zeros(1,n_q_int);
            ga_s_cimax_vals=zeros(1,n_q_int);
            
            om_s_vals=zeros(1,n_q_int);
            om_s_ste_vals=zeros(1,n_q_int);
            om_s_cimin_vals=zeros(1,n_q_int);
            om_s_cimax_vals=zeros(1,n_q_int);
            
            q_vals=qbin(q_integers);
            for i_q_int = 1:numel(q_vals)
                i_q = q_integers(i_q_int);
                fprintf('q = %.3f, ',q_vals(i_q_int));
                cf=real(gmperpmperp(i_q:numel(qbin):end))/real(gmperpmperp(i_q));
                fitob_sym=fit(t(:),cf(:),fitfunc_symm,'StartPoint',[1e-2, 1e-2], 'Lower', [0, 0]);
                ga_s_vals(i_q_int)=fitob_sym.gamma;
                om_s_vals(i_q_int)=fitob_sym.omega_1;
                ci = confint(fitob_sym);
                
                ga_s_cimin_vals(i_q_int)=ci(1,1);
                ga_s_cimax_vals(i_q_int)=ci(2,1);
                ga_s_ste_vals(i_q_int)=(ci(2,1)-ci(1,1))/3.92;

                om_s_cimin_vals(i_q_int)=ci(1,2);
                om_s_cimax_vals(i_q_int)=ci(2,2);
                om_s_ste_vals(i_q_int)=(ci(2,2)-ci(1,2))/3.92;
                
            end

            fprintf('\n');
            i_select = find(ga_s_vals > 1e-9);
            q_select = q_vals(i_select);
            ga_select = ga_s_vals(i_select);
            ga_ste_select = ga_s_ste_vals(i_select);
            
            for i_fittype = 1:numel(fittypes)
                fittype = fittypes(i_fittype);
                for i_fit = 1:n_fit
                    
                    %                     ste_cur = ga_s_ste_vals;
                    ga_cur = ga_select + normrnd(0,ga_ste_select);    
                    if fittype == "lin"
                        fitfunc_pow="a*x^b";
                        [fitob,gof_s,out_s]=fit(q_select(:),ga_cur(:),fitfunc_pow, ...
                            'StartPoint',[ga_select(1)/q_select(1)^2,2], ...
                            'Lower',[0,0]);
                    elseif fittype == "log"
                        fitfunc_pow="b*log(x)+a";
                        [fitob,gof_s,out_s]=fit(q_select(:),log(ga_cur(:)),fitfunc_pow, ...
                            'StartPoint',[log(ga_select(1)/q_select(1)^2),2], ...
                            'Lower',[-Inf,0]);
                    end
                    fitob_sigma_cell{i_fit,i_fittype,i}=fitob;
                    sigma_s_vals(i_fit,i_fittype,i)=fitob.b;
                    
                    ci=confint(fitob);
                    sigma_s_cimin_vals(i_fit,i_fittype,i)=ci(1,2);
                    sigma_s_cimax_vals(i_fit,i_fittype,i)=ci(2,2);
                    sigma_err_s_vals(i_fit,i_fittype,i)=(ci(2,2)-ci(1,2))/3.92;
                end
            
                
            end
            
%             figure
%             set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
%             set(gca,'FontSize',fontsize_axis);
%             plot(T_select,mean(sigma_s_cimin_vals),'xr','DisplayName','ci min'); hold on;
%             plot(T_select,mean(sigma_s_cimax_vals),'xg','DisplayName','ci max'); 
%             errorbar(T_select,mean(sigma_s_vals),mean(sigma_err_s_vals)+std(sigma_s_vals)/sqrt(n_fit),'-k','DisplayName','$\sigma_{\gamma}$ fit'); hold on;
%             
%             hlegend=legend('Location','best','Interpreter','latex',...
%                 'FontName','cmr14','fontsize',fontsize_labels);
%             hxlabel=xlabel('$T$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
%             hylabel=ylabel('$\sigma_{\gamma}$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
%             htitle=title(sprintf('$N=(%d)^2$',sqrtN),'Interpreter','latex',...
%                 'FontName','cmr14','fontsize',fontsize_titles);
%             
%             subtitle_txt=sprintf('Fittype "%s", q selection "%s"',fittype,qmode);
%             hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
%                 'FontName','cmr14','fontsize',fontsize_labels);
%             
%             
%             figname=sprintf('%s/figs/%s_sigma_fullT_gaerr_%sfit_%s_sqrtN_%d',basedir,curmodel,fittype,qmode,sqrtN);
%             fprintf('Creating figure %s\n',figname)
%             exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
%             exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
        
        for i_fittype = 1:numel(fittypes)
            fittype = fittypes(i_fittype);
            shapemat=[n_fit,numel(T_select)];
            figure
            set(gcf,'units','centimeters','OuterPosition',[0 0 columnwidth_cm columnwidth_cm]);
            set(gca,'FontSize',fontsize_axis);
            sig_vals = mean(reshape(sigma_s_cimax_vals(:,i_fittype,:),shapemat));
%             semilogy(T_select, ...
%                 mean(reshape(sigma_s_cimax_vals(:,i_fittype,:)-sigma_s_cimin_vals(:,i_fittype,:),shapemat))  ...
%                     ./ sig_vals, ...
%                 'xr','DisplayName','confidence interval');
            semilogy(T_select, ...
                (max(reshape(sigma_s_cimax_vals(:,i_fittype,:),shapemat)) ...
                    - min(reshape(sigma_s_cimin_vals(:,i_fittype,:),shapemat)))...
                    ./ sig_vals, ...
                'xr','DisplayName','confidence interval');
            hold on;
            semilogy(T_select, ...
                max(reshape(sigma_err_s_vals(:,i_fittype,:),shapemat)) ./ sig_vals, ...
                '*b','DisplayName','standard error');
            semilogy(T_select, ...
                std(reshape(sigma_s_vals(:,i_fittype,:),shapemat))/sqrt(n_fit) ./ sig_vals, ...
                'og','DisplayName','$\gamma$ uncertainty');
            
            hlegend=legend('Location','southwest','Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_labels);
            hxlabel=xlabel('$T$','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
            hylabel=ylabel('relative $\sigma$ error','Interpreter','latex','FontName','cmr14','fontsize',fontsize_labels);
            htitle=title(sprintf('$N=(%d)^2$',sqrtN),'Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_titles);

            ylim_vals=get(gca,'ylim');
            ylim([ylim_vals(1)/1000,ylim_vals(2)]);
            set(gca,'YTick',logspace(-6,1,8));

            subtitle_txt=sprintf('Fittype "%s", q selection "%s"',fittype,qmode);
            hsubtitle=subtitle(subtitle_txt,'Interpreter','latex',...
                'FontName','cmr14','fontsize',fontsize_labels);
            
            
            figname=sprintf('%s/figs/%s_sigma_fullT_errcompare_%sfit_%s_sqrtN_%d',basedir,curmodel,fittype,qmode,sqrtN);
            fprintf('Creating figure %s\n',figname)
            exportgraphics(gcf,sprintf('%s.png',figname),'Resolution',300);
            exportgraphics(gcf,sprintf('%s.pdf',figname),'ContentType','vector');
        end
    end
end