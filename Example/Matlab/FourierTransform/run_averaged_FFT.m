function run_averaged_FFT(rootdir,matfilename,runmax,collectfilename, fitswitch)
%RUN_AVERAGED_FFT Calculates the FFT for each run and averages to give
%error estimates
    collectfilepath=sprintf('%s/%s.mat',rootdir,collectfilename);
    addpath('/home/user/thobi/Code/1911_XYModelBasedOnMHoefler/MatlabScripts/sim_analysis/FitFuncs')
    data_check = 0;
    if (isfile(collectfilepath))
        load(collectfilepath,'Sww_omega_0_collect');
        if ( exist('Sww_omega_0_collect','var') || fitswitch == 0)
            data_check = 1;
        else
            clear Sww_omega_0_collect;
        end
    end
    if (data_check == 1)
        load(collectfilepath,'qbin');
        load(collectfilepath,'averaging_times');
        load(collectfilepath,'om_vals');

        load(collectfilepath,'Sxx_collect');
        load(collectfilepath,'Syy_collect');
        load(collectfilepath,'Smperpmperp_collect');
        load(collectfilepath,'Smparmpar_collect');
        load(collectfilepath,'Stt_collect');
        load(collectfilepath,'Sww_collect');

        if (fitswitch == 1)
            load(collectfilepath,'Sxx_noiselevel_collect');
            load(collectfilepath,'Sxx_h_0_collect');
            load(collectfilepath,'Sxx_omega_0_collect');
            load(collectfilepath,'Sxx_gamma_collect');
            load(collectfilepath,'Sxx_fitob_collect');
            load(collectfilepath,'Syy_noiselevel_collect');
            load(collectfilepath,'Syy_h_0_collect');
            load(collectfilepath,'Syy_omega_0_collect');
            load(collectfilepath,'Syy_gamma_collect');
            load(collectfilepath,'Syy_fitob_collect');
            load(collectfilepath,'Smperpmperp_noiselevel_collect');
            load(collectfilepath,'Smperpmperp_h_0_collect');
            load(collectfilepath,'Smperpmperp_omega_0_collect');
            load(collectfilepath,'Smperpmperp_gamma_collect');
            load(collectfilepath,'Smperpmperp_fitob_collect');
            load(collectfilepath,'Smparmpar_noiselevel_collect');
            load(collectfilepath,'Smparmpar_h_0_collect');
            load(collectfilepath,'Smparmpar_omega_0_collect');
            load(collectfilepath,'Smparmpar_gamma_collect');
            load(collectfilepath,'Smparmpar_fitob_collect');
            load(collectfilepath,'Stt_noiselevel_collect');
            load(collectfilepath,'Stt_h_0_collect');
            load(collectfilepath,'Stt_omega_0_collect');
            load(collectfilepath,'Stt_gamma_collect');
            load(collectfilepath,'Stt_fitob_collect');
            load(collectfilepath,'Sww_noiselevel_collect');
            load(collectfilepath,'Sww_h_0_collect');
            load(collectfilepath,'Sww_omega_0_collect');
            load(collectfilepath,'Sww_gamma_collect');
            load(collectfilepath,'Sww_fitob_collect');
        end


        old_runmax = sum(find(Sxx_collect(:,1)));
        old_size = numel(Sxx_collect(:,1));

        Sxx_collect = [Sxx_collect,zeros(runmax-old_size,numel(qbin)*numel(averaging_times))];
        Syy_collect = [Syy_collect,zeros(runmax-old_size,numel(qbin)*numel(averaging_times))];
        Smperpmperp_collect = [Smperpmperp_collect,zeros(runmax-old_size,numel(qbin)*numel(averaging_times))];
        Smparmpar_collect = [Smparmpar_collect,zeros(runmax-old_size,numel(qbin)*numel(averaging_times))];
        Stt_collect = [Stt_collect,zeros(runmax-old_size,numel(qbin)*numel(averaging_times))];
        Sww_collect = [Sww_collect,zeros(runmax-old_size,numel(qbin)*numel(averaging_times))];

        
        
        if (fitswitch == 1)
            Sxx_noiselevel_collect = [Sxx_noiselevel_collect,zeros(runmax-old_size,numel(qbin))];
            Sxx_h_0_collect = [Sxx_h_0_collect,zeros(runmax-old_size,numel(qbin))];
            Sxx_omega_0_collect = [Sxx_omega_0_collect,zeros(runmax-old_size,numel(qbin))];
            Sxx_gamma_collect = [Sxx_gamma_collect,zeros(runmax-old_size,numel(qbin))];
            Sxx_fitob_collect(old_size+1:runmax,1:numel(qbin)) = {};
            Syy_noiselevel_collect = [Syy_noiselevel_collect,zeros(runmax-old_size,numel(qbin))];
            Syy_h_0_collect = [Syy_h_0_collect,zeros(runmax-old_size,numel(qbin))];
            Syy_omega_0_collect = [Syy_omega_0_collect,zeros(runmax-old_size,numel(qbin))];
            Syy_gamma_collect = [Syy_gamma_collect,zeros(runmax-old_size,numel(qbin))];
            Syy_fitob_collect(old_size+1:runmax,1:numel(qbin)) = {};
            Smperpmperp_noiselevel_collect = [Smperpmperp_noiselevel_collect,zeros(runmax-old_size,numel(qbin))];
            Smperpmperp_h_0_collect = [Smperpmperp_h_0_collect,zeros(runmax-old_size,numel(qbin))];
            Smperpmperp_omega_0_collect = [Smperpmperp_omega_0_collect,zeros(runmax-old_size,numel(qbin))];
            Smperpmperp_gamma_collect = [Smperpmperp_gamma_collect,zeros(runmax-old_size,numel(qbin))];
            Smperpmperp_fitob_collect(old_size+1:runmax,1:numel(qbin)) = {};
            Smparmpar_noiselevel_collect = [Smparmpar_noiselevel_collect,zeros(runmax-old_size,numel(qbin))];
            Smparmpar_h_0_collect = [Smparmpar_h_0_collect,zeros(runmax-old_size,numel(qbin))];
            Smparmpar_omega_0_collect = [Smparmpar_omega_0_collect,zeros(runmax-old_size,numel(qbin))];
            Smparmpar_gamma_collect = [Smparmpar_gamma_collect,zeros(runmax-old_size,numel(qbin))];
            Smparmpar_fitob_collect(old_size+1:runmax,1:numel(qbin)) = {};
            Stt_noiselevel_collect = [Stt_noiselevel_collect,zeros(runmax-old_size,numel(qbin))];
            Stt_h_0_collect = [Stt_h_0_collect,zeros(runmax-old_size,numel(qbin))];
            Stt_omega_0_collect = [Stt_omega_0_collect,zeros(runmax-old_size,numel(qbin))];
            Stt_gamma_collect = [Stt_gamma_collect,zeros(runmax-old_size,numel(qbin))];
            Stt_fitob_collect(old_size+1:runmax,1:numel(qbin)) = {};
            Sww_noiselevel_collect = [Sww_noiselevel_collect,zeros(runmax-old_size,numel(qbin))];
            Sww_h_0_collect = [Sww_h_0_collect,zeros(runmax-old_size,numel(qbin))];
            Sww_omega_0_collect = [Sww_omega_0_collect,zeros(runmax-old_size,numel(qbin))];
            Sww_gamma_collect = [Sww_gamma_collect,zeros(runmax-old_size,numel(qbin))];
            Sww_fitob_collect(old_size+1:runmax,1:numel(qbin)) = {};
        end
    else
        cur_matfile = sprintf('%s/run_1/output/%s.mat',rootdir,matfilename);
        load(cur_matfile,'qbin')
        load(cur_matfile,'averaging_times')

        Sxx_collect = zeros(runmax,numel(qbin)*numel(averaging_times));
        Syy_collect = zeros(runmax,numel(qbin)*numel(averaging_times));
        Smperpmperp_collect = zeros(runmax,numel(qbin)*numel(averaging_times));
        Smparmpar_collect = zeros(runmax,numel(qbin)*numel(averaging_times));
        Stt_collect = zeros(runmax,numel(qbin)*numel(averaging_times));
        Sww_collect = zeros(runmax,numel(qbin)*numel(averaging_times));

        gxx_start = zeros(size(qbin));
        gyy_start = zeros(size(qbin));
        gmperpmperp_start = zeros(size(qbin));
        gmparmpar_start = zeros(size(qbin));
        gtt_start = zeros(size(qbin));
        gww_start = zeros(size(qbin));

        if (fitswitch == 1)
            Sxx_noiselevel_collect = zeros(runmax,numel(qbin));
            Sxx_h_0_collect = zeros(runmax,numel(qbin));
            Sxx_omega_0_collect = zeros(runmax,numel(qbin));
            Sxx_gamma_collect = zeros(runmax,numel(qbin));
            Sxx_fitob_collect = cell(runmax,numel(qbin));
            Syy_noiselevel_collect = zeros(runmax,numel(qbin));
            Syy_h_0_collect = zeros(runmax,numel(qbin));
            Syy_omega_0_collect = zeros(runmax,numel(qbin));
            Syy_gamma_collect = zeros(runmax,numel(qbin));
            Syy_fitob_collect = cell(runmax,numel(qbin));
            Smperpmperp_noiselevel_collect = zeros(runmax,numel(qbin));
            Smperpmperp_h_0_collect = zeros(runmax,numel(qbin));
            Smperpmperp_omega_0_collect = zeros(runmax,numel(qbin));
            Smperpmperp_gamma_collect = zeros(runmax,numel(qbin));
            Smperpmperp_fitob_collect = cell(runmax,numel(qbin));
            Smparmpar_noiselevel_collect = zeros(runmax,numel(qbin));
            Smparmpar_h_0_collect = zeros(runmax,numel(qbin));
            Smparmpar_omega_0_collect = zeros(runmax,numel(qbin));
            Smparmpar_gamma_collect = zeros(runmax,numel(qbin));
            Smparmpar_fitob_collect = cell(runmax,numel(qbin));
            Stt_noiselevel_collect = zeros(runmax,numel(qbin));
            Stt_h_0_collect = zeros(runmax,numel(qbin));
            Stt_omega_0_collect = zeros(runmax,numel(qbin));
            Stt_gamma_collect = zeros(runmax,numel(qbin));
            Stt_fitob_collect = cell(runmax,numel(qbin));
            Sww_noiselevel_collect = zeros(runmax,numel(qbin));
            Sww_h_0_collect = zeros(runmax,numel(qbin));
            Sww_omega_0_collect = zeros(runmax,numel(qbin));
            Sww_gamma_collect = zeros(runmax,numel(qbin));
            Sww_fitob_collect = cell(runmax,numel(qbin));
        end

        old_runmax = 0;
    
    end

    for i_run = old_runmax+1:runmax
        cur_matfile = sprintf('%s/run_%i/output/%s.mat',rootdir,i_run,matfilename);
        disp(cur_matfile)
        load(cur_matfile,'gxx');
        load(cur_matfile,'gyy');
        load(cur_matfile,'gmperpmperp');
        load(cur_matfile,'gmparmpar');
        load(cur_matfile,'gtt');
        load(cur_matfile,'gww');
        for i_q = 1:numel(qbin)
            q_indices = (i_q):length(qbin):numel(qbin)*numel(averaging_times);
            % xx 
            [ft_vals,om_vals]=FT_correlation(averaging_times,gxx(q_indices), 0);
            Sxx_collect(i_run,q_indices) = ft_vals;
            if (fitswitch == 1)
                c = fit_DampedOscillator_FourierSpace(om_vals,real(ft_vals));
                Sxx_noiselevel_collect(i_run,i_q) = c.noiselevel;
                Sxx_h_0_collect(i_run,i_q) = c.h_0;
                Sxx_omega_0_collect(i_run,i_q) = c.omega_0;
                Sxx_gamma_collect(i_run,i_q) = c.gamma;
                Sxx_fitob_collect{i_run,i_q} = c;
            end
            gxx_start(i_q) = gxx_start(i_q) + gxx(q_indices(1));

            % yy
            [ft_vals,om_vals]=FT_correlation(averaging_times,gyy(q_indices), 0);
            Syy_collect(i_run,q_indices) = ft_vals;
            if (fitswitch == 1)
                c = fit_DampedOscillator_FourierSpace(om_vals,real(ft_vals));
                Syy_noiselevel_collect(i_run,i_q) = c.noiselevel;
                Syy_h_0_collect(i_run,i_q) = c.h_0;
                Syy_omega_0_collect(i_run,i_q) = c.omega_0;
                Syy_gamma_collect(i_run,i_q) = c.gamma;
                Syy_fitob_collect{i_run,i_q} = c;
            end
            gyy_start(i_q) = gyy_start(i_q) + gyy(q_indices(1));

            % mperp
            [ft_vals,om_vals]=FT_correlation(averaging_times,gmperpmperp(q_indices), 0);
            Smperpmperp_collect(i_run,q_indices) = ft_vals;
            if (fitswitch == 1)
                c = fit_DampedOscillator_FourierSpace(om_vals,real(ft_vals));
                Smperpmperp_noiselevel_collect(i_run,i_q) = c.noiselevel;
                Smperpmperp_h_0_collect(i_run,i_q) = c.h_0;
                Smperpmperp_omega_0_collect(i_run,i_q) = c.omega_0;
                Smperpmperp_gamma_collect(i_run,i_q) = c.gamma;
                Smperpmperp_fitob_collect{i_run,i_q} = c;
            end
            gmperpmperp_start(i_q) = gmperpmperp_start(i_q) + gmperpmperp(q_indices(1));
            
            % mpar
            [ft_vals,om_vals]=FT_correlation(averaging_times,gmparmpar(q_indices), 0);
            Smparmpar_collect(i_run,q_indices) = ft_vals;
            if (fitswitch == 1)
                c = fit_DampedOscillator_FourierSpace(om_vals,real(ft_vals));
                Smparmpar_noiselevel_collect(i_run,i_q) = c.noiselevel;
                Smparmpar_h_0_collect(i_run,i_q) = c.h_0;
                Smparmpar_omega_0_collect(i_run,i_q) = c.omega_0;
                Smparmpar_gamma_collect(i_run,i_q) = c.gamma;
                Smparmpar_fitob_collect{i_run,i_q} = c;
            end
            gmparmpar_start(i_q) = gmparmpar_start(i_q) + gmparmpar(q_indices(1));

            % tt
            [ft_vals,om_vals]=FT_correlation(averaging_times,gtt(q_indices), 0);
            Stt_collect(i_run,q_indices) = ft_vals;
            if (fitswitch == 1)
                c = fit_DampedOscillator_FourierSpace(om_vals,real(ft_vals));
                Stt_noiselevel_collect(i_run,i_q) = c.noiselevel;
                Stt_h_0_collect(i_run,i_q) = c.h_0;
                Stt_omega_0_collect(i_run,i_q) = c.omega_0;
                Stt_gamma_collect(i_run,i_q) = c.gamma;
                Stt_fitob_collect{i_run,i_q} = c;
            end
            gtt_start(i_q) = gtt_start(i_q) + gtt(q_indices(1));

            % ww
            [ft_vals,om_vals]=FT_correlation(averaging_times,gww(q_indices), 0);
            Sww_collect(i_run,q_indices) = ft_vals;
            if (fitswitch == 1)
                c = fit_DampedOscillator_FourierSpace(om_vals,real(ft_vals));
                Sww_noiselevel_collect(i_run,i_q) = c.noiselevel;
                Sww_h_0_collect(i_run,i_q) = c.h_0;
                Sww_omega_0_collect(i_run,i_q) = c.omega_0;
                Sww_gamma_collect(i_run,i_q) = c.gamma;
                Sww_fitob_collect{i_run,i_q} = c;
            end
            gww_start(i_q) = gww_start(i_q) + gww(q_indices(1));

            % saving data to be able to restart
%             if (fitswitch == 1)
%                 save(collectfilepath,...
%                     'Sxx_collect', 'Syy_collect', 'Smperpmperp_collect', 'Smparmpar_collect', ...
%                     'Stt_collect', 'Sww_collect', ...
%                     'Sxx_noiselevel_collect', 'Sxx_omega_0_collect', 'Sxx_gamma_collect', 'Sxx_fitob_collect', ...
%                     'Syy_noiselevel_collect', 'Syy_omega_0_collect', 'Syy_gamma_collect', 'Syy_fitob_collect', ...
%                     'Smperpmperp_noiselevel_collect', 'Smperpmperp_omega_0_collect', 'Smperpmperp_gamma_collect', 'Smperpmperp_fitob_collect', ...
%                     'Smparmpar_noiselevel_collect', 'Smparmpar_omega_0_collect', 'Smparmpar_gamma_collect', 'Smparmpar_fitob_collect', ....
%                     'Stt_noiselevel_collect', 'Stt_omega_0_collect', 'Stt_gamma_collect', 'Stt_fitob_collect', ....
%                     'Sww_noiselevel_collect', 'Sww_omega_0_collect', 'Sww_gamma_collect', 'Sww_fitob_collect', ...
%                     'om_vals','qbin','runmax','averaging_times');
%             else
%                 save(collectfilepath,...
%                     'Sxx_collect', 'Syy_collect', 'Smperpmperp_collect', 'Smparmpar_collect', ...
%                     'Stt_collect', 'Sww_collect', ...
%                     'om_vals','qbin','runmax','averaging_times');
%             end
        end
    end
    Sxx = mean(Sxx_collect);
    Syy = mean(Syy_collect);
    Smperpmperp = mean(Smperpmperp_collect);
    Smparmpar = mean(Smparmpar_collect);
    Stt = mean(Stt_collect);
    Sww = mean(Sww_collect);

    gxx_start = abs(gxx_start) / runmax;
    gyy_start = abs(gyy_start) / runmax;
    gmperpmperp_start = abs(gmperpmperp_start) / runmax;
    gmparmpar_start = abs(gmparmpar_start) / runmax;
    gtt_start = abs(gtt_start) / runmax;
    gww_start = abs(gww_start) / runmax;


    Sxx_noiselevel = zeros(size(qbin));
    Sxx_h_0 = zeros(size(qbin));
    Sxx_omega_0 = zeros(size(qbin));
    Sxx_gamma = zeros(size(qbin));
    Syy_noiselevel = zeros(size(qbin));
    Syy_h_0 = zeros(size(qbin));
    Syy_omega_0 = zeros(size(qbin));
    Syy_gamma = zeros(size(qbin));
    Smperpmperp_noiselevel = zeros(size(qbin));
    Smperpmperp_h_0 = zeros(size(qbin));
    Smperpmperp_omega_0 = zeros(size(qbin));
    Smperpmperp_gamma = zeros(size(qbin));
    Smparmpar_noiselevel = zeros(size(qbin));
    Smparmpar_h_0 = zeros(size(qbin));
    Smparmpar_omega_0 = zeros(size(qbin));
    Smparmpar_gamma = zeros(size(qbin));
    Stt_noiselevel = zeros(size(qbin));
    Stt_h_0 = zeros(size(qbin));
    Stt_omega_0 = zeros(size(qbin));
    Stt_gamma = zeros(size(qbin));
    Sww_noiselevel = zeros(size(qbin));
    Sww_h_0 = zeros(size(qbin));
    Sww_omega_0 = zeros(size(qbin));
    Sww_gamma = zeros(size(qbin));
    Sxx_fitob = cell(size(qbin));
    Syy_fitob = cell(size(qbin));
    Smperpmperp_fitob = cell(size(qbin));
    Smparmpar_fitob = cell(size(qbin));
    Stt_fitob = cell(size(qbin));
    Sww_fitob = cell(size(qbin));

    for i_q = 1:numel(qbin)
        q_indices = (i_q):length(qbin):numel(qbin)*numel(averaging_times);
        Sxx_fitob{i_q} = fit_DampedOscillator_FourierSpace(om_vals, real(Sxx(q_indices)),gxx_start(i_q));
        Sxx_noiselevel(i_q) = Sxx_fitob{i_q}.noiselevel;
        Sxx_h_0(i_q) = Sxx_fitob{i_q}.h_0;
        Sxx_omega_0(i_q) = Sxx_fitob{i_q}.omega_0;
        Sxx_gamma(i_q) = Sxx_fitob{i_q}.gamma;
        Syy_fitob{i_q} = fit_DampedOscillator_FourierSpace(om_vals, real(Syy(q_indices)),gyy_start(i_q));
        Syy_noiselevel(i_q) = Syy_fitob{i_q}.noiselevel;
        Syy_h_0(i_q) = Syy_fitob{i_q}.h_0;
        Syy_omega_0(i_q) = Syy_fitob{i_q}.omega_0;
        Syy_gamma(i_q) = Syy_fitob{i_q}.gamma;
        Smperpmperp_fitob{i_q} = fit_DampedOscillator_FourierSpace(om_vals, real(Smperpmperp(q_indices)),gmperpmperp_start(i_q));
        Smperpmperp_noiselevel(i_q) = Smperpmperp_fitob{i_q}.noiselevel;
        Smperpmperp_h_0(i_q) = Smperpmperp_fitob{i_q}.h_0;
        Smperpmperp_omega_0(i_q) = Smperpmperp_fitob{i_q}.omega_0;
        Smperpmperp_gamma(i_q) = Smperpmperp_fitob{i_q}.gamma;
        Smparmpar_fitob{i_q} = fit_DampedOscillator_FourierSpace(om_vals, real(Smparmpar(q_indices)),gmparmpar_start(i_q));
        Smparmpar_noiselevel(i_q) = Smparmpar_fitob{i_q}.noiselevel;
        Smparmpar_h_0(i_q) = Smparmpar_fitob{i_q}.h_0;
        Smparmpar_omega_0(i_q) = Smparmpar_fitob{i_q}.omega_0;
        Smparmpar_gamma(i_q) = Smparmpar_fitob{i_q}.gamma;
        Stt_fitob{i_q} = fit_DampedOscillator_FourierSpace(om_vals, real(Stt(q_indices)),gtt_start(i_q));
        Stt_noiselevel(i_q) = Stt_fitob{i_q}.noiselevel;
        Stt_h_0(i_q) = Stt_fitob{i_q}.h_0;
        Stt_omega_0(i_q) = Stt_fitob{i_q}.omega_0;
        Stt_gamma(i_q) = Stt_fitob{i_q}.gamma;
        Sww_fitob{i_q} = fit_DampedOscillator_FourierSpace(om_vals, real(Sww(q_indices)),gww_start(i_q));
        Sww_noiselevel(i_q) = Sww_fitob{i_q}.noiselevel;
        Sww_h_0(i_q) = Sww_fitob{i_q}.h_0;
        Sww_omega_0(i_q) = Sww_fitob{i_q}.omega_0;
        Sww_gamma(i_q) = Sww_fitob{i_q}.gamma;
    end
    if (fitswitch == 1)

        Sxx_noiselevel_mean = mean(Sxx_noiselevel_collect);
        Sxx_h_0_mean = mean(Sxx_h_0_collect);
        Sxx_omega_0_mean = mean(Sxx_omega_0_collect);
        Sxx_gamma_mean = mean(Sxx_gamma_collect);
        Syy_noiselevel_mean = mean(Syy_noiselevel_collect);
        Syy_h_0_mean = mean(Syy_h_0_collect);
        Syy_omega_0_mean = mean(Syy_omega_0_collect);
        Syy_gamma_mean = mean(Syy_gamma_collect);
        Smperpmperp_noiselevel_mean = mean(Smperpmperp_noiselevel_collect);
        Smperpmperp_h_0_mean = mean(Smperpmperp_h_0_collect);
        Smperpmperp_omega_0_mean = mean(Smperpmperp_omega_0_collect);
        Smperpmperp_gamma_mean = mean(Smperpmperp_gamma_collect);
        Smparmpar_noiselevel_mean = mean(Smparmpar_noiselevel_collect);
        Smparmpar_h_0_mean = mean(Smparmpar_h_0_collect);
        Smparmpar_omega_0_mean = mean(Smparmpar_omega_0_collect);
        Smparmpar_gamma_mean = mean(Smparmpar_gamma_collect);
        Stt_noiselevel_mean = mean(Stt_noiselevel_collect);
        Stt_h_0_mean = mean(Stt_h_0_collect);
        Stt_omega_0_mean = mean(Stt_omega_0_collect);
        Stt_gamma_mean = mean(Stt_gamma_collect);
        Sww_noiselevel_mean = mean(Sww_noiselevel_collect);
        Sww_h_0_mean = mean(Sww_h_0_collect);
        Sww_omega_0_mean = mean(Sww_omega_0_collect);
        Sww_gamma_mean = mean(Sww_gamma_collect);
    end


    Sxx_std = std(Sxx_collect);
    Syy_std = std(Syy_collect);
    Smperpmperp_std = std(Smperpmperp_collect);
    Smparmpar_std = std(Smparmpar_collect);
    Stt_std = std(Stt_collect);
    Sww_std = std(Sww_collect);

    if (fitswitch == 1)
        Sxx_noiselevel_std = std(Sxx_noiselevel_collect);
        Sxx_h_0_std = std(Sxx_h_0_collect);
        Sxx_omega_0_std = std(Sxx_omega_0_collect);
        Sxx_gamma_std = std(Sxx_gamma_collect);
        Syy_noiselevel_std = std(Syy_noiselevel_collect);
        Syy_h_0_std = std(Syy_h_0_collect);
        Syy_omega_0_std = std(Syy_omega_0_collect);
        Syy_gamma_std = std(Syy_gamma_collect);
        Smperpmperp_noiselevel_std = std(Smperpmperp_noiselevel_collect);
        Smperpmperp_h_0_std = std(Smperpmperp_h_0_collect);
        Smperpmperp_omega_0_std = std(Smperpmperp_omega_0_collect);
        Smperpmperp_gamma_std = std(Smperpmperp_gamma_collect);
        Smparmpar_noiselevel_std = std(Smparmpar_noiselevel_collect);
        Smparmpar_h_0_std = std(Smparmpar_h_0_collect);
        Smparmpar_omega_0_std = std(Smparmpar_omega_0_collect);
        Smparmpar_gamma_std = std(Smparmpar_gamma_collect);
        Stt_noiselevel_std = std(Stt_noiselevel_collect);
        Stt_h_0_std = std(Stt_h_0_collect);
        Stt_omega_0_std = std(Stt_omega_0_collect);
        Stt_gamma_std = std(Stt_gamma_collect);
        Sww_noiselevel_std = std(Sww_noiselevel_collect);
        Sww_h_0_std = std(Sww_h_0_collect);
        Sww_omega_0_std = std(Sww_omega_0_collect);
        Sww_gamma_std = std(Sww_gamma_collect);
    end


    disp(collectfilepath);
    
    if (fitswitch == 1)
        save(collectfilepath,'Sxx', 'Syy', 'Smperpmperp', 'Smparmpar', 'Stt', 'Sww', ...
            'Sxx_collect', 'Syy_collect', 'Smperpmperp_collect', 'Smparmpar_collect', ...
            'Stt_collect', 'Sww_collect', ...
            'Sxx_std', 'Syy_std', 'Smperpmperp_std', 'Smparmpar_std', 'Stt_std', 'Sww_std', ...
            'Sxx_noiselevel_collect', 'Sxx_h_0_collect', 'Sxx_omega_0_collect', 'Sxx_gamma_collect', 'Sxx_fitob_collect', ...
            'Syy_noiselevel_collect', 'Syy_h_0_collect', 'Syy_omega_0_collect', 'Syy_gamma_collect', 'Syy_fitob_collect', ...
            'Smperpmperp_noiselevel_collect', 'Smperpmperp_h_0_collect', 'Smperpmperp_omega_0_collect', 'Smperpmperp_gamma_collect', 'Smperpmperp_fitob_collect', ...
            'Smparmpar_noiselevel_collect', 'Smparmpar_h_0_collect', 'Smparmpar_omega_0_collect', 'Smparmpar_gamma_collect', 'Smparmpar_fitob_collect', ....
            'Stt_noiselevel_collect', 'Stt_h_0_collect', 'Stt_omega_0_collect', 'Stt_gamma_collect', 'Stt_fitob_collect', ....
            'Sww_noiselevel_collect', 'Sww_h_0_collect', 'Sww_omega_0_collect', 'Sww_gamma_collect', 'Sww_fitob_collect', ...
            'Sxx_noiselevel', 'Sxx_h_0', 'Sxx_omega_0', 'Sxx_gamma', 'Sxx_fitob', 'Syy_noiselevel', 'Syy_h_0', 'Syy_omega_0', 'Syy_gamma', 'Syy_fitob', ...
            'Smperpmperp_noiselevel', 'Smperpmperp_h_0', 'Smperpmperp_omega_0', 'Smperpmperp_gamma', 'Smperpmperp_fitob', ...
            'Smparmpar_noiselevel', 'Smparmpar_h_0', 'Smparmpar_omega_0', 'Smparmpar_gamma', 'Smparmpar_fitob', ...
            'Stt_noiselevel', 'Stt_h_0', 'Stt_omega_0', 'Stt_gamma', 'Stt_fitob', ...
            'Sww_noiselevel', 'Sww_h_0', 'Sww_omega_0', 'Sww_gamma', 'Sww_fitob', ...
            'Sxx_noiselevel_mean', 'Sxx_h_0_mean', 'Sxx_omega_0_mean', 'Sxx_gamma_mean', ...
            'Syy_noiselevel_mean', 'Syy_h_0_mean', 'Syy_omega_0_mean', 'Syy_gamma_mean', ...
            'Smperpmperp_noiselevel_mean', 'Smperpmperp_h_0_mean', 'Smperpmperp_omega_0_mean', 'Smperpmperp_gamma_mean', ...
            'Smparmpar_noiselevel_mean', 'Smparmpar_h_0_mean', 'Smparmpar_omega_0_mean', 'Smparmpar_gamma_mean', ...
            'Stt_noiselevel_mean', 'Stt_h_0_mean', 'Stt_omega_0_mean', 'Stt_gamma_mean', ...
            'Sww_noiselevel_mean', 'Sww_h_0_mean', 'Sww_omega_0_mean', 'Sww_gamma_mean', ...
            'Sxx_noiselevel_std', 'Sxx_h_0_std', 'Sxx_omega_0_std', 'Sxx_gamma_std', ...
            'Syy_noiselevel_std', 'Syy_h_0_std', 'Syy_omega_0_std', 'Syy_gamma_std', ....
            'Smperpmperp_noiselevel_std', 'Smperpmperp_h_0_std', 'Smperpmperp_omega_0_std', 'Smperpmperp_gamma_std', ...
            'Smparmpar_noiselevel_std', 'Smparmpar_h_0_std', 'Smparmpar_omega_0_std', 'Smparmpar_gamma_std', ...
            'Stt_noiselevel_std', 'Stt_h_0_std', 'Stt_omega_0_std', 'Stt_gamma_std', ...
            'Sww_noiselevel_std', 'Sww_h_0_std', 'Sww_omega_0_std', 'Sww_gamma_std', ...
            'om_vals','qbin','runmax','averaging_times');
    else
        save(collectfilepath,'Sxx', 'Syy', 'Smperpmperp', 'Smparmpar', 'Stt', 'Sww', ...
            'Sxx_collect', 'Syy_collect', 'Smperpmperp_collect', 'Smparmpar_collect', ...
            'Stt_collect', 'Sww_collect', ...
            'Sxx_std', 'Syy_std', 'Smperpmperp_std', 'Smparmpar_std', 'Stt_std', 'Sww_std', ...
            'Sxx_noiselevel', 'Sxx_h_0', 'Sxx_omega_0', 'Sxx_gamma', 'Sxx_fitob', 'Syy_noiselevel', 'Syy_h_0', 'Syy_omega_0', 'Syy_gamma', 'Syy_fitob', ...
            'Smperpmperp_noiselevel', 'Smperpmperp_h_0', 'Smperpmperp_omega_0', 'Smperpmperp_gamma', 'Smperpmperp_fitob', 'Smparmpar_noiselevel', 'Smparmpar_h_0', 'Smparmpar_omega_0', 'Smparmpar_gamma', 'Smparmpar_fitob', ...
            'Stt_noiselevel', 'Stt_h_0', 'Stt_omega_0', 'Stt_gamma', 'Stt_fitob', 'Sww_noiselevel', 'Sww_h_0', 'Sww_omega_0', 'Sww_gamma', 'Sww_fitob', ...
            'om_vals','qbin','runmax','averaging_times');
    end
end

