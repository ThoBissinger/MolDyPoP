function [i_T] = find_T_index(dataset_id,T)
%FIND_T_INDEX Finds the T index in the dataset
    if (strcmp(dataset_id,'LepriRuffo_mxy'))
        T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20];
    elseif (strcmp(dataset_id,'LepriRuffo_fmxy'))
        T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .185 .19 .195 .20 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
    elseif (strcmp(dataset_id,'LepriRuffo_xy'))
        % WRONG, BUT PROBABLY UNNECESSARY
        T_vals=[.10, .20, .30, .40, .50, .60, .70, .80, .90, 1.00];
    elseif (strcmp(dataset_id,'LepriRuffo_extended_mxy'))
        T_vals=[.07 .09 .11 .14 .15 .16 .165 .17 .175 .18 .185 .19 .20];
    elseif (strcmp(dataset_id,'eq_mxy'))
        % WRONG, BUT PROBABLY UNNECESSARY
        T_vals=[];
    elseif (strcmp(dataset_id,'eq_xy') || strcmp(dataset_id,'eq_xy_s'))
        % WRONG, BUT PROBABLY UNNECESSARY
        T_vals=[];
    elseif (strcmp(dataset_id,'dynamics_mxy') || strcmp(dataset_id,'dynamics_mxy_fullT') || strcmp(dataset_id,'dynamics_mxy_better_q'))
        T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .25 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
    elseif ( strcmp(dataset_id,'dynamics_mxy_LinearTime') || strcmp(dataset_id,'static_mxy') || strcmp(dataset_id,'dynamics_fmxy_LinearTime') || strcmp(dataset_id,'static_fmxy'))
        T_vals=[.01 .03 .05 .07 .09 .11 .13 .14 .15 .155 .16 .165 .17 .175 .18 .185 .19 .195 .20 .205 .21 .22 .23 .24 .27 .29 .31 .33 .35 .37 .40 .43 .46 .49 .52];
    elseif (strcmp(dataset_id,'dynamics_xy_LinearTime') || strcmp(dataset_id,'static_xy'))
        % WRONG, BUT PROBABLY UNNECESSARY
        T_vals=[];
    elseif (strcmp(dataset_id,'dynamics_xy_s_AdjustedTime') )
         T_vals = [0.85 0.89 0.91 0.93 0.95 1.00];
    end
    i_T = find(T_vals == T);
end

