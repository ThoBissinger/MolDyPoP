/*! \file parameters.cpp
 *
 *  \brief cpp-File to class declaration of parameters. Implements the routines defined
 *  in parameters.h. For details, check there.
 *
 *  \date Created: 2019-04-12
 *  \date Last Updated: 2023-08-02
 *
 */
#include <vector>
#include "parameters.h"

//um xy wg zu bekommen

#include "computations.h"//<-- wichtig





//void parameters::read_from_file(std::ifstream& infile, std::ofstream& stdoutfile, int& errorcode){
int parameters::read_from_file(std::ifstream& infile){
//	set_to_default();
	std::string line;
	while(getline(infile,line)) {
		if (line.rfind("system =", 0) == 0 ){
			system_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("kT =", 0) == 0 ){
			kT_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("randomseed =", 0) == 0 ){
			randomseed_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("sqrtN =", 0) == 0 ){
			sqrtN_ = std::stod(line.erase(0,line.find("=")+2));
			N_ =  sqrtN_ * sqrtN_;
		} else if (line.rfind("N =", 0) == 0 ){
			N_ = std::stod(line.erase(0,line.find("=")+2));
			sqrtN_ = (int)sqrt(N_);
		} else if (line.rfind("L =", 0) == 0 ){
			L_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("mode =", 0) == 0 ){
			mode_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("job_id =", 0) == 0 ){
			if ( (line == "job_id = " ) | (line == "job_id = ''" ) | (line == "job_id = \"\"" ) ){
				job_id_ = "";
			} else {
				job_id_ = line.erase(0,line.find("=")+2);
			}

		} else if (line.rfind("J =", 0) == 0 ){
			J_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("U =", 0) == 0 ){
			U_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("cutoff =", 0) == 0 ){
			cutoff_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("I =", 0) == 0 ){
			I_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("m =", 0) == 0 ){
			m_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("Tmax =", 0) == 0 ){
			Tmax_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("dt =", 0) == 0 ){
			dt_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("vm_v =", 0) == 0 ){
			vm_v_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("vm_eta =", 0) == 0 ){
			vm_eta_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("activity =", 0) == 0 ){
			activity_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("lattice_type =", 0) == 0 ){
			lattice_type_ = line.at(line.find("=")+2);

		} else if (line.rfind("init_mode =", 0) == 0 ){
			init_mode_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("init_file =", 0) == 0 ){
			init_file_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("init_kT =", 0) == 0 ){
			init_kT_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("init_random_displacement =", 0) == 0 ){
			init_random_displacement_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("init_random_angle =", 0) == 0 ){
			init_random_angle_ = std::stod(line.erase(0,line.find("=")+2));


		} else if (line.rfind("N_rbin =", 0) == 0 ){
			N_rbin_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("N_qbin =", 0) == 0 ){
			N_qbin_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("min_binwidth_r =", 0) == 0 ){
			min_binwidth_r_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("qbin_type =", 0) == 0 ){
			qbin_type_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("qmax =", 0) == 0 ){
			qmax_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("min_binwidth_q =", 0) == 0 ){
			min_binwidth_q_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("qsamps_per_bin =", 0) == 0 ){
			qsamps_per_bin_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("n_rsamps =", 0) == 0 ){
			n_rsamps_ = std::stod(line.erase(0,line.find("=")+2));

		} else if (line.rfind("sample_integrator_type =", 0) == 0 ){
			sample_integrator_type_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("ensemble =", 0) == 0 ){
			ensemble_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("nhnp_pi =", 0) == 0 ){
			nhnp_pi_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("nh_eta =", 0) == 0 ){
			nh_eta_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("np_s =", 0) == 0 ){
			np_s_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("nhnp_Q =", 0) == 0 ){
			nhnp_Q_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("nhnp_tau =", 0) == 0 ){
			nhnp_tau_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("mc_steplength_theta =", 0) == 0 ){
			mc_steplength_theta_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("brownian_kT_omega =", 0) == 0 ){
			brownian_kT_omega_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("brownian_kT_p =", 0) == 0 ){
			brownian_kT_p_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("brownian_timestep =", 0) == 0 ){
			brownian_timestep_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("gamma_ld_p =", 0) == 0 ){
			gamma_ld_p_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("gamma_ld_om =", 0) == 0 ){
			gamma_ld_om_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("mc_steplength_r =", 0) == 0 ){
			mc_steplength_r_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("sampling_time_sequence =", 0) == 0 ){
			sampling_time_sequence_ = line.erase(0,line.find("=")+2);

		} else if (line.rfind("eq_mode =", 0) == 0 ){
			eq_mode_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("eq_integrator_type =", 0) == 0 ){
			eq_integrator_type_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("eq_Tmax =", 0) == 0 ){
			eq_Tmax_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_breakcond =", 0) == 0 ){
			eq_breakcond_ = line.erase(0,line.find("=")+2);
		} else if (line.rfind("eq_agreement_threshold =", 0) == 0 ){
			eq_agreement_threshold_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_av_time =", 0) == 0 ){
			eq_av_time_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("tau_berendsen =", 0) == 0 ){
			tau_berendsen_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_anneal_rate =", 0) == 0 ){
			eq_anneal_rate_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_anneal_step =", 0) == 0 ){
			eq_anneal_step_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_Tprintstep =", 0) == 0 ){
			eq_Tprintstep_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_brownian_kT_p =", 0) == 0 ){
			eq_brownian_kT_p_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_brownian_timestep =", 0) == 0 ){
			eq_brownian_timestep_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_brownian_kT_omega =", 0) == 0 ){
			eq_brownian_kT_omega_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_sampswitch =", 0) == 0 ){
			eq_sampswitch_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_Nsamp =", 0) == 0 ){
			eq_Nsamp_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("eq_samp_time_sequence =", 0) == 0 ){
			eq_samp_time_sequence_ = line.erase(0,line.find("=")+2);

		} else if (line.rfind("samplestart =", 0) == 0 ){
			samplestart_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("samplestep =", 0) == 0 ){
			samplestep_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("av_time_spacing =", 0) == 0 ){
			av_time_spacing_ = std::stod(line.erase(0,line.find("=")+2));
        } else if (line.rfind("Nsamp =", 0) == 0 ){
                Nsamp_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("print_snapshots =", 0) == 0 ){
			print_snapshots_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("on_fly_sampling =", 0) == 0 ){
			on_fly_sampling_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("snap_overview_file =", 0) == 0 ){
			snap_overview_file_ = line.erase(0,line.find("=")+2);
//		} else if (line.rfind("printoption =", 0) == 0 ){
//			printoption_ = std::stod(line.erase(0,line.find("=")+2));
		} else if (line.rfind("output_folder =", 0) == 0 ){//added by David
			output_folder_ = line.erase(0,line.find("=")+2);
		}
	}
	if (system_ == "xy" || system_ == "fvm" ){
		L_ = sqrtN_;
	}
//	if ( ( job_id_ != "") && (job_id_.at(0) != '_')){
	if ( job_id_ != ""){
		if (job_id_.at(0) != '_'){
			job_id_ = "_" + job_id_;
		}
	}
	outfilename_ = output_folder_ + "/data" + job_id_ + ".out";
	std::ofstream stdoutfile(outfilename_, std::ios::out);
	int errorcode = correct_values(stdoutfile);
	initialize_bins();
//	if (printoption_ != 3) {
	stdoutfile << "==========================================================\n";
	stdoutfile << "==================  Simulation data ======================\n";
	stdoutfile << "==========================================================\n";
	stdoutfile.precision(7);
	stdoutfile << std::setw(8) << "system:                   "   		<< system_ 						<< std::endl;
	stdoutfile << std::setw(8) << "mode:                     "			<< mode_						<< std::endl;
	stdoutfile << std::setw(8) << "job_id:                   "   		<< job_id_ 						<< std::endl;
	stdoutfile << std::setw(8) << "random seed:              "			<< randomseed_					<< std::endl;

	stdoutfile << std::setw(8) << "N:                        "			<< N_  							<< std::endl;
	stdoutfile << std::setw(8) << "sqrtN:                    "			<< sqrtN_  						<< std::endl;
	stdoutfile << std::setw(8) << "L:                        "			<< L_							<< std::endl;
	stdoutfile << std::setw(8) << "kT:                       "			<< kT_							<< std::endl;
	stdoutfile << std::setw(8) << "m:                        "			<< m_							<< std::endl;
	stdoutfile << std::setw(8) << "I:                        "			<< I_							<< std::endl;
	stdoutfile << std::setw(8) << "J:                        "			<< J_							<< std::endl;
	stdoutfile << std::setw(8) << "U:                        "			<< U_							<< std::endl;
	stdoutfile << std::setw(8) << "cutoff:                   "			<< cutoff_						<< std::endl;
	stdoutfile << std::setw(8) << "m:                        "			<< m_							<< std::endl;
	stdoutfile << std::setw(8) << "dt:                       "			<< dt_							<< std::endl;
	stdoutfile << std::setw(8) << "Tmax:                     "			<< Tmax_						<< std::endl;

	stdoutfile << std::setw(8) << "init_mode:                "			<< init_mode_ 					<< std::endl;
	stdoutfile << std::setw(8) << "init_file:                "			<< init_file_ 					<< std::endl;
	stdoutfile << std::setw(8) << "init_kT:                  " 			<< init_kT_ 					<< std::endl;
	stdoutfile << std::setw(8) << "init_random_displacement: "			<< init_random_displacement_ 	<< std::endl;
	stdoutfile << std::setw(8) << "init_random_angle:        "			<< init_random_angle_ 			<< std::endl;


	stdoutfile << std::setw(8) << "eq_mode:                  "  		<< eq_mode_ 					<< std::endl;
	stdoutfile << std::setw(8) << "eq_Tmax:                  " 			<< eq_Tmax_ 					<< std::endl;
	stdoutfile << std::setw(8) << "eq_breakcond:             "			<< eq_breakcond_				<< std::endl;
	stdoutfile << std::setw(8) << "eq_agreement_threshold:   "   		<< eq_agreement_threshold_ 		<< std::endl;
	stdoutfile << std::setw(8) << "eq_av_time:               "       	<< eq_av_time_ 					<< std::endl;
	stdoutfile << std::setw(8) << "tau_berendsen:            "    		<< tau_berendsen_ 				<< std::endl;
	stdoutfile << std::setw(8) << "eq_anneal_rate:           "   		<< eq_anneal_rate_		 		<< std::endl;
	stdoutfile << std::setw(8) << "eq_anneal_step:           "   		<< eq_anneal_step_ 				<< std::endl;
	stdoutfile << std::setw(8) << "eq_brownian_kT_omega:     "			<< eq_brownian_kT_omega_		<< std::endl;
	stdoutfile << std::setw(8) << "eq_brownian_kT_p:         "			<< eq_brownian_kT_p_			<< std::endl;
	stdoutfile << std::setw(8) << "eq_brownian_timestep:     "			<< eq_brownian_timestep_		<< std::endl;
	if (! std::isnan(eq_Tprintstep_)){
		stdoutfile << std::setw(8) << "eq_Tprintstep:            "			<< eq_Tprintstep_				<< std::endl;

	}
	stdoutfile << std::setw(8) << "eq_sampswitch:            "			<< eq_sampswitch_				<< std::endl;
	stdoutfile << std::setw(8) << "eq_Nsamp:                 "			<< eq_Nsamp_					<< std::endl;
	stdoutfile << std::setw(8) << "eq_samp_time_sequence:    "			<< eq_samp_time_sequence_		<< std::endl;

	stdoutfile << std::setw(8) << "N_rbin:                   "			<< N_rbin_						<< std::endl;
	stdoutfile << std::setw(8) << "N_qbin:                   "			<< N_qbin_						<< std::endl;
	stdoutfile << std::setw(8) << "min_binwidth_r:           "			<< min_binwidth_r_				<< std::endl;
	stdoutfile << std::setw(8) << "qbin_type:                "			<< qbin_type_					<< std::endl;
	stdoutfile << std::setw(8) << "qmax:                     "			<< qmax_						<< std::endl;
	stdoutfile << std::setw(8) << "min_binwidth_q:           "			<< min_binwidth_q_				<< std::endl;
	stdoutfile << std::setw(8) << "qsamps_per_bin:           "			<< qsamps_per_bin_				<< std::endl;
	stdoutfile << std::setw(8) << "n_rsamps:                 "			<< n_rsamps_ 					<< std::endl;
	stdoutfile << std::setw(8) << "lattice_type:             "			<< lattice_type_				<< std::endl;
	stdoutfile << std::setw(8) << "activity:                 "			<< activity_					<< std::endl;
	stdoutfile << std::setw(8) << "vm_v:                     "			<< vm_v_						<< std::endl;
	stdoutfile << std::setw(8) << "vm_eta:                   "			<< vm_eta_						<< std::endl;
	stdoutfile << std::setw(8) << "eq_integrator_type:       "			<< eq_integrator_type_			<< std::endl;
	stdoutfile << std::setw(8) << "sample_integrator_type:   "			<< sample_integrator_type_		<< std::endl;
	stdoutfile << std::setw(8) << "ensemble:                 "			<< ensemble_	  				<< std::endl;
	stdoutfile << std::setw(8) << "nhnp_pi:                  "    		<< nhnp_pi_ 					<< std::endl;
	stdoutfile << std::setw(8) << "nhnp_Q:                   "    		<< nhnp_Q_ 						<< std::endl;
	stdoutfile << std::setw(8) << "nhnp_tau:                 "    		<< nhnp_tau_ 					<< std::endl;
	stdoutfile << std::setw(8) << "nh_eta:                   "    		<< nh_eta_ 						<< std::endl;
	stdoutfile << std::setw(8) << "mc_steplength_theta:      "			<< mc_steplength_theta_			<< std::endl;
	stdoutfile << std::setw(8) << "brownian_kT_omega:        "			<< brownian_kT_omega_			<< std::endl;
	stdoutfile << std::setw(8) << "brownian_kT_p:            "			<< brownian_kT_p_				<< std::endl;
	stdoutfile << std::setw(8) << "brownian_timestep:        "			<< brownian_timestep_			<< std::endl;
	stdoutfile << std::setw(8) << "gamma_ld_p:               "			<< gamma_ld_p_					<< std::endl;
	stdoutfile << std::setw(8) << "gamma_ld_om:              "			<< gamma_ld_om_					<< std::endl;
	stdoutfile << std::setw(8) << "mc_steplength_r:          "			<< mc_steplength_r_				<< std::endl;


	stdoutfile << std::setw(8) << "sampling_time_sequence:   "   		<< sampling_time_sequence_		<< std::endl;
	stdoutfile << std::setw(8) << "samplestart:              "			<< samplestart_  				<< std::endl;
	stdoutfile << std::setw(8) << "samplestep:               "			<< samplestep_  				<< std::endl;
	stdoutfile << std::setw(8) << "av_time_spacing:          "			<< av_time_spacing_				<< std::endl;
	stdoutfile << std::setw(8) << "Nsamp:                    "			<< Nsamp_						<< std::endl;
	stdoutfile << std::setw(8) << "print_snapshots:          "			<< print_snapshots_  			<< std::endl;
	stdoutfile << std::setw(8) << "on_fly_sampling:          "			<< on_fly_sampling_  			<< std::endl;
	stdoutfile << std::setw(8) << "snap_overview_file:       "			<< snap_overview_file_			<< std::endl;

//		stdoutfile << std::setw(8) << "printoption:              "			<< printoption_  			<< std::endl;
	stdoutfile << "==========================================================\n\n";
//	}
	stdoutfile.close();

	return errorcode;
}


int parameters::correct_values(std::ofstream& stdoutfile){
	int errorcode = 0;
	// system_ input error check
	if ((system_ != "xy") && (system_ != "mxy") && (system_ != "fmxy") && (system_ != "vm") && (system_ != "fvm")){
		stdoutfile << "== ERROR" << std::endl;
		stdoutfile << "== System type " << system_ << " unknown." << std::endl;
		stdoutfile << " = Known types:" << std::endl;
		stdoutfile << " = 'xy'    : XY model" << std::endl;
		stdoutfile << " = 'mxy'   : mobile XY model (including Bore type models)" << std::endl;
		stdoutfile << " = 'fmxy'  : frozen mobile XY model (reads out initial positions and freezes them, incorporating mxy spin interactions)" << std::endl;
		stdoutfile << " = 'vm'    : Vicsek model" << std::endl;
		stdoutfile << " = 'fvm'   : Frozen Vicsek model" << std::endl;
		errorcode = 1;
	} else if ((system_ != "xy") && (system_ != "mxy") && (system_ != "fmxy") && (system_ != "vm") && (system_ != "fvm")) {
		stdoutfile << "== ERROR" << std::endl;
		stdoutfile << "== System type " << system_ << " not yet implemented." << std::endl;
		stdoutfile << " = Aborting calculation." << std::endl;
		errorcode = 1;
	}
	// Setting degrees of freedom (necessary for temperature calculations)
	if (system_ == "xy" || system_ == "fvm" || (system_ == "fmxy") ){
		dof_ = N_;
	} else if (system_ == "mxy" || system_ == "vm" ) {
		dof_ = 3 * N_;
	}
	// lattice_type_ input error check for lattice based systems
	if ( system_ == "xy" || system_ == "fvm") {
		if ( (lattice_type_ != 't' ) && (lattice_type_ != 's' )) {
			stdoutfile << "== ERROR" << std::endl;
			stdoutfile << "== Lattice type " << lattice_type_ << " unknown." << std::endl;
			stdoutfile << "== System  " << system_ << " requires a valid lattice type for computation." << std::endl;
			stdoutfile << " = Known types:" << std::endl;
			stdoutfile << " = 's'     : square lattice" << std::endl;
			stdoutfile << " = 't'     : trigonal lattice" << std::endl;
			stdoutfile << " = Aborting calculation." << std::endl;
			errorcode = 1;
		}
	}
	if (mode_ == "integ") { // allows for shortening. Just avoids a bit of confusion.
		mode_ = "integrate";
	}
	if (mode_ == "eq") { // allows for shortening. Just avoids a bit of confusion.
		mode_ = "equilibrate";
	}
	if ((mode_ != "none") && (mode_ != "integrate") && (mode_ != "integrate_cont") && (mode_ != "samp") && (mode_ != "test") && (mode_ != "equilibrate")){
		stdoutfile << "=  ERROR" << std::endl;
		stdoutfile << "== mode " << mode_ << " unknown." << std::endl;
		errorcode = 1;
	} else if ( mode_ != "none"){
		if ( sqrtN_ * sqrtN_ != N_){
			stdoutfile << "=  ERROR" << std::endl;
			stdoutfile << "== sqrtN * sqrtN = " << sqrtN_ * sqrtN_ << " != " << N_ << " = N." << std::endl;
			errorcode = 1;
		}
		if ((mode_ == "integrate") | (mode_ == "equilibrate")) {
			if ( (eq_mode_ == "Brownian") ||  (eq_mode_ == "brownian") ){
				eq_mode_ = "brownian";
				if ( eq_brownian_kT_omega_ < 0){
					eq_brownian_kT_omega_ = kT_;
				}
				if ( eq_brownian_kT_p_ < 0){
					eq_brownian_kT_p_ = kT_;
				}
				if ( eq_brownian_timestep_ <= 0){
					stdoutfile << "=  ERROR" << std::endl;
					stdoutfile << "== Inappropriate value of eq_brownian_timestep, eq_brownian_timestep = " << eq_brownian_timestep_ << std::endl;
					stdoutfile << "== Value must be positive.\n"
							<< "== Run will be aborted.\n\n";
					errorcode = 1;
				}
			}
			if ( (eq_mode_ == "berendsen") && (tau_berendsen_ <= 0)) {
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Inappropriate value of tau_berendsen, tau = " << tau_berendsen_ << std::endl;
				stdoutfile << "== Value must be positive.\n"
						<< "== Run will be aborted.\n\n";
				errorcode = 1;
			}
			if ( (eq_mode_ == "anneal") && ( (eq_anneal_rate_ <= 0) | (eq_anneal_rate_ >= 1) )) {
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Inappropriate value of anneal_rate, anneal_rate = " << eq_anneal_rate_ << std::endl;
				stdoutfile << "== Value must be between 0 and 1.\n"
						<< "== Run will be aborted.\n\n";
				errorcode = 1;
			}
			if ( (eq_mode_ == "anneal") && (eq_anneal_step_ <= 0) ) {
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Inappropriate value of anneal_step, anneal_step = " << eq_anneal_step_ << std::endl;
				stdoutfile << "== Value must be positive.\n"
						<< "== Run will be aborted.\n\n";
				errorcode = 1;
			}
			if ( (eq_mode_ != "berendsen") && (eq_mode_ != "anneal") && (eq_mode_ != "brownian")) {
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Unknown equilibration mode " << eq_mode_ << std::endl;
				stdoutfile << "== Possible values for eq_mode are\n"
						<<	"== 'berendsen'    (Berendsen thermostat)\n"
						<<	"== 'anneal'       (annealing)\n"
						<<	"== 'brownian'     (Brownian dynamics)\n"
						<<	"== Run will be aborted\n\n";
				errorcode = 1;
			}
			if ( (eq_integrator_type_ != "rk4") && (eq_integrator_type_ != "lf")
					&& (eq_integrator_type_ != "mc") && (eq_integrator_type_ != "vm")
					&& (eq_integrator_type_ != "langevin")){
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Unknown equilibration integrator type " << eq_integrator_type_ << std::endl;
				stdoutfile << "== Possible values are\n"
						<<	"== 'rk4'		(fourth order Runge-Kutta)\n"
						<<	"== 'lf'		(leapfrog)\n"
						<<	"== 'langevin'	(Langevin dynamics integration)\n"
						<<	"== 'vm'		(Vicsek model with angular noise)\n"
						<<	"== 'mc'		(Monte Carlo, Metropolis-Hastings. NOT IMPLEMENTED)\n"
						<<	"== Run will be aborted\n\n";
				errorcode = 1;
			}

		} else if ((mode_ == "integrate") | (mode_ == "integrate_cont")) {
			if (samplestep_ < 0){
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== samplestep value " << samplestep_ << " negative." << std::endl;
				errorcode = 1;
			}
			if ( (sample_integrator_type_ != "rk4") && (sample_integrator_type_ != "lf") && (sample_integrator_type_ != "mc")
					&& (sample_integrator_type_ != "nh") && (sample_integrator_type_ != "np")
					&& (sample_integrator_type_ != "vm") && (sample_integrator_type_ != "langevin")){
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Unknown sampling integrator type " << eq_integrator_type_ << std::endl;
				stdoutfile << "== Possible values are\n"
						<<	"== 'rk4'		(fourth order Runge-Kutta)\n"
						<<	"== 'lf'		(leapfrog)\n"
						<<	"== 'langevin'	(Langevin dynamics integration)\n"
						<<	"== 'mc'		(Monte Carlo, Metropolis-Hastings. NOT IMPLEMENTED)\n"
						<<	"== 'nh'		(Nose-Hoover)\n"
						<<	"== 'np'		(Nose-Poincare)\n"
						<<	"== 'vm'		(Vicsek model with angular noise)\n"
						<<	"== Run will be aborted\n\n";
				errorcode = 1;
			}
			if ( (ensemble_ == "nvt") && ((sample_integrator_type_ == "lf") || (sample_integrator_type_ == "rk4"))){
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== The integrator " << sample_integrator_type_
						<< " is incompatible with the nvt ensemble."<< std::endl;
				stdoutfile << "== Possible values are\n"
						<<	"== 'nh'		(Nose-Hoover)\n"
						<<	"== 'np'		(Nose-Poincare)\n"
						<<	"== 'mc'		(Monte Carlo, Metropolis-Hastings. NOT IMPLEMENTED)\n"
						<<	"== Run will be aborted\n\n";
				errorcode = 1;
			}
			if ( (ensemble_ == "nve") && ((sample_integrator_type_ == "nh") || (sample_integrator_type_ == "np") || (sample_integrator_type_ == "langevin"))){
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== The integrator " << sample_integrator_type_
						<< " is incompatible with the nve ensemble."<< std::endl;
				stdoutfile << "== Possible values are\n"
						<<	"== 'rk4'		(fourth order Runge-Kutta)\n"
						<<	"== 'lf'		(leapfrog)\n"
						<<	"== Run will be aborted\n\n";
				errorcode = 1;
			}
			if ( (ensemble_ == "Brownian") ||  (ensemble_ == "brownian") ){
				ensemble_ = "brownian";
				if ( brownian_kT_omega_ < 0) {
					brownian_kT_omega_ = kT_;
				}
				if ( brownian_kT_p_ < 0) {
					brownian_kT_p_ = kT_;
				}
			}
			if ( (sampling_time_sequence_ != "lin") && (sampling_time_sequence_ != "log")){
				stdoutfile << "=  ERROR" << std::endl;
				stdoutfile << "== Unknown sampling time sequence " << sampling_time_sequence_ << std::endl;
				stdoutfile << "== Possible values are\n"
						<<	"== 'lin'		(linear times sequence)\n"
						<<	"== 'log'		(logarithmic time sequence)\n"
						<<	"== Run will be aborted\n\n";
				errorcode = 1;
			}
		}
	}
	if (isnan(qmax_) || qmax_ < 0){
		stdoutfile 	<< "=  WARNING" << std::endl;
		stdoutfile 	<< "== qmax: value " << qmax_ << " not feasible.\n"
					<< "== qmax: choosing 2 * Pi instead." << std::endl;
		qmax_ = 2 * M_PI;
	}
	if ( (qbin_type_ != "all") && (qbin_type_ != "mult")){
		stdoutfile << "=  ERROR" << std::endl;
		stdoutfile << "== Unknown qbin_type, qbin_type = " << qbin_type_ << std::endl;
		stdoutfile << "== Possible values are\n"
				<<	  "== 'all'         (all possible values of q)\n"
				<<	  "== 'mult'        (only integer multiples of 2pi/L_x and 2pi/L_y)\n"
				<<    "== Run will be aborted.\n\n";
		errorcode = 1;
	}
	if ( gamma_ld_p_ < 0) {
			stdoutfile << "=  ERROR" << std::endl;
			stdoutfile << "== gamma_ld_p = " << gamma_ld_p_ << " is negative.\n"
					<<    "== Please set to 0 if not required, otherwise choose positive value.\n"
					<<    "== Run will be aborted.\n\n";
			errorcode = 1;
		}
	if ( gamma_ld_om_ < 0) {
		stdoutfile << "=  ERROR" << std::endl;
		stdoutfile << "== gamma_ld_om = " << gamma_ld_om_ << " is negative.\n"
					<<    "== Please set to 0 if not required, otherwise choose positive value.\n"
					<<    "== Run will be aborted.\n\n";
		errorcode = 1;
	}
	// TODO
	//struct stat st;
	  //      stat(output_folder_.c_str(), &st);
	//if(st.st_mode & S_IFDIR){
  //		mkdir(output_folder_.c_str());//added by david
	//}


	return errorcode;
}


void parameters::initialize_bins(){

	initialize_rbin(N_rbin_);
	initialize_qbin();
	//	double r_max = .5 * L_, r_min = 1;
	//	int Nlin_r = 0;
//	double curpos = r_min;
//	// Ensure that the the logarithmic bins (in space) are not smaller than the minimal bin size min_bindwith_r_.
//		double exp_factor = std::exp(1./ (N_rbin_ - Nlin_r) * std::log(r_max / curpos));
//		while ( ( (exp_factor - 1) * curpos < min_binwidth_r_ ) && ( Nlin_r < N_rbin_ ) ){
//			Nlin_r++;
//			curpos += min_binwidth_r_;
//			exp_factor = std::exp(1./ (N_rbin_ - Nlin_r) * std::log(r_max / curpos));
//
//		}

}

void parameters::initialize_rbin(int N_rbin){
//	rbin_ = lin_log_bin(1, .5 * L_, min_binwidth_r_, Nlin, Nlog);
/*	if ( system_ == "xy" || system_ == "fvm" ) {
		rbin_ = log_bin(1, .5 * L_, min_binwidth_r_, N_rbin);
	} else if ( system_ == "mxy" || system_ == "vm" || (system_ == "fmxy") ) {
		rbin_ = log_bin(0, .5 * L_, min_binwidth_r_, N_rbin);
	}
	N_rbin_ = rbin_.size(); // Binning function may be off by 1 in size, stored here.
	*/
	rbin_.clear();
	double bin_r_max = .5 * L_;
	if ( system_ == "xy" || system_ == "fvm" ) {
		if (lattice_type_ == 't'){
			bin_r_max = .5 * .5 * sqrt(3) * L_; // System size in Y direction different from the one in X direction
		} else if (lattice_type_ == 's'){
			bin_r_max = .5 * L_;
		}
	} else if ( system_ == "mxy" || system_ == "vm" || (system_ == "fmxy") ) {
		bin_r_max = .5 * L_;
	}
	double bin_step = bin_r_max / N_rbin;
	for (int i = 0; i < N_rbin; i++){
		rbin_.push_back((i+1) * bin_step);
	}
}


void parameters::initialize_qbin(){
	qbin_.clear();
	topology::Vector2d linwidth_vec;
	if ( ( system_ == "xy" || system_ == "fvm") && lattice_type_ == 't') {
		linwidth_vec = topology::Vector2d(2 * M_PI / L_, 2 * M_PI / (sqrt(3) /2 * L_ ) ) ;
	} else if ( ( system_ == "mxy" || system_ == "vm" || (system_ == "fmxy") ) || lattice_type_ == 's') {
		linwidth_vec = topology::Vector2d(2 * M_PI / L_, 2 * M_PI / L_ ) ;
	}
	std::vector<double> distances;
	if (qbin_type_ == "all") {
		for (int i = 0; i * linwidth_vec.get_x() < qmax_ ; i++){
			for (int j = 0; j * linwidth_vec.get_y() < qmax_ ; j++){
				if ( i + j != 0){
					distances.push_back(std::sqrt(norm2(topology::Vector2d( i * linwidth_vec.get_x(), j * linwidth_vec.get_y()))));
				}
			}
		}
	} else if (qbin_type_ == "mult") {
		for (int i = 1; i * linwidth_vec.get_x() < qmax_ ; i++){
			distances.push_back(i * linwidth_vec.get_x());
		}
		if (linwidth_vec.get_x() != linwidth_vec.get_y() ) {
			for (int i = 1; i * linwidth_vec.get_y() < qmax_ ; i++){
				distances.push_back(i * linwidth_vec.get_y());
			}
			std::sort(distances.begin(),distances.end());
		}
	}
	std::sort(distances.begin(),distances.end()); // Sorts the vector
	// Next step removes multiple mentions of the same distance.
	// std::unique has wrong output, so this is manually
	std::vector<double> distances_unique;
	distances_unique.push_back(distances[0]);
	for (size_t i = 1; i < distances.size(); i++){
		if (distances[i] != distances_unique.back()){
			distances_unique.push_back(distances[i]);
		}
	}
	if (qbin_type_ == "all") {
		size_t i_last = 0;
		for (size_t i = 0; i < distances_unique.size() ; i++){
			if ( ( ( i == 0) || (distances_unique[i] - distances_unique[i_last] > min_binwidth_q_ * std::min(linwidth_vec.get_x(),linwidth_vec.get_y()) ) )
					&& (distances_unique[i] < qmax_) ){
				qbin_.push_back(1.0001 * distances_unique[i]); // to avoid problems with rounding errors in bin, we multiply by 1.00001.
				i_last = i;
			}
		}
		// That one is problematic, should be removed (but it is part of older implementations)
		//if (qbin_.back() < qmax_){
		//	qbin_.push_back(qmax_);
		//}
	} else if (qbin_type_ == "mult") {
		for (size_t i = 0; i < distances_unique.size() ; i++){
			qbin_.push_back(distances[i]);
		}
	}

}


void parameters::scale_tau(double scale_factor){
	tau_berendsen_ *= scale_factor;
}
