/*! \file routines.cpp
 *
 *  \brief Cpp-File to the namespace routines. Implements the functions from
 *  routines.h
 *
 *  Detailed descriptions can be found in routines.h. This file contains some further
 *  comments on the precise way the funcitons are implemented. It is recommended to read
 *  the code here in order to better understand what the functions do.
 (
 *  \date Created: 2019-12-11
 *  \date Last Updated: 2023-08-06
 *
  *  \author Thomas Bissinger
  *
 */

#include "routines.h"
// =====================================================================
// =====================================================================
// NAMESPACE ROUTINES
// =====================================================================
// =====================================================================


//   ===================================================================================================
//   integration routine
//   ===================================================================================================
int routines::integration(parameters par){
	std::ifstream samp_infile("input/samp_config.in", std::ios::in);
	std::ifstream eq_samp_infile("input/eq_samp_config.in", std::ios::in);
	std::ofstream stdoutfile(par.outfilename(), std::ios::out | std::fstream::app);
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	initprint("Initialization",stdoutfile);
	stdoutfile << " = " << par.N() << " Particles in a periodic square box of length " << par.L() << std::endl << std::endl;
	std::string eq_sampout_name = par.output_folder() + "/eq_sampling_output" + par.job_id() + ".m";//added by David
	std::ofstream eq_sampoutfile;


	if (par.randomseed() == -1) {
		srand(time(NULL));
	} else {
		srand(par.randomseed());
	}
	//   ===================================================================================================
	//   Initialization
	//   ===================================================================================================
	//   1. group initialization
	group G(par);
	G.initialize(par);


	//   2. integrator initialization
	integrator integ(par.dt(),par.eq_integrator_type());
	stdoutfile << " = Integrator initialized." << std::endl << std::endl;



	//   3. sampler(s) initialization
	sampler samp(par);
	sampler eq_samp(par);
	samp.switches_from_file(samp_infile,stdoutfile);
	if ( par.eq_sampswitch() ){
		eq_samp.switches_from_file(eq_samp_infile,stdoutfile);
	} else {
		eq_samp.all_switches_off();
	}
	if ( ! std::isnan(par.eq_Tprintstep())){
		eq_samp.print_snapshots_on();
	}

	//   ===================================================================================================
	//   Equilibration
	//   ===================================================================================================
	double av_time = par.eq_av_time();
	double agreement_threshold = par.eq_agreement_threshold();
	int eq_returnval = 0;
	std::string break_cond = par.eq_breakcond();


	if ((par.mode() == "integrate") || (par.mode() == "equilibrate")) {
		double tau_berendsen = par.tau_berendsen();
		double eq_Tmax = par.eq_Tmax();
		double current_temperaturediff = 0;
		double current_T = 0;
		double t = 0;
		while (eq_returnval != 1){
			std::ofstream eqfile;
			eqfile.open(par.output_folder() + "/equilibration_final" + par.job_id() + ".out", std::ios::out );//added by David
			eq_returnval = routines::equilibrate(G, par, eq_samp,
									eq_Tmax, t, break_cond, stdoutfile);
			G.print_group(eqfile);
			if (eq_returnval != 1){
				if ( par.system() == "xy" || par.system() == "mxy" || par.system() == "fmxy" ){
					current_T = G.calc_temperature();
				}
				current_temperaturediff = std::abs(current_T - par.kT());

				// Calculate new tmax for berendsen or annealing equilibration by assuming exponential decay. Factor 1.5 is to account
				// for deviations from ideal exponential decay. Agreement threshold gives percentage.
				if ( par.eq_mode() == "berendsen" ) {
					eq_Tmax = 1.5 * std::log(current_temperaturediff / (agreement_threshold * par.kT()) ) * tau_berendsen;
				} else if ( par.eq_mode() == "anneal" ) {
					eq_Tmax = 1.5 * std::log(current_temperaturediff / (agreement_threshold * par.kT()) ) * par.eq_anneal_step() / (1 - par.eq_anneal_rate() );
					// Approximately, annealing corresponds to exponential decay exp(-t / tau) with tau = - anneal_step / log(anneal_rate) = anneal_step / (1 - anneal_rate).
				}
				// eq_Tmax = - 2 * std::log(current_T / (current_T + agreement_threshold) ) * tau_berendsen; // Factor 2 is to keep
				// TODO Proper thermostat repeat conditions.

				if ( eq_Tmax > .5 * par.eq_Tmax() ){
					eq_Tmax = .5 * par.eq_Tmax();
				} else if (	eq_Tmax < 20 * av_time ) {
					if ( eq_Tmax > 3 * av_time ) {
						eq_Tmax = 20 * av_time; 						// There should still be enough averaging possible in the equilibration time interval
					} else {
						eq_returnval = 1;								// If the current value is close enough to the desired one, we can count the sim as equilibrated.
					}
				}
				if ( eq_returnval != 1){
					stdoutfile << "++ Desired temperature not reached. " << std::endl;
					if ( break_cond != "time_hard" ){
						stdoutfile << "++ Continuing equilibration, using new Tmax = " << eq_Tmax << std::endl << std::endl;
					} else {
						stdoutfile << "++ Terminating equilibration because eq_breakcond = 'time_hard' " << std::endl << std::endl;
						eq_returnval = 1;
					}
				} else {
					stdoutfile << "++ Temperature close enough to desired value. Assuming equilibration." << std::endl << std::endl;
				}
			}
		}
	}
	if ( par.eq_sampswitch() ){
		eq_sampoutfile.open(eq_sampout_name, std::ios::out );
		eq_sampoutfile << "N = " << par.N() << ";\n"
					"Tmax = " << par.eq_Tmax() << ";" << std::endl;
		eq_samp.print_matlab(eq_sampoutfile);
		eq_samp.print_averages_matlab(eq_sampoutfile);
		eq_sampoutfile.close();
	}

	//   ===================================================================================================
	//   Integration
	//   ===================================================================================================
	if ( (par.mode() == "integrate") | (par.mode() == "integrate_cont" )){
		routines::integrate_snapshots(G, par, samp, stdoutfile);
	}
	elapsed = std::chrono::high_resolution_clock::now() - start;
	std::cout << "Run complete! Time consumed: " << elapsed.count() << " s" << std::endl;
	stdoutfile << "Run complete! Time consumed: " << elapsed.count() << " s" << std::endl;
	stdoutfile << "Exit value: " << 1 << std::endl;
	stdoutfile.close();
	return 1;

}


//   ===================================================================================================
//   sampling routine
//   ===================================================================================================
int routines::sampling(parameters par){
	std::ifstream samp_infile("input/samp_config.in", std::ios::in);
	std::ofstream sampoutfile;
	std::ofstream stdoutfile(par.outfilename(), std::ios::out | std::fstream::app);
	std::ifstream snaplist(par.snap_overview_file(), std::ios::in);


	stdoutfile << " + Reading snapshots from file: " << par.snap_overview_file() << std::endl;
	std::cout << " + Reading snapshots from file: " << par.snap_overview_file() << std::endl;

	std::string sampout_name = par.output_folder() + "/sampling_output" + par.job_id() + ".m";//added by David



	stdoutfile << "  + Writing sampling data to " << sampout_name << std::endl;
	stdoutfile << std::endl;


	group G(par), G_initial;
	G.initialize_random(par.kT());

	sampler samp(par);
	samp.switches_from_file(samp_infile,stdoutfile);

	auto start_total = std::chrono::high_resolution_clock::now();
	auto start_sampling = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_sampling = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_total;

	std::string snapshotname;
	double t, t_0 = 0;
	int sampcount = 0;
	while(snaplist >> t >> snapshotname) {
		G.read_from_snapshot(snapshotname);
		if ( par.system() == "mxy" || par.system() == "vm" || par.system() == "fmxy" ) {
			G.fill_partition();
			if ( par.system() == "fmxy" ){
				G.generate_neighbor_list();
			}
		}
		if ( sampcount == 0){
			G_initial = G;
			t_0 = t;
		}

		start_sampling = std::chrono::high_resolution_clock::now();
		samp.sample(G,G_initial,t - t_0);
		elapsed_sampling += std::chrono::high_resolution_clock::now() - start_sampling;
		sampcount++;

		stdoutfile << "    + " << sampcount << ": Sampling performed at t = " << t << std::endl;
		std::cout << "    + " << sampcount << ": Sampling performed at t = " << t << std::endl;
	}

	sampoutfile.open(sampout_name, std::ios::out );
	sampoutfile << "N = " << par.N() << ";\n"
			"Tmax = " << par.Tmax() << ";" << std::endl;
	samp.print_matlab(sampoutfile);
	samp.print_averages_matlab(sampoutfile);



	elapsed_total = std::chrono::high_resolution_clock::now() - start_total;
	stdoutfile << "  + Time consumption diagnostics for sampling:\n";
	stdoutfile << "  + Total time: " << elapsed_total.count() << " s" << std::endl;
	stdoutfile << "  + Sampling time: " << elapsed_sampling.count() << " s" << std::endl;

	std::cout << "Sampling complete! Time consumed: " << elapsed_total.count() << " s" << std::endl;
	stdoutfile << "Exit value: " << 1 << std::endl;
	sampoutfile.close();
	stdoutfile.close();


	return 1;
}

















//   ===================================================================================================
//   equilibration routine
//   ===================================================================================================
int routines::equilibrate(group& G, const parameters& par, sampler& samp,
		const double Tmax, double& t,
		const std::string breakcond, std::ofstream& stdoutfile){
	std::string routine_name = "equilibrate";
	initprint(routine_name,stdoutfile);
	std::string snapfile_basename = par.output_folder()+"/eq_snapshot";
	std::string snapout_name = snapfile_basename + "_overview.out";
	std::string snapfile_name;
	std::ofstream snapoutfile;
	snapoutfile.open(snapout_name, std::ios::out);

	initprint(routine_name,stdoutfile);

	int exitval = 0;


	double agreement_threshold = par.eq_agreement_threshold();
	stdoutfile		<< "  ++ eq_mode:                '" << par.eq_mode() << "'\n";
	stdoutfile		<< "  ++ Tmax:                   " << Tmax << "\n";
	stdoutfile		<< "  ++ Start time:             " << t << "\n";
	if ( par.eq_mode() == "berendsen" ){
		stdoutfile 	<< "  ++ tau_berendsen:          " << par.tau_berendsen() << "\n"
					<< "  ++ eq_av_time:             " << par.eq_av_time() << "\n";
	} else if ( par.eq_mode() == "anneal" ){
		stdoutfile 	<< "  ++ anneal_rate:            " << par.eq_anneal_rate() << "\n"
					<< "  ++ anneal_step:            " << par.eq_anneal_step() << "\n";
	}
	stdoutfile		<< "  ++ Break condition:        '" << breakcond << "'\n"
					<< "  ++ Agreement threshold:    " << agreement_threshold << "\n"
					<< "  ++ eq_samplingswitch:      " << par.eq_sampswitch() << "\n";
	if( par.eq_sampswitch()){
		stdoutfile  << "  ++ eq_samp_time_sequence:  " << par.eq_samp_time_sequence() << "\n"
					<< "  ++ eq_Nsamp:               " << par.eq_Nsamp() << "\n";

	}

	stdoutfile << std::endl;
	if ((breakcond != "time") && (breakcond != "time_hard") && (breakcond != "temperature") && (breakcond != "any")){
		stdoutfile << "!!! ERROR !!! Unknown break condition '" << breakcond << "'!\n" <<
				"!!! Known break conditions are\n" <<
				"!!! 'temperature'     :  Breaks if equilibrium is reached\n" <<
				"!!! 'time'            :  Breaks if eq_Tmax is reached\n" <<
				"!!! 'time_hard'       :  Breaks if eq_Tmax is reached (and does not perform reruns)\n" <<
				"!!! 'any'             :  Either of the above\n\n" <<
				"!!! ERROR. EXITING !!!" << std::endl;
		return -1;
	}
	if (Tmax <= 0){
		stdoutfile 	<< "!!! ERROR !!! eq_Tmax non-positive!\n"
					<< "!!! ERROR. EXITING.\n\n";
		return -1;
	}


	group G_initial;
	if ( par.on_fly_sampling() ){
		G_initial = G; // This assumes G has already been initialized outside.
	}

	auto start_total = std::chrono::high_resolution_clock::now();
	auto start_integration = std::chrono::high_resolution_clock::now();
	auto start_sampling = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_integration = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_sampling = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_sample_print = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_total;


	double t_0 = t;
	double tstart = t;
	std::vector<double> sampletimes;
	if ( par.eq_samp_time_sequence() == "lin"){
		double samptime = t_0;
		while ( samptime <= t_0 + Tmax ){
			sampletimes.push_back(samptime);
			samptime += par.eq_av_time();
		}
	} if ( par.eq_samp_time_sequence() == "log"){
		sampletimes = log_bin(1, Tmax, par.dt(), par.eq_Nsamp() );
	}

	int sampcount = 0;

	bool breakswitch = 0; // Switch that can be used to terminate the loop.
	std::vector<double> Av_Times, Av_Temperatures;

	integrator integ;
	integ.integrator_eq(par);
	integ.initialize(par,G);

	double t_last = 0, sum_temperature = 0;
	double T_cur;

	int count_av_agreements = 0; // Counts agreement between consecutive averages (used to estimate equilibrium)
	int sufficient_agreement = 3; // Threshold for claiming equilibrium
	double eq_time = -1;

	// For printing snapshots during equilibration (when par.eq_Tprintstep is set)
	std::vector<double> temperatures_for_printing;

	std::vector<double> momdot = {};
	double t_brownian_last = 0;

	while ( !breakswitch){ // Loop ends after break condition is reached.
		                   // Some storage may be expensive, so it is only done when it is needed for printing.

		start_integration = std::chrono::high_resolution_clock::now();
		G = integ.integrate(G,momdot);
		elapsed_integration += std::chrono::high_resolution_clock::now() - start_integration;

		///////////////////////////////////////////////////////////
		// Annealing or Berendsen for Hamiltonian Models
		///////////////////////////////////////////////////////////
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" ) {
			if (par.eq_mode() == "berendsen" ) {
				T_cur = integ.berendsen_thermostat(G,par.kT(),par.tau_berendsen());
				sum_temperature += T_cur;
			}
			if ( (par.eq_mode() == "brownian") && (t - t_brownian_last > par.brownian_timestep() )){
				G.set_temperature_w(par.eq_brownian_kT_omega());
				G.set_temperature_p(par.eq_brownian_kT_p());
				G.mom_to_zero();
				t_brownian_last = t;
			}
			if ( ( ( par.eq_mode() == "berendsen" ) && (t - t_last >= par.eq_av_time() ) )
					|| ( ( par.eq_mode() == "anneal" ) && (t - t_last >= par.eq_anneal_step() ) )
					){
				if ( par.eq_mode() == "berendsen" ) {
					Av_Temperatures.push_back(sum_temperature * par.dt() / (t - t_last)); // Averages temperature in vector.
					sum_temperature = 0; // Reinitializing
				} else if ( par.eq_mode() == "anneal" ) {
					Av_Temperatures.push_back(G.calc_temperature()); // No averaging
					if (Av_Temperatures.back() > par.kT()) {
						G.scale_mom(par.eq_anneal_rate());
					} else if ( (std::abs(Av_Temperatures.back()/par.kT() - 1) > agreement_threshold) | (agreement_threshold > (1 - par.eq_anneal_rate()) ) ){
						// Also the possibility of annealing up near the desired temperature, but only if
						// a) the current temperature is not within the agreement region or
						// b) the annealing rate is smaller than the agreement threshold.
						G.scale_mom(1.0 / par.eq_anneal_rate());
					}
				}
				if ( ! par.eq_sampswitch() ){
					std::cout << "  eq " << t << " " << Av_Temperatures.back() << std::endl;
				}

				if ( Av_Temperatures.size() > 3 ){ // Avoids too long vectors
					Av_Temperatures.erase(Av_Temperatures.begin());
				}
				t_last = t;
				if ( std::abs(Av_Temperatures.back()/par.kT() - 1) < agreement_threshold){
					count_av_agreements++;
					if ( count_av_agreements == sufficient_agreement){
						if (eq_time == -1 ){
							eq_time = t;
						}
						if ( (breakcond == "temperature") || (breakcond == "any")){ // Changes
							breakswitch = 1;
						}
					}
				} else {
					count_av_agreements = 0;
				}
			}
		}
		if ( par.eq_sampswitch() && sampcount < (int)sampletimes.size() && t >= sampletimes[sampcount]){ // comparison important for calculations that don't terminate at Tmax. Further sampling
			start_sampling = std::chrono::high_resolution_clock::now();
			samp.sample(G,G_initial,t); // Must be changed if one wants to calculate correct ACF (G_initial would have to be supplied externally). But that's not that important in this setting. Just FYI.
			elapsed_sampling += std::chrono::high_resolution_clock::now() - start_sampling;
			if ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" ) {
				stdoutfile << "    ++ eq ++ " << sampcount << " t  " << t << "  T  " << samp.get_last_temperature() << std::endl;
				std::cout  << "    ++ eq ++ " << sampcount << " t  " << t << "  T  " << samp.get_last_temperature() << std::endl;
			} else {
				stdoutfile << "    ++ eq ++ " << sampcount << " t  " << t << std::endl;
				std::cout  << "    ++ eq ++ " << sampcount << " t  " << t << std::endl;
			}
			sampcount++;
		}
		// Possibility to print intermediate snapshots
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" ) {
			if ( ! std::isnan(par.eq_Tprintstep()) && (Av_Temperatures.size() > 1) ){
				if ( temperatures_for_printing.size() == 0 ){
					for (double T = par.kT(); T < Av_Temperatures.back(); T += par.eq_Tprintstep()){
						temperatures_for_printing.push_back(T);
					}
				} else if ( (Av_Temperatures.back() < temperatures_for_printing.back()) && (Av_Temperatures.back() > par.kT())) { // Only prints snapshots above the desired temperature
					snapfile_name = snapfile_basename + "_T_" + std::to_string(temperatures_for_printing.back()) + ".out";
					samp.sample_snapshots(G,temperatures_for_printing.back(),snapfile_name,snapoutfile);
					temperatures_for_printing.pop_back();
				}
			}
		}

		t += par.dt();
		if ( ( G.get_group_type() == "vm" || G.get_group_type() == "fvm" || breakcond == "time" || breakcond == "time_hard" || breakcond == "any" ) && (t > t_0 + Tmax) ){
				breakswitch = 1; // Breaks if break condition is time
		}
	}

	if (eq_time != -1){
		stdoutfile << "  O Equilibration reached at t = " << eq_time << std::endl;
	} else {
		stdoutfile << "  O Equilibration not reached after t = " << t << std::endl;
	}
	terminateprint(routine_name,stdoutfile);
	snapoutfile.close();
	if (eq_time != -1){
		exitval = 1;
	} else if ( G.get_group_type() == "vm" || G.get_group_type() == "fvm") {
		exitval = 1;
	}


	elapsed_total = std::chrono::high_resolution_clock::now() - start_total;
	stdoutfile << "  + Time consumption diagnostics for " << routine_name << ":\n";
	stdoutfile << "  + Total time: " << elapsed_total.count() << " s" << std::endl;
	stdoutfile << "  + Integration time: " << elapsed_integration.count() << " s" << std::endl;
	stdoutfile << "  + Sampling time: " << elapsed_sampling.count() << " s" << std::endl;
	stdoutfile << "  + Total time per particle per system second: " << elapsed_total.count()/(par.N() * (t - tstart)) << " s" << std::endl;
	return exitval;
}








/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sampling routines
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void routines::integrate_snapshots(group& G, const parameters& par, sampler& samp,
		std::ofstream& stdoutfile){
	std::string routine_name = "sample_nvt";
	std::string sampout_name = par.output_folder() + "/sampling_output" + par.job_id() + ".m";//added by David
	std::string posout_name = par.output_folder() + "/positions_output.out";//added by David
	std::string snapfile_basename = par.output_folder() + "/snapshot" + par.job_id();//added by David
	std::string snapout_name = par.snap_overview_file();
	std::string snapfile_name;
	initprint(routine_name,stdoutfile);
	std::ofstream sampoutfile;
	std::ofstream snapoutfile;
	std::ofstream posoutfile;
	if( par.on_fly_sampling()){
		stdoutfile << "  + Writing sampling data to " << sampout_name << std::endl;

	}
	if ( par.print_snapshots() ){
		if ( par.mode() == "integrate_cont" ){
			snapoutfile.open(snapout_name, std::ios_base::app);
		} else {
			snapoutfile.open(snapout_name, std::ios::out);
			posoutfile.open(posout_name, std::ios::out);
		}
		stdoutfile << "  + Writing snapshots to " << snapout_name << std::endl;
		stdoutfile << "  + Writing positions to " << posout_name << std::endl;
	}
	stdoutfile << std::endl;

	group G_initial;
	if ( par.on_fly_sampling() ){
		G_initial = G; // This assumes G has already been initialized outside.
	}

	auto start_total = std::chrono::high_resolution_clock::now();
	auto start_integration = std::chrono::high_resolution_clock::now();
	auto start_sampling = std::chrono::high_resolution_clock::now();
	auto start_sample_print = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_integration = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_sampling = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_sample_print = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_total;
	double dt = par.dt();


	double t = 0;
	double tstart = 0;
	double t_0 = - par.av_time_spacing();
	std::vector<double> sampletimes;
	if ( par.sampling_time_sequence() == "lin"){
		double samptime = par.samplestart();
		while ( samptime <= par.Tmax() ){
			sampletimes.push_back(samptime);
			samptime += par.samplestep();
		}
	} if ( par.sampling_time_sequence() == "log"){
		sampletimes = log_bin(par.samplestart(), par.Tmax(), par.dt(), par.Nsamp() );
	}

	int sampcount = 0;
	int snapcount = 0;
	if ( par.mode() == "integrate_cont"){
		std::ifstream snaplist(par.snap_overview_file(), std::ios::in);
		stdoutfile << " + Mode is '" << par.mode() << "'." << std::endl;
		stdoutfile << " + Reading snapshots from file: " << par.snap_overview_file() << std::endl;
		std::string snapshotname;
		while(snaplist >> tstart >> snapshotname) {
			snapcount++;
		}
		if (snapshotname.rfind(par.output_folder() + "/", 0) != 0){
			snapshotname=par.output_folder() + "/" + snapshotname;
		}
		G.read_from_snapshot(snapshotname);
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "vm" ){
			G.fill_partition();
		}
		tstart += dt;
		stdoutfile << " + Group initialized from " << snapshotname << "." << std::endl;
		stdoutfile << " + snapcount is  " << snapcount <<"." << std::endl;
		stdoutfile << " + Simulation is continued at time " << tstart <<"." << std::endl;
		if ( snapcount >= (int)sampletimes.size()){
			stdoutfile << "\n ! WARNING!" << std::endl;
			stdoutfile << " ! Previous runs has more snapshots than specified in input!\n"
					<< " ! No computation will be performed!" << std::endl << std::endl;
		}
		if ( tstart > par.Tmax() ){
			stdoutfile << "\n ! WARNING!" << std::endl;
			stdoutfile << " ! Timestamp of last snapshots exceeds Tmax specified in input!\n"
					<< " ! No computation will be performed!" << std::endl << std::endl;
		}

	} else {
		if (par.system() == "xy"){
			start_sample_print = std::chrono::high_resolution_clock::now();
			G.print_r(posoutfile);
			elapsed_sample_print += std::chrono::high_resolution_clock::now() - start_sample_print;
		}
	}

	std::vector<double> momdot = {}; // required for leapfrog


	integrator integ;
	integ.integrator_sample(par);
	integ.initialize(par,G);

	t=tstart;

	group G_MSD_last = G;

	////////////////////////////////////////////////////////////////////////////////////
	// Need to calculate refresh rate of MSD to make sure periodic boundaries are handled properly.
	////////////////////////////////////////////////////////////////////////////////////
	std::vector<double> MSD_vec;
	double t_MSD_refresh = 0, t_MSD_last = 0;
	double t_brownian_last = 0;
	double maxvel = 1; // some initialization, just in case.
	if (par.system() == "xy" || par.system() == "mxy" ){
		std::vector<double> velocity_vec = G.time_derivative_coord();
		maxvel = std::max(*std::max(begin(velocity_vec),end(velocity_vec)), - *std::min(begin(velocity_vec),end(velocity_vec)));
	} else if (par.system() == "vm" ){
		maxvel = par.vm_v();
	} else if (par.system() == "fvm" ){
		maxvel = 1;
	}
	if (par.system() == "xy" || par.system() == "fvm" || par.system() == "fmxy"){
		t_MSD_refresh = M_PI / maxvel;
		MSD_vec = std::vector<double>(par.N(),0.0);
	} else if (par.system() == "mxy" || par.system() == "vm"){
		t_MSD_refresh = std::min(M_PI / maxvel, par.L() / maxvel);
		MSD_vec = std::vector<double>(3 * par.N(),0.0);
	}

	while ( ( t <= par.Tmax() ) && ( snapcount < (int)sampletimes.size() ) ){
				// Loop ends at Tmax. Could be that no sampling is performed there (which is desired,
				// because we can then run the code for equilibration purposes).

		////////////////////////////////////////////////////////////////////////////////////////////
		// INTEGRATION
		////////////////////////////////////////////////////////////////////////////////////////////
		start_integration = std::chrono::high_resolution_clock::now();
		G = integ.integrate(G,momdot);
		elapsed_integration += std::chrono::high_resolution_clock::now() - start_integration;
		t += dt;

		////////////////////////////////////////////////////////////////////////////////////////////
		// MSD CALCULATION
		////////////////////////////////////////////////////////////////////////////////////////////
		if (t - t_MSD_last >= t_MSD_refresh && t >= sampletimes[0]){
			G.accumulative_MSD(MSD_vec,G_MSD_last);
			G_MSD_last = G;
			t_MSD_last = t;
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		// BROWNIAN THERMOSTAT
		////////////////////////////////////////////////////////////////////////////////////////////
		if ( (par.ensemble() == "brownian") && (t - t_brownian_last > par.brownian_timestep() )){
			G.set_temperature_p(par.brownian_kT_p());
			G.set_temperature_w(par.brownian_kT_omega());
			G.mom_to_zero();
			t_brownian_last = t;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		// SAMPLING
		////////////////////////////////////////////////////////////////////////////////////////////
		if ( t >= sampletimes[snapcount] ){
			if (snapcount == 0){
				t = sampletimes[0]; // Sets time to first value in sampletimes (in case of incongruencies)
				t_0 = t;
				t_MSD_last = t;
			}
			if ( par.on_fly_sampling() ){
				start_sampling = std::chrono::high_resolution_clock::now();
//				samp.refresh_qvals(par.qfullmax(), par.qsamps_per_bin());
				samp.sample(G,G_initial,t - t_0);
				G.accumulative_MSD(MSD_vec,G_MSD_last);
				G_MSD_last = G;
				t_MSD_last = t;
				samp.sample_MSD(MSD_vec);
				elapsed_sampling += std::chrono::high_resolution_clock::now() - start_sampling;
				stdoutfile << "    ++ " << sampcount << ": Sampling performed at t = " << t << std::endl;
				std::cout << "    ++ " << sampcount << ": Sampling performed at t = " << t << std::endl;
				sampcount++;
			}
			snapfile_name = snapfile_basename + "_" + std::to_string(snapcount) + ".out";
			start_sample_print = std::chrono::high_resolution_clock::now();
			samp.sample_snapshots(G,t,snapfile_name,snapoutfile);
			elapsed_sample_print += std::chrono::high_resolution_clock::now() - start_sample_print;
			if ( par.print_snapshots() ){
				stdoutfile << "    + " << snapcount << ": Snapshot stored at t = " << t << std::endl;
				std::cout << "    + " << snapcount << ": Snapshot stored at t = " << t << std::endl;
			}
			snapcount++;
		}
	}

	start_sample_print = std::chrono::high_resolution_clock::now();
	if ( par.on_fly_sampling() ){
		sampoutfile.open(sampout_name, std::ios::out );
		sampoutfile << "N = " << G.get_N() << ";\n"
					"Tmax = " << par.Tmax() << ";" << std::endl;
		samp.print_matlab(sampoutfile);
		samp.print_averages_matlab(sampoutfile);
		sampoutfile.close();
	}
	elapsed_sample_print += std::chrono::high_resolution_clock::now() - start_sample_print;

	// At least the final geometry should be printed
	if ( ! par.print_snapshots() ){
		snapfile_name = snapfile_basename + "_final.out";
		std::ofstream final_snapfile(snapfile_name, std::ios::out );
		//final_snapfile.open(snapfile_name, std::ios::out );
		G.print_group(final_snapfile);
		stdoutfile << "    +  " << ": Final snapshot printed at t = " << t << std::endl;
		std::cout  << "    +  " << ": Final snapshot printed at t = " << t << std::endl;
	}
	if (snapoutfile.is_open()){
		snapoutfile.close();
	}
	elapsed_total = std::chrono::high_resolution_clock::now() - start_total;
	stdoutfile << "  + Time consumption diagnostics for " << routine_name << ":\n";
	stdoutfile << "  + Total time: " << elapsed_total.count() << " s" << std::endl;
	stdoutfile << "  + Integration time: " << elapsed_integration.count() << " s" << std::endl;
	stdoutfile << "  + Sampling time: " << elapsed_sampling.count() << " s" << std::endl;
	stdoutfile << "  + Printing time: " << elapsed_sample_print.count() << " s" << std::endl;
	stdoutfile << "  + Total time per particle per system second: " << elapsed_total.count()/(par.N() * (t - tstart)) << " s" << std::endl;

	terminateprint(routine_name,stdoutfile);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Internal routines (input/output etc)
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void routines::initprint(std::string routine_name, std::ofstream& outfile){
	outfile << "/////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
	outfile << 	"++ Routine " << routine_name << ":\n"
						<< "++ initialized \n";
}

void routines::terminateprint(std::string routine_name, std::ofstream& outfile){
	outfile << 	"++ Routine " << routine_name << ":\n"
							<< "++ terminated \n";
	outfile << "/////////////////////////////////////////////////////////////////////////////////////////////////////////\n\n";
}
