/*! \file sampler.cpp
 *
 *  \brief Cpp-file to class declaration of sampler. Implements the routines
 *  declared insampler.h.
 *
 *  \author Thomas Bissinger

 *  \date Created: 2020-01-20
 *  \date Last Updated: 2023-08-06
 *
 *
 */
#include "sampler.h"



////////////////////////////////////////////////////////////////////////////////////
// Initialization
////////////////////////////////////////////////////////////////////////////////////

void sampler::set_parameters(parameters par){
	par_ = par;
}


void sampler::set_qvals(std::vector<topology::Vector2d> qvals){
	qvals_ = qvals;
}


void sampler::switches_from_vector(const std::vector<bool>& boolvec){

	int i = 0;
	store_static_ = boolvec[i]; i++;

	store_vortices_ = boolvec[i]; i++;


	store_mxq_ = boolvec[i]; i++;
	store_myq_ = boolvec[i]; i++;
	store_eq_ = boolvec[i]; i++;
	store_wq_ = boolvec[i]; i++;
	store_rq_ = boolvec[i]; i++;
	store_jparq_ = boolvec[i]; i++;
	store_jperpq_ = boolvec[i]; i++;
	store_lq_ = boolvec[i]; i++;
	store_teq_ = boolvec[i]; i++;

	store_SCF_ = boolvec[i]; i++;
	store_TCF_ = boolvec[i]; i++;

	store_TransCoeff_ = boolvec[i]; i++;
	store_MemoryKernels_ = boolvec[i]; i++;

	print_snapshots_ = boolvec[i]; i++;

	if (par_.system() == "xy" ){
		store_rq_ = 0;
		store_jparq_ = 0;
		store_jperpq_ = 0;
		store_lq_ = 0;
		store_TransCoeff_ = 0;
		store_MemoryKernels_ = 0;
	} else if (par_.system() == "mxy" ){
		store_vortices_ = 0;
	} else if (par_.system() == "fmxy" ){
		store_vortices_ = 0;
		store_rq_ = 0;
		store_jparq_ = 0;
		store_jperpq_ = 0;
		store_lq_ = 0;
		store_TransCoeff_ = 0;
		store_MemoryKernels_ = 0;
	} else if (par_.system() == "vm" || par_.system() == "fvm"){
		if (par_.system() == "vm") {
			store_vortices_ = 0;
		}
		store_wq_ = 0;
		store_eq_ = 0;
		store_rq_ = 0;
		store_jparq_ = 0;
		store_jperpq_ = 0;
		store_lq_ = 0;
		store_TransCoeff_ = 0;
		store_MemoryKernels_ = 0;
	}

}

void sampler::switches_from_file(std::ifstream& infile, std::ofstream& outfile){
	std::string line;
	std::vector<std::string> var_names = {"static", "vortices", //"H", "W", "M", "Theta", "temperature", "vortices",
			"mxq", "myq", "eq", "wq",
			"rq", "jparq", "jperpq", "lq",
			"teq",
			"SCF", "TCF", "TransCoeff", "MemoryKernels", "refresh_q",
			"print_snapshots"};
	std::vector<bool> vars(var_names.size(),0);

	while(getline(infile,line)) {
		for (size_t i = 0; i < var_names.size(); i++){
			if (line.rfind(var_names[i] + " = y", 0) == 0 ){
				vars[i] = 1;
			}
		}
	}
	switches_from_vector(vars);

	if (par_.mode() != "samp"){
		check_on_fly_sampling();
	}
	outfile << "==========================================================\n";
	outfile << "================  Sampling variables =====================\n";
	outfile << "==========================================================\n";
	for (size_t i = 0; i < var_names.size(); i++){
		outfile << std::setw(15) << std::left << var_names[i] << " :  " << vars[i] << std::endl;
	}
	outfile << "==========================================================\n\n";
}


void sampler::all_switches_on(){
	std::vector<bool> boolvec(NumberOfSwitches_,1);
	switches_from_vector(boolvec);
}

void sampler::all_switches_off(){
	std::vector<bool> boolvec(NumberOfSwitches_,0);
	switches_from_vector(boolvec);
}

void sampler::check_on_fly_sampling(){
	if ( !par_.on_fly_sampling() ){
		all_switches_off();
		if (par_.print_snapshots()){
			print_snapshots_ = 1;
		}
	}
}

void sampler::print_snapshots_on(){
	print_snapshots_ = true;
}

void sampler::print_snapshots_off(){
	print_snapshots_ = false;
}

////////////////////////////////////////////////////////////////////////////////////
// get functions
////////////////////////////////////////////////////////////////////////////////////
std::vector<bool> sampler::get_switches() const {
	std::vector<bool> boolvec;
	boolvec.push_back(store_static_);
	boolvec.push_back(store_vortices_);

	boolvec.push_back(store_mxq_);
	boolvec.push_back(store_myq_);
	boolvec.push_back(store_eq_);
	boolvec.push_back(store_wq_);
	boolvec.push_back(store_teq_);

	boolvec.push_back(store_SCF_);
	boolvec.push_back(store_TCF_);
	boolvec.push_back(store_TransCoeff_);
	boolvec.push_back(store_MemoryKernels_);

	boolvec.push_back(refresh_q_);

	boolvec.push_back(print_snapshots_);

	return boolvec;
}



////////////////////////////////////////////////////////////////////////////////////
// binning functions and qvalue handling
////////////////////////////////////////////////////////////////////////////////////

void sampler::refresh_qvals(const topology::Vector2d& boxsize) {
	if (qvals_.empty() || refresh_q_){
		qvals_.clear();
		topology::Vector2d qgridwidth = 2.0 * M_PI * topology::Vector2d(1 / boxsize.get_x(), 1 / boxsize.get_y());
		if (par_.qbin_type() == "all" ){
			int index = 0;
			if (par_.qfullmax() > par_.qbin().back()){
				index = par_.qbin().size() - 1;
			} else {
				while(par_.qbin()[index] < par_.qfullmax() ){
					index++;
				}
			}
			double qmax2 = par_.qbin()[index] * par_.qbin()[index];
			for(double qx = 0; qx * qx < qmax2 ; qx += 2 * M_PI / boxsize.get_x()){
				for(double qy = 0; qx * qx + qy * qy < qmax2 ; qy += 2 * M_PI / boxsize.get_y()){
					if ( ( qx != 0 ) | (qy != 0) ){ // Avoiding the zero vector since the fluctuations are defined to be zero at vector 0.
						qvals_.push_back(topology::Vector2d(qx,qy));
					}
				}
			}
			std::vector<topology::Vector2d> newq;
			for(size_t i = index; i < par_.qbin().size() - 1; i++){
				newq = topology::qvalues_within_radius(par_.qbin()[i], par_.qbin()[i+1], qgridwidth, par_.qsamps_per_bin());
				qvals_.insert(qvals_.end(), newq.begin(),newq.end());
			}
		} else if (par_.qbin_type() == "mult" ){
			for (int i = 1; i * 2 * M_PI / boxsize.get_x() <= par_.qmax(); i++){
				qvals_.push_back(topology::Vector2d(i * 2 * M_PI / boxsize.get_x(),0));
			}
			for (int i = 1; i * 2 * M_PI / boxsize.get_y() <= par_.qmax(); i++){
				qvals_.push_back(topology::Vector2d(0,i * 2 * M_PI / boxsize.get_y()));
			}
		}
	}
}

std::vector<double> sampler::bin_qvals_to_q(std::vector<double> vals) const {
	std::vector<double> binned_vals(par_.qbin().size(),0.0);
	std::vector<int> bincount(par_.qbin().size(),0);
	int index;
	for(size_t i = 0; i < qvals_.size(); i++){
		index = find_bin(par_.qbin(), std::sqrt(qvals_[i].norm2()));
		binned_vals[index] += vals[i];
		bincount[index]++;
	}
	for(size_t j = 0; j < par_.qbin().size(); j++){
		if (bincount[j] != 0){
			binned_vals[j] /= bincount[j];
		}
	}
	return binned_vals;
}


std::vector<std::complex<double> > sampler::bin_qvals_to_q(std::vector<std::complex<double> > vals) const {
	std::vector<std::complex<double> > binned_vals(par_.qbin().size(),0.0);
	std::vector<int> bincount(par_.qbin().size(),0);
	int index;
	for(size_t i = 0; i < qvals_.size(); i++){
		index = find_bin(par_.qbin(), std::sqrt(qvals_[i].norm2()));
		binned_vals[index] += vals[i];
		bincount[index]++;
	}
	for(size_t j = 0; j < par_.qbin().size(); j++){
		if (bincount[j] != 0){
			binned_vals[j] /= bincount[j];
		}
	}
	return binned_vals;
}

////////////////////////////////////////////////////////////////////////////////////
// sampling functions
////////////////////////////////////////////////////////////////////////////////////

void sampler::sample(const group& G, const group& G_initial, double t){
	nsamp_++;
	averaging_times_.push_back(t);
	refresh_qvals(G.get_boxsize());
	sample_static(G,t);
	sample_TCF(G,G_initial,t);
}



void sampler::sample_static(const group& G, double t){
	std::vector<double> store_double;
	std::vector<std::complex<double> > store_complex;
	if (store_static_){
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" ){
			H_.push_back(G.calc_energy()/G.get_N());
			H_2_.push_back(G.sum_e_squared()/G.get_N());
			Hint_2_.push_back(G.sum_eint_squared()/G.get_N());
			Hkin_2_.push_back(G.sum_ekin_squared()/G.get_N());
			W_.push_back(G.sum_w()/G.get_N());
			W_2_.push_back(W_.back() * W_.back());
			temperature_.push_back(G.calc_temperature());
			temperature_squared_.push_back(std::pow(temperature_.back(),2));
			temperature_omega_.push_back(G.calc_temperature_w());
			temperature_omega_squared_.push_back(std::pow(temperature_omega_.back(),2));
			if ( G.get_group_type() == "mxy" ) {
				temperature_p_.push_back(G.calc_temperature_p());
				temperature_p_squared_.push_back(std::pow(temperature_p_.back(),2));
			}
			store_double = G.calc_helicity(1/par_.kT());
			Upsilon_.push_back(store_double[0]);
			H_x_.push_back(store_double[1]);
			H_y_.push_back(store_double[2]);
			I_x_.push_back(store_double[3]);
			I_y_.push_back(store_double[4]);
			I_x_2_.push_back( I_x_.back() * I_x_.back() );
			I_y_2_.push_back( I_y_.back() * I_y_.back() );
		}
		M_.push_back(G.sum_s()/G.get_N());
		M_2_.push_back(topology::norm2(M_.back()));
		M_4_.push_back(M_2_.back() * M_2_.back());
		absM_.push_back(std::sqrt(M_2_.back()));
		M_angle_.push_back(std::atan2(M_.back().get_y(),M_.back().get_x()));
		Theta_.push_back(G.sum_theta()/G.get_N());
		Theta_2_.push_back(std::pow(Theta_.back(),2));
		Theta_4_.push_back(std::pow(Theta_.back(),4));
		Theta_rel_to_M_.push_back(mod(Theta_.back() -  M_angle_.back(),-M_PI,M_PI));
		Theta_rel_to_M_2_.push_back(std::pow(Theta_rel_to_M_.back(),2));
		Theta_rel_to_M_4_.push_back(std::pow(Theta_rel_to_M_.back(),4));
		if ( G.get_group_type() == "mxy" ) {
			P_.push_back(G.sum_p()/G.get_N());
			P_2_.push_back(G.sum_p_squared()/G.get_N());
			P_4_.push_back(G.sum_p_4()/G.get_N());
		}
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "vm" ) {
			coordination_number_.push_back(G.calc_neighbor_mean(0, 0, 0, 0, 0, 0, 0) / G.get_N());
		}

	}
	if ( store_vortices_ && ( G.get_group_type() == "xy" || G.get_group_type() == "fvm" ) ){
		abs_vortices_.push_back(G.calc_vortexdensity_unsigned());
		signed_vortices_.push_back(G.calc_vortexdensity_signed());
	}

	if ( store_mxq_ ){
		mxq_cur_ = G.calc_fieldfluct(qvals_,"mxq");
		store_complex = bin_qvals_to_q(mxq_cur_);
		mxq_.insert(mxq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_myq_ ){
		myq_cur_ = G.calc_fieldfluct(qvals_,"myq");
		store_complex = bin_qvals_to_q(myq_cur_);
		myq_.insert(myq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_mxq_ && store_myq_){
		topology::Vector2d curM;
		if (store_static_){
			curM = M_.back();
		} else {
			curM = G.sum_s()/G.get_N();
		}
		double xcomp = curM.get_x()/std::sqrt(curM.norm2());
		double ycomp = curM.get_y()/std::sqrt(curM.norm2());
		mparq_cur_ = vector_sum(vector_scale(mxq_cur_,xcomp),vector_scale(myq_cur_,ycomp));
		mperpq_cur_ = vector_sum(vector_scale(mxq_cur_,-ycomp),vector_scale(myq_cur_,xcomp));

		store_complex = bin_qvals_to_q(mparq_cur_);
		mparq_.insert(mparq_.end(), store_complex.begin(), store_complex.end());

		store_complex = bin_qvals_to_q(mperpq_cur_);
		mperpq_.insert(mperpq_.end(), store_complex.begin(), store_complex.end());

	}
	if ( store_wq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" )){
		wq_cur_ = G.calc_fieldfluct(qvals_,"wq");
		store_complex = bin_qvals_to_q(wq_cur_);
		wq_.insert(wq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_eq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" )){
		eq_cur_ = G.calc_fieldfluct(qvals_,"eq");
		store_complex = bin_qvals_to_q(eq_cur_);
		eq_.insert(eq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_rq_ ){
		rq_cur_ = G.calc_fieldfluct(qvals_,"rq");
		store_complex = bin_qvals_to_q(rq_cur_);
		rq_.insert(rq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_jparq_ && G.get_group_type() == "mxy" ){
		jparq_cur_ = G.calc_fieldfluct(qvals_,"jparq");
		store_complex = bin_qvals_to_q(jparq_cur_);
		jparq_.insert(jparq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_jperpq_ && par_.system() == "mxy" ){
		jperpq_cur_ = G.calc_fieldfluct(qvals_,"jperpq");
		store_complex = bin_qvals_to_q(jperpq_cur_);
		jperpq_.insert(jperpq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_lq_ && G.get_group_type() == "mxy" ){
		lq_cur_ = G.calc_fieldfluct(qvals_,"lq");
		store_complex = bin_qvals_to_q(lq_cur_);
		lq_.insert(lq_.end(), store_complex.begin(), store_complex.end());
	}
	if ( store_teq_){
		teq_cur_ = G.calc_fieldfluct(qvals_,"teq");
		store_complex = bin_qvals_to_q(teq_cur_);
		teq_.insert(teq_.end(), store_complex.begin(), store_complex.end());
	}

	if ( store_SCF_){
		std::vector<double> SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"S");
		SCF_Spin_.insert(SCF_Spin_.end(), SCF_storage.begin(), SCF_storage.end());

		SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"S_par");
		SCF_Spin_par_.insert(SCF_Spin_par_.end(), SCF_storage.begin(), SCF_storage.end());

		SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"S_perp");
		SCF_Spin_perp_.insert(SCF_Spin_perp_.end(), SCF_storage.begin(), SCF_storage.end());


		SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"anglediff");
		SCF_anglediff_.insert(SCF_anglediff_.end(), SCF_storage.begin(), SCF_storage.end());

		if ( G.get_group_type() == "mxy" || G.get_group_type() == "vm" ) {
			SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"g");
			SCF_g_.insert(SCF_g_.end(), SCF_storage.begin(), SCF_storage.end());
		}
		if ( G.get_group_type() == "mxy" ) {
			SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"P");
			SCF_P_.insert(SCF_P_.end(), SCF_storage.begin(), SCF_storage.end());

		}

		if ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" ){

			SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"Ekin");
			SCF_Ekin_.insert(SCF_Ekin_.end(), SCF_storage.begin(), SCF_storage.end());

			SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"Eint");
			SCF_Eint_.insert(SCF_Eint_.end(), SCF_storage.begin(), SCF_storage.end());

			SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"W");
			SCF_W_.insert(SCF_W_.end(), SCF_storage.begin(), SCF_storage.end());

			SCF_storage = G.calc_SCF_averaged(par_.rbin(),par_.n_rsamps(),"E");
			SCF_E_.insert(SCF_E_.end(), SCF_storage.begin(), SCF_storage.end());
		}

		if ( store_mxq_){
			store_double = bin_qvals_to_q(vector_pw_norm(mxq_cur_));
			chimxq_.insert(chimxq_.end(), store_double.begin(), store_double.end());
			if ( store_myq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(mxq_cur_,myq_cur_));
				SCFq_xy_.insert(SCFq_xy_.end(), store_complex.begin(), store_complex.end());

				store_double = bin_qvals_to_q(vector_pw_norm(mparq_cur_));
				chimparq_.insert(chimparq_.end(), store_double.begin(), store_double.end());

				store_double = bin_qvals_to_q(vector_pw_norm(mperpq_cur_));
				chimperpq_.insert(chimperpq_.end(), store_double.begin(), store_double.end());

				store_complex = bin_qvals_to_q(vector_pw_mult(mparq_cur_,mperpq_cur_));
				SCFq_mparmperp_.insert(SCFq_mparmperp_.end(), store_complex.begin(), store_complex.end());
			}
			if ( store_wq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" ) ){
				store_complex = bin_qvals_to_q(vector_pw_mult(mxq_cur_,wq_cur_));
				SCFq_xw_.insert(SCFq_xw_.end(), store_complex.begin(), store_complex.end());
			}
			if ( store_eq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" ) ){
				store_complex = bin_qvals_to_q(vector_pw_mult(mxq_cur_,eq_cur_));
				SCFq_xe_.insert(SCFq_xe_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_myq_){
			store_double = bin_qvals_to_q(vector_pw_norm(myq_cur_));
			chimyq_.insert(chimyq_.end(), store_double.begin(), store_double.end());
			if ( store_wq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" )){
				store_complex = bin_qvals_to_q(vector_pw_mult(myq_cur_,wq_cur_));
				SCFq_yw_.insert(SCFq_yw_.end(), store_complex.begin(), store_complex.end());
			}
			if ( store_eq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" )){
				store_complex = bin_qvals_to_q(vector_pw_mult(myq_cur_,eq_cur_));
				SCFq_ye_.insert(SCFq_ye_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_wq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" )){
			store_double = bin_qvals_to_q(vector_pw_norm(wq_cur_));
			chiwq_.insert(chiwq_.end(), store_double.begin(), store_double.end());
			if ( store_eq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(wq_cur_,eq_cur_));
				SCFq_we_.insert(SCFq_we_.end(), store_complex.begin(), store_complex.end());
			}			}
		if ( store_eq_ && ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy")){
			store_double = bin_qvals_to_q(vector_pw_norm(eq_cur_));
			chieq_.insert(chieq_.end(), store_double.begin(), store_double.end());
			if ( store_rq_ ) {
				store_complex = bin_qvals_to_q(vector_pw_mult(rq_cur_,eq_cur_));
				SCFq_re_.insert(SCFq_re_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_rq_){
			store_double = bin_qvals_to_q(vector_pw_norm(rq_cur_));
			chirq_.insert(chirq_.end(), store_double.begin(), store_double.end());
		}
		if ( store_jparq_ && G.get_group_type() == "mxy" ){
			store_double = bin_qvals_to_q(vector_pw_norm(jparq_cur_));
			chijparq_.insert(chijparq_.end(), store_double.begin(), store_double.end());
		}
		if ( store_jperpq_ && G.get_group_type() == "mxy" ){
			store_double = bin_qvals_to_q(vector_pw_norm(jperpq_cur_));
			chijperpq_.insert(chijperpq_.end(), store_double.begin(), store_double.end());
		}
		if ( store_lq_ && G.get_group_type() == "mxy" ){
			store_double = bin_qvals_to_q(vector_pw_norm(lq_cur_));
			chilq_.insert(chilq_.end(), store_double.begin(), store_double.end());
		}

		if ( store_teq_){
			store_double = bin_qvals_to_q(vector_pw_norm(teq_cur_));
			chiteq_.insert(chiteq_.end(), store_double.begin(), store_double.end());
		}
	}
	if ( store_TransCoeff_ && G.get_group_type() == "mxy" ) { // Transport coefficients only defined for mxy group
		TransCoeff_J1_.push_back(			G.calc_neighbor_mean(0, 0, 0, 0, 1, 0, 0) / G.get_N());
		TransCoeff_J1_cos1_.push_back(		G.calc_neighbor_mean(0, 0, 1, 0, 1, 0, 0) / G.get_N());
		TransCoeff_J1_cos2_.push_back( 		G.calc_neighbor_mean(0, 0, 2, 0, 1, 0, 0) / G.get_N());
		TransCoeff_J1_cos2_r2_.push_back( 	G.calc_neighbor_mean(0, 2, 2, 0, 1, 0, 0) / G.get_N());
		TransCoeff_J1_cos1_r2_.push_back( 	G.calc_neighbor_mean(0, 2, 1, 0, 1, 0, 0) / G.get_N());
		TransCoeff_J1_sin1_te_.push_back( 	G.calc_neighbor_mean(1, 0, 0, 1, 1, 0, 0) / G.get_N());
		TransCoeff_Up_.push_back( 			G.calc_neighbor_mean(0, 0, 0, 0, 0, 1, 0) / G.get_N());
		TransCoeff_Up_rinv_.push_back( 		G.calc_neighbor_mean(0, -1, 0, 0, 0, 1, 0) / G.get_N());
		TransCoeff_Upp_.push_back( 			G.calc_neighbor_mean(0, 0, 0, 0, 0, 0, 1) / G.get_N());
		TransCoeff_Up_cos1_.push_back( 		G.calc_neighbor_mean(0, 0, 1, 0, 0, 1, 0) / G.get_N());
		TransCoeff_Up_rinv_cos1_.push_back(	G.calc_neighbor_mean(0, -1, 1, 0, 0, 1, 0) / G.get_N());
		TransCoeff_Upp_cos1_.push_back( 	G.calc_neighbor_mean(0, 0, 1, 0, 0, 1, 0) / G.get_N());
		TransCoeff_cos1_.push_back( 		G.calc_neighbor_mean(0, 0, 1, 0, 0, 0, 0) / G.get_N());
		TransCoeff_Up_rinv_te2_.push_back( 	G.calc_neighbor_mean(2, -1, 0, 0, 0, 1, 0) / G.get_N());
		TransCoeff_Upp_te2_.push_back( 		G.calc_neighbor_mean(2, 0, 0, 0, 0, 0, 1) / G.get_N());
	}

}



void sampler::sample_TCF(const group& G, const group& G_initial, const double t){
	if ( store_TCF_ ){
		TCF_times_.push_back(t);
		topology::Vector2d curM;
//		refresh_qvals(par_.qfullmax(), par_.qsamps_per_bin());
		std::vector<double> store_double;
		std::vector<std::complex<double> > store_complex;

		ACF_Spin_.push_back(G.calc_ACF_sp(G_initial,"S"));
		ACF_anglediff_.push_back(G.calc_ACF_sp(G_initial,"anglediff"));
		ACF_Sx_.push_back(G.calc_ACF_sp(G_initial,"Sx"));
		ACF_Sy_.push_back(G.calc_ACF_sp(G_initial,"Sy"));
		ACF_q0_M_.push_back(G.calc_ACF_q0(G_initial,"M"));
		ACF_q0_absM_.push_back(G.calc_ACF_q0(G_initial,"absM"));
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "xy" || G.get_group_type() == "fmxy" ){
			ACF_W_.push_back(G.calc_ACF_sp(G_initial,"W"));
			ACF_E_.push_back(G.calc_ACF_sp(G_initial,"E"));
			ACF_Eint_.push_back(G.calc_ACF_sp(G_initial,"Eint"));
			ACF_Ekin_.push_back(G.calc_ACF_sp(G_initial,"Ekin"));
		}
		if ( G.get_group_type() == "mxy" ){
			ACF_P_.push_back(G.calc_ACF_sp(G_initial,"P"));
			ACF_Ppar_.push_back(G.calc_ACF_sp(G_initial,"Ppar"));
			ACF_Pperp_.push_back(G.calc_ACF_sp(G_initial,"Pperp"));
		}
//		ACF_Sx_.push_back(G.calc_ACF_sp(G_initial,"Sx"));

//		std::vector<std::complex<double> > mxq_initial, myq_initial, wq_initial, eq_initial;
//		std::chrono::duration<double> elapsed;
//		auto start = std::chrono::high_resolution_clock::now();
		if ( store_mxq_ && ( mxq_initial_.empty()  || refresh_q_ ) ){
			mxq_initial_ = G_initial.calc_fieldfluct(qvals_,"mxq");
		}
		if ( store_myq_ && ( myq_initial_.empty()  || refresh_q_ ) ){
			myq_initial_ = G_initial.calc_fieldfluct(qvals_,"myq");
		}
		if ( ( store_mxq_ && store_myq_ ) && (mparq_initial_.empty() || refresh_q_) ) {
			if (store_static_){
				curM = M_.back();
			} else {
				curM = G.sum_s()/G.get_N();
			}
			double xcomp = curM.get_x()/std::sqrt(curM.norm2());
			double ycomp = curM.get_y()/std::sqrt(curM.norm2());
			mparq_initial_ = vector_sum(vector_scale(mxq_initial_,xcomp),vector_scale(myq_initial_,ycomp));
			mperpq_initial_ = vector_sum(vector_scale(mxq_initial_,-ycomp),vector_scale(myq_initial_,xcomp));
		}
		if ( store_wq_ && ( wq_initial_.empty()  || refresh_q_ ) ){
			wq_initial_ = G_initial.calc_fieldfluct(qvals_,"wq");
		}
		if ( store_eq_ && ( eq_initial_.empty()  || refresh_q_ ) ){
			eq_initial_ = G_initial.calc_fieldfluct(qvals_,"eq");
		}
		if ( store_rq_ && ( rq_initial_.empty()  || refresh_q_ ) ){
			rq_initial_ = G_initial.calc_fieldfluct(qvals_,"rq");
		}
		if ( store_jparq_ && ( jparq_initial_.empty()  || refresh_q_ ) ){
			jparq_initial_ = G_initial.calc_fieldfluct(qvals_,"jparq");
		}
		if ( store_jperpq_ && ( jperpq_initial_.empty()  || refresh_q_ ) ){
			jperpq_initial_ = G_initial.calc_fieldfluct(qvals_,"jperpq");
		}
		if ( store_lq_ && ( lq_initial_.empty()  || refresh_q_ ) ){
			lq_initial_ = G_initial.calc_fieldfluct(qvals_,"lq");
		}
//		std::cout << "TCF:  " << std::endl;
//		elapsed = std::chrono::high_resolution_clock::now() - start;
//		std::cout << "Computation time TFF:  " << elapsed.count() << std::endl;

//		std::vector<std::complex<double> > mxq_new = G_new.calc_fieldfluct(qvals_,"mxq");
//		std::vector<std::complex<double> > myq_new = G_new.calc_fieldfluct(qvals_,"myq");
//		std::vector<std::complex<double> > wq_new = G_new.calc_fieldfluct(qvals_,"wq");
//		std::vector<std::complex<double> > eq_new = G_new.calc_fieldfluct(qvals_,"eq");
		if ( store_mxq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(mxq_initial_,mxq_cur_));
			gxx_.insert(gxx_.end(), store_complex.begin(), store_complex.end());
			if ( store_myq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(mxq_initial_,myq_cur_));
				gxy_.insert(gxy_.end(), store_complex.begin(), store_complex.end());

				store_complex = bin_qvals_to_q(vector_pw_mult(mparq_initial_,mparq_cur_));
				gmparmpar_.insert(gmparmpar_.end(), store_complex.begin(), store_complex.end());

				store_complex = bin_qvals_to_q(vector_pw_mult(mparq_initial_,mperpq_cur_));
				gmparmperp_.insert(gmparmperp_.end(), store_complex.begin(), store_complex.end());

				store_complex = bin_qvals_to_q(vector_pw_mult(mperpq_initial_,mperpq_cur_));
				gmperpmperp_.insert(gmperpmperp_.end(), store_complex.begin(), store_complex.end());
			}
			if ( store_wq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(mxq_initial_,wq_cur_));
				gxw_.insert(gxw_.end(), store_complex.begin(), store_complex.end());
			}
			if ( store_eq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(mxq_initial_,eq_cur_));
				gxe_.insert(gxe_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_myq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(myq_initial_,myq_cur_));
			gyy_.insert(gyy_.end(), store_complex.begin(), store_complex.end());
			if ( store_wq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(myq_initial_,wq_cur_));
				gyw_.insert(gyw_.end(), store_complex.begin(), store_complex.end());
			}
			if ( store_eq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(myq_initial_,eq_cur_));
				gye_.insert(gye_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_wq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(wq_initial_,wq_cur_));
			gww_.insert(gww_.end(), store_complex.begin(), store_complex.end());
			if ( store_eq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(wq_initial_,eq_cur_));
				gwe_.insert(gwe_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_eq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(eq_initial_,eq_cur_));
			gee_.insert(gee_.end(), store_complex.begin(), store_complex.end());
			if ( store_rq_){
				store_complex = bin_qvals_to_q(vector_pw_mult(rq_initial_,eq_cur_));
				gre_.insert(gre_.end(), store_complex.begin(), store_complex.end());
			}
		}
		if ( store_teq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(teq_initial_,teq_cur_));
			gtt_.insert(gtt_.end(), store_complex.begin(), store_complex.end());
		}
		if ( store_rq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(rq_initial_,rq_cur_));
			grr_.insert(grr_.end(), store_complex.begin(), store_complex.end());
		}
		if ( store_jparq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(jparq_initial_,jparq_cur_));
			gjparjpar_.insert(gjparjpar_.end(), store_complex.begin(), store_complex.end());
		}
		if ( store_jperpq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(jperpq_initial_,jperpq_cur_));
			gjperpjperp_.insert(gjperpjperp_.end(), store_complex.begin(), store_complex.end());
		}
		if ( store_lq_){
			store_complex = bin_qvals_to_q(vector_pw_mult(lq_initial_,lq_cur_));
			gll_.insert(gll_.end(), store_complex.begin(), store_complex.end());
		}

	}



		if ( store_MemoryKernels_ && G.get_group_type() == "mxy" ){
			std::vector<std::complex<double> > store_complex;
			if (  mxq_initial_.empty()  || ( refresh_q_ && ! store_TCF_ ) ){
				mxq_initial_ = G_initial.calc_fieldfluct(qvals_,"mxq");
			}
			if (  myq_initial_.empty()  || ( refresh_q_ && ! store_TCF_ ) ){
				myq_initial_ = G_initial.calc_fieldfluct(qvals_,"myq");
			}
			if (  wq_initial_.empty()  || ( refresh_q_ && ! store_TCF_ ) ){
				wq_initial_ = G_initial.calc_fieldfluct(qvals_,"wq");
			}
			if (  jparq_initial_.empty()  || ( refresh_q_ && ! store_TCF_ ) ){
				jparq_initial_ = G_initial.calc_fieldfluct(qvals_,"jparq");
			}
			if ( ! store_mxq_ || ! store_TCF_ ){
				mxq_cur_ = G.calc_fieldfluct(qvals_,"mxq");
			}
			if ( ! store_myq_ || ! store_TCF_ ){
				myq_cur_ = G.calc_fieldfluct(qvals_,"myq");
			}
			if ( ! store_wq_ || ! store_TCF_ ){
				wq_cur_ = G.calc_fieldfluct(qvals_,"wq");
			}
			if ( ! store_jparq_ || ! store_TCF_ ){
				jparq_cur_ = G.calc_fieldfluct(qvals_,"jparq");
			}
			if (  convol_wmx_initial_.empty()  || refresh_q_ ) {
				convol_wmx_initial_ = G_initial.calc_fieldfluct_convolution(qvals_,"wq","mxq");
				convol_wmy_initial_ = G_initial.calc_fieldfluct_convolution(qvals_,"wq","myq");
				convol_jparmx_initial_ = G_initial.calc_fieldfluct_convolution(qvals_,"jparq","myq");
				convol_jparmy_initial_ = G_initial.calc_fieldfluct_convolution(qvals_,"jparq","mxq");
			}

			convol_wmx_cur_ = G.calc_fieldfluct_convolution(qvals_,"wq","mxq");
			convol_wmy_cur_ = G.calc_fieldfluct_convolution(qvals_,"wq","myq");
			convol_jparmx_cur_ = G.calc_fieldfluct_convolution(qvals_,"jparq","myq");
			convol_jparmy_cur_ = G.calc_fieldfluct_convolution(qvals_,"jparq","mxq");


			store_complex = bin_qvals_to_q(vector_pw_mult(wq_initial_,mxq_initial_,wq_cur_,mxq_cur_));
			K_wmx_.insert(K_wmx_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(wq_initial_,myq_initial_,wq_cur_,myq_cur_));
			K_wmy_.insert(K_wmy_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(jparq_initial_,mxq_initial_,jparq_cur_,mxq_cur_));
			K_jmx_.insert(K_jmx_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(jparq_initial_,myq_initial_,jparq_cur_,myq_cur_));
			K_jmy_.insert(K_jmy_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(wq_initial_,mxq_initial_,jparq_cur_,myq_cur_));
			K_wmx_jmy_.insert(K_wmx_jmy_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(wq_initial_,myq_initial_,jparq_cur_,mxq_cur_));
			K_wmy_jmx_.insert(K_wmy_jmx_.end(), store_complex.begin(), store_complex.end());

			store_complex = bin_qvals_to_q(vector_pw_mult(convol_wmx_cur_,convol_wmx_cur_));
			K_convol_wmx_.insert(K_convol_wmx_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(convol_wmy_cur_,convol_wmy_cur_));
			K_convol_wmy_.insert(K_convol_wmy_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(convol_jparmx_cur_,convol_jparmx_cur_));
			K_convol_jmx_.insert(K_convol_jmx_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(convol_jparmy_cur_,convol_jparmy_cur_));
			K_convol_jmy_.insert(K_convol_jmy_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(convol_wmx_cur_,convol_jparmy_cur_));
			K_convol_wmx_jmy_.insert(K_convol_wmx_jmy_.end(), store_complex.begin(), store_complex.end());
			store_complex = bin_qvals_to_q(vector_pw_mult(convol_wmy_cur_,convol_jparmx_cur_));
			K_convol_wmy_jmx_.insert(K_convol_wmy_jmx_.end(), store_complex.begin(), store_complex.end());
		}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SAMPLE MSD
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sampler::sample_MSD(const std::vector<double>& MSD){
	if ( store_TCF_ ) {
		double ang_MSD_aux = 0;
		double pos_MSD_aux = 0;
		if ( par_.system() == "xy" || par_.system() == "fmxy" || par_.system() == "fvm" ) {
			for (int i = 0; i < par_.N(); i++){
				ang_MSD_aux += MSD[i] * MSD[i];
			}
			ACF_ang_MSD_.push_back(ang_MSD_aux / par_.N());
		} else if ( ( par_.system() == "mxy" ) || ( par_.system() == "vm" ) ) {
			for (int i = 0; i < par_.N(); i++){
				ang_MSD_aux += MSD[i] * MSD[i];
				pos_MSD_aux += MSD[par_.N() + i] * MSD[par_.N() + i] + MSD[2 * par_.N() + i] * MSD[2 * par_.N() + i];
			}
			ACF_ang_MSD_.push_back(ang_MSD_aux / par_.N());
			ACF_MSD_.push_back(pos_MSD_aux / par_.N());
		}
	}
}


/* template<class GROUP>
void sampler::sample_TransCoeff(const GROUP& G, const double t){
	if ( par_.system() == "xy" && store_TransCoeff_ ) {
		TransCoeff_times_.push_back(t);
		if (t == 0){
			tau_initial_ = G.calc_current("tau");
			tau_.push_back(tau_initial_);
			je_initial_ = G.calc_current("je");
			je_.push_back(je_initial_);
		} else {
			tau_.push_back(G.calc_current("tau"));
			je_.push_back(G.calc_current("je"));
		}
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 2; j++){
				Gamma_TransCoeff_.push_back(tau_initial_[i]*tau_.back()[j]);
				kappa_TransCoeff_.push_back(je_initial_[i]*je_.back()[j]);
			}
		}
	} else if ( par_.system() == "mxy" && store_TransCoeff_ ) {
		TransCoeff_J1_ =			G.calc_neighbor_mean(0, 0, 0, 0, 1, 0, 0);
		TransCoeff_J1_cos1_ =		G.calc_neighbor_mean(0, 0, 1, 0, 1, 0, 0);
		TransCoeff_J1_cos2_ = 		G.calc_neighbor_mean(0, 0, 2, 0, 1, 0, 0);
		TransCoeff_J1_cos2_r2_ = 	G.calc_neighbor_mean(0, 2, 2, 0, 1, 0, 0);
		TransCoeff_J1_cos1_r2_ = 	G.calc_neighbor_mean(0, 2, 1, 0, 1, 0, 0);
		TransCoeff_J1_sin1_te_ = 	G.calc_neighbor_mean(1, 0, 0, 1, 1, 0, 0);
		TransCoeff_Up_ = 			G.calc_neighbor_mean(0, 0, 0, 0, 0, 1, 0);
		TransCoeff_Up_rinv_ = 		G.calc_neighbor_mean(0, -1, 0, 0, 0, 1, 0);
		TransCoeff_Upp_ = 			G.calc_neighbor_mean(0, 0, 0, 0, 0, 0, 1);
		TransCoeff_Up_cos1_ = 		G.calc_neighbor_mean(0, 0, 1, 0, 0, 1, 0);
		TransCoeff_Up_rinv_cos1_ = 	G.calc_neighbor_mean(0, -1, 1, 0, 0, 1, 0);
		TransCoeff_Upp_cos1_ = 		G.calc_neighbor_mean(0, 0, 1, 0, 0, 1, 0);
		TransCoeff_cos1_ = 			G.calc_neighbor_mean(0, 0, 1, 0, 0, 0, 0);
		TransCoeff_Up_rinv_te2_ = 	G.calc_neighbor_mean(2, -1, 0, 0, 0, 1, 0);
		TransCoeff_Upp_te2_ = 		G.calc_neighbor_mean(2, 0, 0, 0, 0, 0, 1);
	}

}
template void sampler::sample_TransCoeff(const xygroup& G, const double t);
template void sampler::sample_TransCoeff(const mxygroup& G, const double t);


void sampler::sample_TransCoeff(const xygroup& G, const double t){
	if ( par_.system() == "xy" && store_TransCoeff_ ) {
		TransCoeff_times_.push_back(t);
		if (t == 0){
			tau_initial_ = G.calc_current("tau");
			tau_.push_back(tau_initial_);
			je_initial_ = G.calc_current("je");
			je_.push_back(je_initial_);
		} else {
			tau_.push_back(G.calc_current("tau"));
			je_.push_back(G.calc_current("je"));
		}
		for (int i = 0; i < 2; i++){
			for (int j = 0; j < 2; j++){
				Gamma_TransCoeff_.push_back(tau_initial_[i]*tau_.back()[j]);
				kappa_TransCoeff_.push_back(je_initial_[i]*je_.back()[j]);
			}
		}
	}
}

void sampler::sample_TransCoeff(const mxygroup& G, const double t){
	if ( par_.system() == "mxy" && store_TransCoeff_ ) {
		TransCoeff_J1_.push_back(			G.calc_neighbor_mean(0, 0, 0, 0, 1, 0, 0));
		TransCoeff_J1_cos1_.push_back(		G.calc_neighbor_mean(0, 0, 1, 0, 1, 0, 0));
		TransCoeff_J1_cos2_.push_back( 		G.calc_neighbor_mean(0, 0, 2, 0, 1, 0, 0));
		TransCoeff_J1_cos2_r2_.push_back( 	G.calc_neighbor_mean(0, 2, 2, 0, 1, 0, 0));
		TransCoeff_J1_cos1_r2_.push_back( 	G.calc_neighbor_mean(0, 2, 1, 0, 1, 0, 0));
		TransCoeff_J1_sin1_te_.push_back( 	G.calc_neighbor_mean(1, 0, 0, 1, 1, 0, 0));
		TransCoeff_Up_.push_back( 			G.calc_neighbor_mean(0, 0, 0, 0, 0, 1, 0));
		TransCoeff_Up_rinv_.push_back( 		G.calc_neighbor_mean(0, -1, 0, 0, 0, 1, 0));
		TransCoeff_Upp_.push_back( 			G.calc_neighbor_mean(0, 0, 0, 0, 0, 0, 1));
		TransCoeff_Up_cos1_.push_back( 		G.calc_neighbor_mean(0, 0, 1, 0, 0, 1, 0));
		TransCoeff_Up_rinv_cos1_.push_back( 	G.calc_neighbor_mean(0, -1, 1, 0, 0, 1, 0));
		TransCoeff_Upp_cos1_.push_back( 		G.calc_neighbor_mean(0, 0, 1, 0, 0, 1, 0));
		TransCoeff_cos1_.push_back( 			G.calc_neighbor_mean(0, 0, 1, 0, 0, 0, 0));
		TransCoeff_Up_rinv_te2_.push_back( 	G.calc_neighbor_mean(2, -1, 0, 0, 0, 1, 0));
		TransCoeff_Upp_te2_.push_back( 		G.calc_neighbor_mean(2, 0, 0, 0, 0, 0, 1));
	}

}
*/


/////////////////////////////
// sampling: averages
/////////////////////////////
void sampler::average(){
	if ( store_static_){
//		H_.push_back(vector_mean(H_));
		if ( par_.system() == "mxy" || par_.system() == "xy" || par_.system() == "fmxy" ){
			H_.push_back(selective_vector_mean(H_,averaging_times_, par_.av_time_spacing()));
			H_2_.push_back(selective_vector_mean(H_2_,averaging_times_, par_.av_time_spacing()));
			Hint_2_.push_back(selective_vector_mean(Hint_2_,averaging_times_, par_.av_time_spacing()));
			Hkin_2_.push_back(selective_vector_mean(Hkin_2_,averaging_times_, par_.av_time_spacing()));
			W_.push_back(selective_vector_mean(W_,averaging_times_, par_.av_time_spacing()));
			W_2_.push_back(selective_vector_mean(W_2_,averaging_times_, par_.av_time_spacing()));
			temperature_.push_back(selective_vector_mean(temperature_,averaging_times_, par_.av_time_spacing()));
			temperature_squared_.push_back(selective_vector_mean(temperature_squared_,averaging_times_, par_.av_time_spacing()));
			temperature_omega_.push_back(selective_vector_mean(temperature_omega_,averaging_times_, par_.av_time_spacing()));
			temperature_omega_squared_.push_back(selective_vector_mean(temperature_omega_squared_,averaging_times_, par_.av_time_spacing()));
			if ( par_.system() == "mxy" ){
				temperature_p_.push_back(selective_vector_mean(temperature_p_,averaging_times_, par_.av_time_spacing()));
				temperature_p_squared_.push_back(selective_vector_mean(temperature_p_squared_,averaging_times_, par_.av_time_spacing()));
			}
			Upsilon_.push_back(selective_vector_mean(Upsilon_,averaging_times_, par_.av_time_spacing()));
			H_x_.push_back(selective_vector_mean(H_x_,averaging_times_, par_.av_time_spacing()));
			H_y_.push_back(selective_vector_mean(H_y_,averaging_times_, par_.av_time_spacing()));
			I_x_.push_back(selective_vector_mean(I_x_,averaging_times_, par_.av_time_spacing()));
			I_y_.push_back(selective_vector_mean(I_y_,averaging_times_, par_.av_time_spacing()));
			I_x_2_.push_back(selective_vector_mean(I_x_2_,averaging_times_, par_.av_time_spacing()));
			I_y_2_.push_back(selective_vector_mean(I_y_2_,averaging_times_, par_.av_time_spacing()));

		}
		M_.push_back(selective_vector_mean(M_,averaging_times_, par_.av_time_spacing()));
		M_2_.push_back(selective_vector_mean(M_2_,averaging_times_, par_.av_time_spacing()));
		M_4_.push_back(selective_vector_mean(M_4_,averaging_times_, par_.av_time_spacing()));
		absM_.push_back(selective_vector_mean(absM_,averaging_times_, par_.av_time_spacing()));
		M_angle_.push_back(selective_vector_mean(M_angle_,averaging_times_, par_.av_time_spacing()));
		Theta_.push_back(selective_vector_mean(Theta_,averaging_times_, par_.av_time_spacing()));
		Theta_2_.push_back(selective_vector_mean(Theta_2_,averaging_times_, par_.av_time_spacing()));
		Theta_4_.push_back(selective_vector_mean(Theta_4_,averaging_times_, par_.av_time_spacing()));
		Theta_rel_to_M_.push_back(selective_vector_mean(Theta_rel_to_M_,averaging_times_, par_.av_time_spacing()));
		Theta_rel_to_M_2_.push_back(selective_vector_mean(Theta_rel_to_M_2_,averaging_times_, par_.av_time_spacing()));
		Theta_rel_to_M_4_.push_back(selective_vector_mean(Theta_rel_to_M_4_,averaging_times_, par_.av_time_spacing()));
		if ( par_.system() == "mxy" ){
			P_.push_back(selective_vector_mean(P_,averaging_times_, par_.av_time_spacing()));
			P_2_.push_back(selective_vector_mean(P_2_,averaging_times_, par_.av_time_spacing()));
			P_4_.push_back(selective_vector_mean(P_4_,averaging_times_, par_.av_time_spacing()));
		}
		if ( par_.system() == "mxy" || par_.system() == "vm" ) {
			coordination_number_.push_back(selective_vector_mean(coordination_number_,averaging_times_, par_.av_time_spacing()));
		}
	}
	if ( store_vortices_ && (par_.system() == "xy" || par_.system() == "fvm") ){
		abs_vortices_.push_back(selective_vector_mean(abs_vortices_,averaging_times_, par_.av_time_spacing()));
		signed_vortices_.push_back(selective_vector_mean(signed_vortices_,averaging_times_, par_.av_time_spacing()));
	}
	if ( par_.system() == "mxy" && store_TransCoeff_ ) {
			TransCoeff_J1_.push_back(selective_vector_mean(TransCoeff_J1_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_J1_cos1_.push_back(selective_vector_mean(TransCoeff_J1_cos1_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_J1_cos2_.push_back(selective_vector_mean(TransCoeff_J1_cos2_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_J1_cos2_r2_.push_back(selective_vector_mean(TransCoeff_J1_cos2_r2_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_J1_cos1_r2_.push_back(selective_vector_mean(TransCoeff_J1_cos1_r2_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_J1_sin1_te_.push_back(selective_vector_mean(TransCoeff_J1_sin1_te_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Up_.push_back(selective_vector_mean(TransCoeff_Up_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Up_rinv_.push_back(selective_vector_mean(TransCoeff_Up_rinv_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Upp_.push_back(selective_vector_mean(TransCoeff_Upp_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Up_cos1_.push_back(selective_vector_mean(TransCoeff_Up_cos1_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Up_rinv_cos1_.push_back(selective_vector_mean(TransCoeff_Up_rinv_cos1_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Upp_cos1_.push_back(selective_vector_mean(TransCoeff_Upp_cos1_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_cos1_.push_back(selective_vector_mean(TransCoeff_cos1_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Up_rinv_te2_.push_back(selective_vector_mean(TransCoeff_Up_rinv_te2_,averaging_times_, par_.av_time_spacing()));
			TransCoeff_Upp_te2_.push_back(selective_vector_mean(TransCoeff_Upp_te2_,averaging_times_, par_.av_time_spacing()));
	}

	for (size_t k =0; k < par_.qbin().size(); k++){
		if ( store_mxq_){
			mxq_.push_back(selective_column_mean(mxq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_myq_){
			myq_.push_back(selective_column_mean(myq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_myq_){
			mparq_.push_back(selective_column_mean(mparq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
			mperpq_.push_back(selective_column_mean(mperpq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_eq_){
			eq_.push_back(selective_column_mean(eq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_wq_){
			wq_.push_back(selective_column_mean(wq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_teq_){
			teq_.push_back(selective_column_mean(teq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_rq_){
			rq_.push_back(selective_column_mean(rq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_jparq_){
			jparq_.push_back(selective_column_mean(jparq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_jperpq_){
			jperpq_.push_back(selective_column_mean(jperpq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}
		if ( store_lq_){
			lq_.push_back(selective_column_mean(lq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
		}

		if ( store_SCF_){
			if ( store_mxq_){
				chimxq_.push_back(selective_column_mean(chimxq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				if ( store_myq_){
					SCFq_xy_.push_back(selective_column_mean(SCFq_xy_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
					SCFq_mparmperp_.push_back(selective_column_mean(SCFq_mparmperp_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
					chimparq_.push_back(selective_column_mean(chimparq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
					chimperpq_.push_back(selective_column_mean(chimperpq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
				if (store_wq_){
					SCFq_xw_.push_back(selective_column_mean(SCFq_xw_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
				if (store_eq_){
					SCFq_xe_.push_back(selective_column_mean(SCFq_xe_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
			}
			if ( store_myq_){
				chimyq_.push_back(selective_column_mean(chimyq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				if (store_wq_){
					SCFq_yw_.push_back(selective_column_mean(SCFq_yw_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
				if (store_eq_){
					SCFq_ye_.push_back(selective_column_mean(SCFq_ye_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
			}
			if ( store_wq_){
				chiwq_.push_back(selective_column_mean(chiwq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				if (store_eq_){
					SCFq_we_.push_back(selective_column_mean(SCFq_we_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
			}
			if ( store_eq_){
				chieq_.push_back(selective_column_mean(chieq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				if (store_rq_){
					SCFq_re_.push_back(selective_column_mean(SCFq_re_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
				}
			}
			if ( store_teq_){
				chiteq_.push_back(selective_column_mean(chiteq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
			}
			if ( store_rq_){
				chirq_.push_back(selective_column_mean(chirq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
			}
			if ( store_jparq_){
				chijparq_.push_back(selective_column_mean(chijparq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
			}
			if ( store_jperpq_){
				chijperpq_.push_back(selective_column_mean(chijperpq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
			}
			if ( store_lq_){
				chilq_.push_back(selective_column_mean(chilq_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.qbin().size(),k));
			}
		}
	}
	for (size_t l =0; l < par_.rbin().size(); l++){
		if ( store_SCF_ ){
			SCF_Spin_.push_back(selective_column_mean(SCF_Spin_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			SCF_Spin_par_.push_back(selective_column_mean(SCF_Spin_par_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			SCF_Spin_perp_.push_back(selective_column_mean(SCF_Spin_perp_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			SCF_anglediff_.push_back(selective_column_mean(SCF_anglediff_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			if ( par_.system() == "mxy" || par_.system() == "xy" || par_.system() == "fmxy" ){
				SCF_W_.push_back(selective_column_mean(SCF_W_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
				SCF_E_.push_back(selective_column_mean(SCF_E_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
				SCF_Ekin_.push_back(selective_column_mean(SCF_Ekin_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
				SCF_Eint_.push_back(selective_column_mean(SCF_Eint_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			}
			if (par_.system() == "mxy"){
				SCF_P_.push_back(selective_column_mean(SCF_P_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
				SCF_g_.push_back(selective_column_mean(SCF_g_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			}
			if (par_.system() == "vm"){
				SCF_g_.push_back(selective_column_mean(SCF_g_,averaging_times_, par_.av_time_spacing(),averaging_times_.size(),par_.rbin().size(),l));
			}

		}
	}
//	if ( par_.system() == "xy" && store_TransCoeff_ ) {
//		tau_.push_back(selective_vector_mean(tau_,TransCoeff_times_,par_.av_time_spacing()));
//		je_.push_back(selective_vector_mean(je_,TransCoeff_times_,par_.av_time_spacing()));
//		for (size_t i = 0; i < averaging_times_.size(); i++){
//			Gamma_TransCoeff_.push_back(column_mean(Gamma_TransCoeff_,averaging_times_.size(),4,i));
//			kappa_TransCoeff_.push_back(column_mean(kappa_TransCoeff_,averaging_times_.size(),4,i));
//		}
//	}
}

/////////////////////////////
// sampling: snapshots
/////////////////////////////

void sampler::sample_snapshots(const group& G, const double t, std::string snapfile_name, std::ofstream& snapoverview){
	nsnap_++;
	if ( print_snapshots_ ){
		std::ofstream snapfile;
		snapfile.open(snapfile_name, std::ios::out );
		snapoverview << t << " " << snapfile_name << std::endl;
		G.print_group(snapfile);
		snapfile.close();
	}
}
////////////////////////////////////////////////////////////////////////////////////
// printing functions
////////////////////////////////////////////////////////////////////////////////////

void sampler::print_matlab(std::ofstream& outfile){
	print_stdvector(averaging_times_, outfile, "averaging_times = [", "];\n");
	print_stdvector(par_.rbin(), outfile, "rbin = [", "];\n");
	print_stdvector(par_.qbin(), outfile, "qbin = [", "];\n");
	if ( store_static_){
		if ( par_.system() == "mxy" || par_.system() == "xy" || par_.system() == "fmxy" ){
			print_stdvector(H_, outfile, "H = [", "];\n");
			print_stdvector(H_2_, outfile, "H_2 = [", "];\n");
			print_stdvector(Hkin_2_, outfile, "Hkin_2 = [", "];\n");
			print_stdvector(Hint_2_, outfile, "Hint_2 = [", "];\n");
			print_stdvector(temperature_, outfile, "temperature = [", "];\n");
			print_stdvector(temperature_squared_, outfile, "temperature_squared = [", "];\n");
			print_stdvector(temperature_omega_, outfile, "temperature_omega = [", "];\n");
			print_stdvector(temperature_omega_squared_, outfile, "temperature_omega_squared = [", "];\n");
			if ( par_.system() == "mxy" ){
				print_stdvector(temperature_p_, outfile, "temperature_p = [", "];\n");
				print_stdvector(temperature_p_squared_, outfile, "temperature_p_squared = [", "];\n");
			}
			print_stdvector(Upsilon_, outfile, "Upsilon = [", "];\n");
			print_stdvector(H_x_, outfile, "H_x = [", "];\n");
			print_stdvector(H_y_, outfile, "H_y = [", "];\n");
			print_stdvector(I_x_, outfile, "I_x = [", "];\n");
			print_stdvector(I_y_, outfile, "I_y = [", "];\n");
			print_stdvector(I_x_2_, outfile, "I_x_2 = [", "];\n");
			print_stdvector(I_y_2_, outfile, "I_y_2 = [", "];\n");


		}
		print_stdvector(W_, outfile, "W = [", "];\n");
		print_stdvector(W_2_, outfile, "W_2 = [", "];\n");
		topology::print_matlab(M_, "M", outfile);
		print_stdvector(M_2_, outfile, "M_2 = [", "];\n");
		print_stdvector(M_4_, outfile, "M_4 = [", "];\n");
		print_stdvector(absM_, outfile, "absM = [", "];\n");
		print_stdvector(M_angle_, outfile, "M_angle = [", "];\n");
		print_stdvector(Theta_, outfile, "Theta = [", "];\n");
		print_stdvector(Theta_2_, outfile, "Theta_2 = [", "];\n");
		print_stdvector(Theta_4_, outfile, "Theta_4 = [", "];\n");
		print_stdvector(Theta_rel_to_M_, outfile, "Theta_rel_to_M = [", "];\n");
		print_stdvector(Theta_rel_to_M_2_, outfile, "Theta_rel_to_M_2 = [", "];\n");
		print_stdvector(Theta_rel_to_M_4_, outfile, "Theta_rel_to_M_4 = [", "];\n");
		if ( par_.system() == "mxy" ){
			topology::print_matlab(P_, "P", outfile);
			print_stdvector(P_2_, outfile, "P_2 = [", "];\n");
			print_stdvector(P_4_, outfile, "P_4 = [", "];\n");
		}
		if ( par_.system() == "mxy" || par_.system() == "vm" ){
			print_stdvector(coordination_number_, outfile, "coordination_number = [", "];\n");
		}
	}
	if ( (par_.system() == "xy" || par_.system() == "fvm") && store_vortices_){
		print_stdvector(abs_vortices_, outfile, "abs_vortices = [", "];\n");
		print_stdvector(signed_vortices_, outfile, "signed_vortices = [", "];\n");
	}

	if ( par_.system() == "mxy" && store_TransCoeff_ ) {
			print_stdvector(TransCoeff_J1_, outfile, "TransCoeff_J1 = [", "];\n");
			print_stdvector(TransCoeff_J1_cos1_, outfile, "TransCoeff_J1_cos1 = [", "];\n");
			print_stdvector(TransCoeff_J1_cos2_, outfile, "TransCoeff_J1_cos2 = [", "];\n");
			print_stdvector(TransCoeff_J1_cos2_r2_, outfile, "TransCoeff_J1_cos2_r2 = [", "];\n");
			print_stdvector(TransCoeff_J1_cos1_r2_, outfile, "TransCoeff_J1_cos1_r2 = [", "];\n");
			print_stdvector(TransCoeff_J1_sin1_te_, outfile, "TransCoeff_J1_sin1_te = [", "];\n");
			print_stdvector(TransCoeff_Up_, outfile, "TransCoeff_Up = [", "];\n");
			print_stdvector(TransCoeff_Up_rinv_, outfile, "TransCoeff_Up_rinv = [", "];\n");
			print_stdvector(TransCoeff_Upp_, outfile, "TransCoeff_Upp = [", "];\n");
			print_stdvector(TransCoeff_Up_cos1_, outfile, "TransCoeff_Up_cos1 = [", "];\n");
			print_stdvector(TransCoeff_Up_rinv_cos1_, outfile, "TransCoeff_Up_rinv_cos1 = [", "];\n");
			print_stdvector(TransCoeff_Upp_cos1_, outfile, "TransCoeff_Upp_cos1 = [", "];\n");
			print_stdvector(TransCoeff_cos1_, outfile, "TransCoeff_cos1 = [", "];\n");
			print_stdvector(TransCoeff_Up_rinv_te2_, outfile, "TransCoeff_Up_rinv_te2 = [", "];\n");
			print_stdvector(TransCoeff_Upp_te2_, outfile, "TransCoeff_Upp_te2 = [", "];\n");
	}
	/* Quite pointless. Time evolution is stored in TCFs, and these are only
	 * relevant for averages.
	if ( store_mxq_){
		print_stdvector(mxq_, outfile, "mxq = [", "];\n", 0, mxq_.size(), 1);
	}
	if ( store_myq_){
		print_stdvector(myq_, outfile, "myq = [", "];\n", 0, myq_.size(), 1);
	}
	if ( store_eq_){
		print_stdvector(eq_, outfile, "eq = [", "];\n", 0, eq_.size(), 1);
	}
	if ( store_wq_){
		print_stdvector(wq_, outfile, "wq = [", "];\n", 0, wq_.size(), 1);
	}
	if ( store_teq_){
		print_stdvector(teq_, outfile, "teq = [", "];\n", 0, teq_.size(), 1);
	}
	*/

	/* Printing SCF is pretty pointless and just costs a lot of storage capacity
	 *
	if ( store_SCF_ ){
		if ( store_mxq_){
			print_stdvector(chimxq_, outfile, "chimxq = [", "];\n", 0, chimxq_.size(), 1);
			if ( store_myq_){
				print_stdvector(SCFq_xy_, outfile, "SCFq_xy = [", "];\n", 0, chimxq_.size(), 1);
			}
			if ( store_wq_){
				print_stdvector(SCFq_xw_, outfile, "SCFq_xw = [", "];\n", 0, chimxq_.size(), 1);
			}
			if ( store_eq_){
				print_stdvector(SCFq_xe_, outfile, "SCFq_xe = [", "];\n", 0, chimxq_.size(), 1);
			}
		}
		if ( store_myq_){
			print_stdvector(chimyq_, outfile, "chimyq = [", "];\n", 0, chimyq_.size(), 1);
			if ( store_wq_){
				print_stdvector(SCFq_yw_, outfile, "SCFq_yw = [", "];\n", 0, chimxq_.size(), 1);
			}
			if ( store_eq_){
				print_stdvector(SCFq_ye_, outfile, "SCFq_ye = [", "];\n", 0, chimxq_.size(), 1);
			}
		}
		if ( store_wq_){
			print_stdvector(chiwq_, outfile, "chiwq = [", "];\n", 0, chiwq_.size(), 1);
			if ( store_eq_){
				print_stdvector(SCFq_we_, outfile, "SCFq_we = [", "];\n", 0, chimxq_.size(), 1);
			}
		}
		if ( store_eq_){
			print_stdvector(chieq_, outfile, "chieq = [", "];\n", 0, chieq_.size(), 1);
		}

		if ( store_teq_){
			print_stdvector(chiteq_, outfile, "chiteq = [", "];\n", 0, chiteq_.size(), 1);
		}
		print_stdvector(SCF_Spin_, outfile, "SCF_Spin = [", "];\n", 0,SCF_Spin_.size(), 1);
		print_stdvector(SCF_anglediff_, outfile, "SCF_anglediff = [", "];\n", 0,SCF_anglediff_.size(), 1);
		if (par_.system() == "mxy"){
			print_stdvector(SCF_g_, outfile, "SCF_g = [", "];\n", 0,SCF_g_.size(), 1);
			print_stdvector(SCF_P_, outfile, "SCF_P = [", "];\n", 0,SCF_P_.size(), 1);
			print_stdvector(SCF_W_, outfile, "SCF_W = [", "];\n", 0,SCF_W_.size(), 1);
			print_stdvector(SCF_Ekin_, outfile, "SCF_Ekin = [", "];\n", 0,SCF_Ekin_.size(), 1);
			print_stdvector(SCF_Eint_, outfile, "SCF_Eint = [", "];\n", 0,SCF_Eint_.size(), 1);
		}
	}
	*/
	if ( store_TCF_ ){
		print_stdvector(TCF_times_, outfile, "TCF_times = [", "];\n");

		print_stdvector(ACF_Spin_, outfile, "ACF_Spin = [", "];\n");
		print_stdvector(ACF_anglediff_, outfile, "ACF_anglediff = [", "];\n");
		print_stdvector(ACF_Sx_, outfile, "ACF_Sx = [", "];\n");
		print_stdvector(ACF_Sy_, outfile, "ACF_Sy = [", "];\n");

		print_stdvector(ACF_q0_M_, outfile, "ACF_q0_M = [", "];\n");
		print_stdvector(ACF_q0_absM_, outfile, "ACF_q0_absM = [", "];\n");
		print_stdvector(ACF_ang_MSD_, outfile, "ACF_ang_MSD = [", "];\n");

		if ( par_.system() == "mxy" || par_.system() == "xy" || par_.system() == "fmxy" ){
			print_stdvector(ACF_W_, outfile, "ACF_W = [", "];\n");
			print_stdvector(ACF_E_, outfile, "ACF_E = [", "];\n");
			print_stdvector(ACF_Eint_, outfile, "ACF_Eint = [", "];\n");
			print_stdvector(ACF_Ekin_, outfile, "ACF_Ekin = [", "];\n");
		}


		if (par_.system() == "mxy"){
			print_stdvector(ACF_P_, outfile, "ACF_P = [", "];\n");
			print_stdvector(ACF_Ppar_, outfile, "ACF_Ppar = [", "];\n");
			print_stdvector(ACF_Pperp_, outfile, "ACF_Pperp = [", "];\n");
		}
		if ( par_.system() == "mxy" || par_.system() == "vm" ){
			print_stdvector(ACF_MSD_, outfile, "ACF_MSD = [", "];\n");
		}

		if ( store_mxq_){
			print_stdvector(gxx_, outfile, "gxx = [", "];\n", 0, gxx_.size(), 1);
			if ( store_myq_){
				print_stdvector(gxy_, outfile, "gxy = [", "];\n", 0, gxx_.size(), 1);

				print_stdvector(gmparmpar_, outfile, "gmparmpar = [", "];\n", 0, gxx_.size(), 1);
				print_stdvector(gmparmperp_, outfile, "gmparmperp = [", "];\n", 0, gxx_.size(), 1);
				print_stdvector(gmperpmperp_, outfile, "gmperpmperp = [", "];\n", 0, gxx_.size(), 1);
			}
			if ( store_wq_){
				print_stdvector(gxw_, outfile, "gxw = [", "];\n", 0, gxx_.size(), 1);
			}
			if ( store_eq_){
				print_stdvector(gxe_, outfile, "gxe = [", "];\n", 0, gxx_.size(), 1);
			}
		}
		if ( store_myq_){
			print_stdvector(gyy_, outfile, "gyy = [", "];\n", 0, gyy_.size(), 1);
			if ( store_wq_){
				print_stdvector(gyw_, outfile, "gyw = [", "];\n", 0, gyy_.size(), 1);
			}
			if ( store_eq_){
				print_stdvector(gye_, outfile, "gye = [", "];\n", 0, gyy_.size(), 1);
			}
		}
		if ( store_wq_){
			print_stdvector(gww_, outfile, "gww = [", "];\n", 0, gww_.size(), 1);
			if ( store_eq_){
				print_stdvector(gwe_, outfile, "gwe = [", "];\n", 0, gwe_.size(), 1);
			}
		}
		if ( store_eq_){
			print_stdvector(gee_, outfile, "gee = [", "];\n", 0, gee_.size(), 1);
			if ( store_rq_){
				print_stdvector(gre_, outfile, "gre = [", "];\n", 0, gre_.size(), 1);
			}
		}
		if ( store_teq_){
			print_stdvector(gtt_, outfile, "gtt = [", "];\n", 0, gtt_.size(), 1);
		}
		if ( store_rq_){
			print_stdvector(grr_, outfile, "grr = [", "];\n", 0, grr_.size(), 1);
		}
		if ( store_jparq_){
			print_stdvector(gjparjpar_, outfile, "gjparjpar = [", "];\n", 0, gjparjpar_.size(), 1);
		}
		if ( store_jperpq_){
			print_stdvector(gjperpjperp_, outfile, "gjperpjperp = [", "];\n", 0, gjperpjperp_.size(), 1);
		}
		if ( store_lq_){
			print_stdvector(gll_, outfile, "gll = [", "];\n", 0, gll_.size(), 1);
		}
	}
	if ( par_.system() == "xy" && store_TransCoeff_ ) {
		print_stdvector(TransCoeff_times_, outfile, "TransCoeff_times = [", "];\n");
		topology::print_matlab(tau_, "tau", outfile);
		topology::print_matlab(je_, "je", outfile);
//		print_stdvector(tau_, outfile, "tau = [", "];\n", 0, tau_.size(), 1);
//		print_stdvector(je_, outfile, "je = [", "];\n", 0, je_.size(), 1);
		print_stdvector(Gamma_TransCoeff_, outfile, "Gamma_TransCoeff = [", "];\n", 0, Gamma_TransCoeff_.size(), 1);
		print_stdvector(kappa_TransCoeff_, outfile, "kappa_TransCoeff = [", "];\n", 0, kappa_TransCoeff_.size(), 1);
		print_stdvector(Gamma_TransCoeff_new_, outfile, "Gamma_TransCoeff_new = [", "];\n", 0, Gamma_TransCoeff_new_.size(), 1);
		print_stdvector(kappa_TransCoeff_new_, outfile, "kappa_TransCoeff_new = [", "];\n", 0, kappa_TransCoeff_new_.size(), 1);
	}
	if ( par_.system() == "mxy" && store_MemoryKernels_ ) {
		print_stdvector(K_wmx_, outfile, "K_wmx = [", "];\n", 0, K_wmx_.size(), 1);
		print_stdvector(K_wmy_, outfile, "K_wmy = [", "];\n", 0, K_wmy_.size(), 1);
		print_stdvector(K_jmx_, outfile, "K_jparmx = [", "];\n", 0, K_jmx_.size(), 1);
		print_stdvector(K_jmy_, outfile, "K_jparmy = [", "];\n", 0, K_jmy_.size(), 1);
		print_stdvector(K_wmx_jmy_, outfile, "K_wmx_jmy = [", "];\n", 0, K_wmx_jmy_.size(), 1);
		print_stdvector(K_wmy_jmx_, outfile, "K_wmy_jmx = [", "];\n", 0, K_wmy_jmx_.size(), 1);

		print_stdvector(K_convol_wmx_, outfile, "K_convol_wmx = [", "];\n", 0, K_convol_wmx_.size(), 1);
		print_stdvector(K_convol_wmy_, outfile, "K_convol_wmy = [", "];\n", 0, K_convol_wmy_.size(), 1);
		print_stdvector(K_convol_jmx_, outfile, "K_convol_jparmx = [", "];\n", 0, K_convol_jmx_.size(), 1);
		print_stdvector(K_convol_jmy_, outfile, "K_convol_jparmy = [", "];\n", 0, K_convol_jmy_.size(), 1);
		print_stdvector(K_convol_wmx_jmy_, outfile, "K_convol_wmx_jmy = [", "];\n", 0, K_convol_wmx_jmy_.size(), 1);
		print_stdvector(K_convol_wmy_jmx_, outfile, "K_convol_wmy_jmx = [", "];\n", 0, K_convol_wmy_jmx_.size(), 1);

	}
}
/*
 *
 *
 *
 * for (size_t k =0; k < par_.qbin().size(); k++){
		if ( store_mxq_){
			print_stdvector(mxq_, outfile, "mxq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,mxq_.size(), par_.qbin().size());
		}
		if ( store_myq_){
			print_stdvector(myq_, outfile, "myq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,myq_.size(), par_.qbin().size());
		}
		if ( store_eq_){
			print_stdvector(eq_, outfile, "eq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,eq_.size(), par_.qbin().size());
		}
		if ( store_wq_){
			print_stdvector(wq_, outfile, "wq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,wq_.size(), par_.qbin().size());
		}
		if ( store_teq_){
			print_stdvector(teq_, outfile, "teq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,teq_.size(), par_.qbin().size());
		}

		if ( store_chimxq_){
			print_stdvector(chimxq_, outfile, "chimxq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,chimxq_.size(), par_.qbin().size());
		}
		if ( store_chimyq_){
			print_stdvector(chimyq_, outfile, "chimyq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,chimyq_.size(), par_.qbin().size());
		}
		if ( store_chieq_){
			print_stdvector(chieq_, outfile, "chieq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,chieq_.size(), par_.qbin().size());
		}
		if ( store_chiwq_){
			print_stdvector(chiwq_, outfile, "chiwq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,chiwq_.size(), par_.qbin().size());
		}
		if ( store_chiteq_){
			print_stdvector(chiteq_, outfile, "chiteq(" + std::to_string((int)k+1) + ",:) = [", "];\n",k,chiteq_.size(), par_.qbin().size());
		}
	}
	for (size_t l =0; l < par_.rbin().size(); l++){
		if ( store_SCF_ ){
			print_stdvector(SCF_Spin_, outfile, "SCF_Spin(" + std::to_string((int)l+1) + ",:) = [", "];\n",l,SCF_Spin_.size(), par_.rbin().size());
			print_stdvector(SCF_anglediff_, outfile, "SCF_anglediff(" + std::to_string((int)l+1) + ",:) = [", "];\n",l,SCF_anglediff_.size(), par_.rbin().size());
		}
	}

	for (size_t l =0; l < par_.rbin().size(); l++){
		if ( store_SCF_ ){
			print_stdvector(SCF_Spin_, outfile, "SCF_Spin(" + std::to_string((int)l+1) + ",:) = [", "];\n",l,SCF_Spin_.size(), par_.rbin().size());
			print_stdvector(SCF_anglediff_, outfile, "SCF_anglediff(" + std::to_string((int)l+1) + ",:) = [", "];\n",l,SCF_anglediff_.size(), par_.rbin().size());
		}
	}
 */



void sampler::print_averages_matlab(std::ofstream& outfile){
	average();
	if ( store_static_){
		if ( par_.system() == "mxy" || par_.system() == "xy" || par_.system() == "fmxy" ){
			outfile << "H_av = " << H_.back() << ";" << std::endl;
			outfile << "H_2_av = " << H_2_.back() << ";" << std::endl;
			outfile << "Hkin_2_av = " << Hkin_2_.back() << ";" << std::endl;
			outfile << "Hint_2_av = " << Hint_2_.back() << ";" << std::endl;
			outfile << "H_var = " << H_2_.back() - H_.back() * H_.back() << ";" << std::endl;
			outfile << "C_0 = " << (H_2_.back() - H_.back() * H_.back())/( temperature_.back() * temperature_.back()) << ";" << std::endl;
			outfile << "W_av = " << W_.back() << ";" << std::endl;
			outfile << "W_2_av = " << W_2_.back() << ";" << std::endl;
			outfile << "W_var = " << W_2_.back() - W_.back() * W_.back() << ";" << std::endl;
			outfile << "temperature_av = " << temperature_.back() << ";" << std::endl;
			outfile << "temperature_squared_av = " << temperature_squared_.back() << ";" << std::endl;
			outfile << "temperature_var = " << temperature_squared_.back() - std::pow(temperature_.back(),2) << ";" << std::endl;
			outfile << "temperature_omega_av = " << temperature_omega_.back() << ";" << std::endl;
			outfile << "temperature_omega_squared_av = " << temperature_omega_squared_.back() << ";" << std::endl;
			outfile << "temperature_omega_var = " << temperature_omega_squared_.back() - std::pow(temperature_omega_.back(),2) << ";" << std::endl;
			if ( par_.system() == "mxy" ){
				outfile << "temperature_p_av = " << temperature_p_.back() << ";" << std::endl;
				outfile << "temperature_p_squared_av = " << temperature_p_squared_.back() << ";" << std::endl;
				outfile << "temperature_p_var = " << temperature_p_squared_.back() - std::pow(temperature_p_.back(),2) << ";" << std::endl;
			}
			outfile << "Upsilon_av = " << Upsilon_.back() << ";" << std::endl;
			outfile << "H_x_av = " << H_x_.back() << ";" << std::endl;
			outfile << "H_y_av = " << H_y_.back() << ";" << std::endl;
			outfile << "I_x_av = " << I_x_.back() << ";" << std::endl;
			outfile << "I_y_av = " << I_y_.back() << ";" << std::endl;
			outfile << "I_x_2_av = " << I_x_2_.back() << ";" << std::endl;
			outfile << "I_y_2_av = " << I_y_2_.back() << ";" << std::endl;


		}
		outfile << "M_av = [" << (M_.back()).get_x() << ", " << (M_.back()).get_y() << "];" << std::endl;
		outfile << "M_2_av = " << M_2_.back() << ";" << std::endl;
		outfile << "M_4_av = " << M_4_.back() << ";" << std::endl;
		outfile << "M_var = " << M_2_.back() - M_.back().norm2() << ";" << std::endl;
		outfile << "Binder_cum = " << 1 - M_4_.back() / (3 * M_2_.back() * M_2_.back()) << ";" << std::endl;
		outfile << "absM_av = " << absM_.back() << ";" << std::endl;
		outfile << "absM_var = " << M_2_.back() - absM_.back() * absM_.back() << ";" << std::endl;
		outfile << "M_angle_av = " << M_angle_.back() << ";" << std::endl;
		outfile << "Theta_av = " << Theta_.back() << ";" << std::endl;
		outfile << "Theta_2_av = " << Theta_2_.back() << ";" << std::endl;
		outfile << "Theta_4_av = " << Theta_4_.back() << ";" << std::endl;
		outfile << "Theta_rel_to_M_av = " << Theta_rel_to_M_.back() << ";" << std::endl;
		outfile << "Theta_rel_to_M_2_av = " << Theta_rel_to_M_2_.back() << ";" << std::endl;
		outfile << "Theta_rel_to_M_4_av = " << Theta_rel_to_M_4_.back() << ";" << std::endl;
		if ( par_.system() == "mxy" ){
			outfile << "P_av = [" << (P_.back()).get_x() << ", " << (P_.back()).get_y() << "];" << std::endl;
			outfile << "P_2_av = " << P_2_.back() << ";" << std::endl;
			outfile << "P_4_av = " << P_4_.back() << ";" << std::endl;
		}
		if ( par_.system() == "mxy" || par_.system() == "vm" ){
			outfile << "coordination_number_av = " << coordination_number_.back() << ";" << std::endl;
		}
	}
	if ( (par_.system() == "xy" || par_.system() == "fvm") && store_vortices_){
		outfile << "abs_vortices_av = " << abs_vortices_.back() << ";" << std::endl;
		outfile << "signed_vortices_av = " << signed_vortices_.back() << ";" << std::endl;
	}
	if ( par_.system() == "mxy" && store_TransCoeff_ ) {
		outfile << "TransCoeff_J1_av = " << TransCoeff_J1_.back() << ";" << std::endl;
		outfile << "TransCoeff_J1_cos1_av = " << TransCoeff_J1_cos1_.back() << ";" << std::endl;
		outfile << "TransCoeff_J1_cos2_av = " << TransCoeff_J1_cos2_.back() << ";" << std::endl;
		outfile << "TransCoeff_J1_cos2_r2_av = " << TransCoeff_J1_cos2_r2_.back() << ";" << std::endl;
		outfile << "TransCoeff_J1_cos1_r2_av = " << TransCoeff_J1_cos1_r2_.back() << ";" << std::endl;
		outfile << "TransCoeff_J1_sin1_te_av = " << TransCoeff_J1_sin1_te_.back() << ";" << std::endl;
		outfile << "TransCoeff_Up_av = " << TransCoeff_Up_.back() << ";" << std::endl;
		outfile << "TransCoeff_Up_rinv_av = " << TransCoeff_Up_rinv_.back() << ";" << std::endl;
		outfile << "TransCoeff_Upp_av = " << TransCoeff_Upp_.back() << ";" << std::endl;
		outfile << "TransCoeff_Up_cos1_av = " << TransCoeff_Up_cos1_.back() << ";" << std::endl;
		outfile << "TransCoeff_Up_rinv_cos1_av = " << TransCoeff_Up_rinv_cos1_.back() << ";" << std::endl;
		outfile << "TransCoeff_Upp_cos1_av = " << TransCoeff_Upp_cos1_.back() << ";" << std::endl;
		outfile << "TransCoeff_cos1_av = " << TransCoeff_cos1_.back() << ";" << std::endl;
		outfile << "TransCoeff_Up_rinv_te2_av = " << TransCoeff_Up_rinv_te2_.back() << ";" << std::endl;
		outfile << "TransCoeff_Upp_te2_av = " << TransCoeff_Upp_te2_.back() << ";" << std::endl;
	}

	if ( store_mxq_){
		print_stdvector(mxq_, outfile, "mxq_av = [", "];\n",mxq_.size()-par_.qbin().size(),mxq_.size(), 1);
	}
	if ( store_myq_){
		print_stdvector(myq_, outfile, "myq_av = [", "];\n",myq_.size()-par_.qbin().size(),myq_.size(), 1);
	}
	if ( store_mxq_ && store_myq_){
		print_stdvector(mparq_, outfile, "mparq_av = [", "];\n",mxq_.size()-par_.qbin().size(),mxq_.size(), 1);
		print_stdvector(mperpq_, outfile, "mperpq_av = [", "];\n",mxq_.size()-par_.qbin().size(),mxq_.size(), 1);
	}
	if ( store_eq_){
		print_stdvector(eq_, outfile, "eq_av = [", "];\n",eq_.size()-par_.qbin().size(),eq_.size(), 1);
	}
	if ( store_wq_){
		print_stdvector(wq_, outfile, "wq_av = [", "];\n",wq_.size()-par_.qbin().size(),wq_.size(), 1);
	}
	if ( store_teq_){
		print_stdvector(teq_, outfile, "teq_av = [", "];\n",teq_.size()-par_.qbin().size(),teq_.size(), 1);
	}
	if ( store_rq_){
		print_stdvector(rq_, outfile, "rq_av = [", "];\n",rq_.size()-par_.qbin().size(),rq_.size(), 1);
	}
	if ( store_jparq_){
		print_stdvector(jparq_, outfile, "jparq_av = [", "];\n",jparq_.size()-par_.qbin().size(),jparq_.size(), 1);
	}
	if ( store_jperpq_){
		print_stdvector(jperpq_, outfile, "jperpq_av = [", "];\n",jperpq_.size()-par_.qbin().size(),jperpq_.size(), 1);
	}
	if ( store_lq_){
		print_stdvector(lq_, outfile, "lq_av = [", "];\n",lq_.size()-par_.qbin().size(),lq_.size(), 1);
	}

	if ( store_SCF_ ){
		if ( store_mxq_){
			print_stdvector(chimxq_, outfile, "chimxq_av = [", "];\n",chimxq_.size()-par_.qbin().size(),chimxq_.size(), 1);
			if ( store_myq_){
				print_stdvector(SCFq_xy_, outfile, "SCFq_xy_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);

				print_stdvector(chimparq_, outfile, "chimparq_av = [", "];\n",chimxq_.size()-par_.qbin().size(),chimxq_.size(), 1);
				print_stdvector(chimperpq_, outfile, "chimperpq_av = [", "];\n",chimxq_.size()-par_.qbin().size(),chimxq_.size(), 1);
				print_stdvector(SCFq_mparmperp_, outfile, "SCFq_mparmperp_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			}
			if ( store_wq_){
				print_stdvector(SCFq_xw_, outfile, "SCFq_xw_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			}
			if ( store_eq_){
				print_stdvector(SCFq_xe_, outfile, "SCFq_xe_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			}
		}

		if ( store_myq_){
			print_stdvector(chimyq_, outfile, "chimyq_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			if ( store_wq_){
				print_stdvector(SCFq_yw_, outfile, "SCFq_yw_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			}
			if ( store_eq_){
				print_stdvector(SCFq_ye_, outfile, "SCFq_ye_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			}
		}
		if ( store_wq_){
			print_stdvector(chiwq_, outfile, "chiwq_av = [", "];\n",chiwq_.size()-par_.qbin().size(),chiwq_.size(), 1);
			if ( store_eq_){
				print_stdvector(SCFq_we_, outfile, "SCFq_we_av = [", "];\n",chimyq_.size()-par_.qbin().size(),chimyq_.size(), 1);
			}
		}
		if ( store_eq_){
			print_stdvector(chieq_, outfile, "chieq_av = [", "];\n",chieq_.size()-par_.qbin().size(),chieq_.size(), 1);
			if ( store_rq_){
				print_stdvector(SCFq_re_, outfile, "SCFq_re_av = [", "];\n",SCFq_re_.size()-par_.qbin().size(),SCFq_re_.size(), 1);
			}
		}
		if ( store_rq_){
			print_stdvector(chirq_, outfile, "chirq_av = [", "];\n",chirq_.size()-par_.qbin().size(),chirq_.size(), 1);
		}
		if ( store_jparq_){
			print_stdvector(chijparq_, outfile, "chijparq_av = [", "];\n",chijparq_.size()-par_.qbin().size(),chijparq_.size(), 1);
		}
		if ( store_jperpq_){
			print_stdvector(chijperpq_, outfile, "chijperpq_av = [", "];\n",chijperpq_.size()-par_.qbin().size(),chijperpq_.size(), 1);
		}
		if ( store_lq_){
			print_stdvector(chilq_, outfile, "chilq_av = [", "];\n",chilq_.size()-par_.qbin().size(),chilq_.size(), 1);
		}
		if ( store_teq_){
			print_stdvector(chiteq_, outfile, "chiteq_av = [", "];\n",chiteq_.size()-par_.qbin().size(),chiteq_.size(), 1);
		}

		print_stdvector(SCF_Spin_, outfile, "SCF_Spin_av = [", "];\n",SCF_Spin_.size()-par_.rbin().size(),SCF_Spin_.size(), 1);
		print_stdvector(SCF_Spin_par_, outfile, "SCF_Spin_par_av = [", "];\n",SCF_Spin_par_.size()-par_.rbin().size(),SCF_Spin_par_.size(), 1);
		print_stdvector(SCF_Spin_perp_, outfile, "SCF_Spin_perp_av = [", "];\n",SCF_Spin_perp_.size()-par_.rbin().size(),SCF_Spin_perp_.size(), 1);
		print_stdvector(SCF_anglediff_, outfile, "SCF_anglediff_av = [", "];\n",SCF_anglediff_.size()-par_.rbin().size(),SCF_anglediff_.size(), 1);
		print_stdvector(SCF_g_, outfile, "SCF_g_av = [", "];\n",SCF_g_.size()-par_.rbin().size(),SCF_g_.size(), 1);
		if ( par_.system() == "mxy" || par_.system() == "xy" || par_.system() == "fmxy" ){
			print_stdvector(SCF_W_, outfile, "SCF_W_av = [", "];\n",SCF_W_.size()-par_.rbin().size(),SCF_W_.size(), 1);
			print_stdvector(SCF_E_, outfile, "SCF_E_av = [", "];\n",SCF_E_.size()-par_.rbin().size(),SCF_E_.size(), 1);
			print_stdvector(SCF_Ekin_, outfile, "SCF_Ekin_av = [", "];\n",SCF_Ekin_.size()-par_.rbin().size(),SCF_Ekin_.size(), 1);
			print_stdvector(SCF_Eint_, outfile, "SCF_Eint_av = [", "];\n",SCF_Eint_.size()-par_.rbin().size(),SCF_Eint_.size(), 1);
		}
		if ( par_.system() == "mxy" || par_.system() == "vm" ){
			print_stdvector(SCF_g_, outfile, "SCF_g_av = [", "];\n",SCF_g_.size()-par_.rbin().size(),SCF_g_.size(), 1);
		}
		if ( par_.system() == "mxy" ){
			print_stdvector(SCF_P_, outfile, "SCF_P_av = [", "];\n",SCF_P_.size()-par_.rbin().size(),SCF_P_.size(), 1);
		}
	}
	if ( par_.system() == "xy" && store_TransCoeff_ ) {
		std::vector<double> av_vec(2);
		av_vec[0] = tau_.back().get_x();
		av_vec[1] = tau_.back().get_y();
		print_stdvector(av_vec,outfile,"tau_av = [", "];\n");
		av_vec[0] = je_.back().get_x();
		av_vec[1] = je_.back().get_y();
		print_stdvector(av_vec,outfile,"je_av = [", "];\n");
//		outfile << "tau_av = [" << tau_.back().get_x << ", " << tau_.back().get_y() << "];" << std::endl;
//		outfile << "tau_av = [" << tau_.back().get_x << ", " << tau_.back().get_y() << "];" << std::endl;
//		outfile << "je_av = " << je_.back() << ";" << std::endl;
//		print_stdvector(Gamma_TransCoeff_, outfile, "Gamma_TransCoeff_av = [", "];\n", Gamma_TransCoeff_.size() - averaging_times_.size(),Gamma_TransCoeff_.size(), 1);
//		print_stdvector(kappa_TransCoeff_, outfile, "kappa_TransCoeff_av = [", "];\n", kappa_TransCoeff_.size() - averaging_times_.size(),Gamma_TransCoeff_.size(), 1);
	}
}
