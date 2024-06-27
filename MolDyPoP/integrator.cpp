/*! \file integrator.cpp
 *
 *  \brief cpp-file to class declaration of integrator. Implements the routines
 *  declared in integrator.h
 *
 *  Mostly comment-free. Some referenes to the equations given in the papers referenced
 *  in integrator.h.
 *
 *  \author Thomas Bissinger, with contributions by Mathias Höfler

 *  \date Created: 2019-04-12
 *  \date Last Updated: 2023-08-06
 *
 *
 */
#include <vector>
#include "integrator.h"


///////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTOR
///////////////////////////////////////////////////////////////////////////////////////////
integrator::integrator(double dtin, std::string type){
   dt_ = dtin ;
   integrator_type_ = type;
   if ((type != "rk4") && (type != "lf") && (type != "vm") && (type != "langevin")){
	   std::cout << "ERROR during initialization of integrator. Integrator type '" << type  << "' unknown.\n"
			   	<< "Options are\n"
				<<	"== 'rk4'		(fourth order Runge-Kutta)\n"
				<<	"== 'lf'		(leapfrog)\n"
				<<	"== 'langevin'	(Langevin dynamics integration)\n"
				<<	"== 'vm'		(Vicsek model with angular noise)\n"
				<<	"== 'mc'		(Monte Carlo, Metropolis-Hastings. NOT IMPLEMENTED)\n"
				<<	"== Run will be aborted\n\n";
	   abort();
   }
}


void integrator::integrator_eq(const parameters& par){
	integrator_type_ = par.eq_integrator_type();
}
void integrator::integrator_sample(const parameters& par){
	integrator_type_ = par.sample_integrator_type();
}
void integrator::initialize_parameters(const parameters& par){
	kT_ = par.kT();
	dt_ = par.dt();
	g_ = par.dof();
	activity_ = par.activity();

	if ( (integrator_type_ == "nh") | (integrator_type_ == "np")){
		pi_ = par.nhnp_pi();
		Q_ = par.nhnp_Q();
		tau_ = par.nhnp_tau();
		if (integrator_type_ == "nh"){
			eta_ = par.nh_eta();
		} else if (integrator_type_ == "np"){
			s_ = par.np_s();
		}
		if ( Q_ <= 0){
			Q_ = g_ * kT_ * tau_ * tau_; 	// Formula by Evans and Holian, below eq (21),
										 	// J. Chem. Phys. 83, 4069 (1985)
											// https://aip.scitation.org/doi/pdf/10.1063/1.449071
		}
	} else if ( integrator_type_ == "langevin" ) {
		gamma_ld_om_ = par.gamma_ld_om();
		gamma_ld_p_ = par.gamma_ld_p();
	}
}



void integrator::initialize(const parameters& par, const group& G){
	initialize_parameters(par);
	if ( ( G.get_group_type() == "mxy" || G.get_group_type() == "fmxy" || G.get_group_type() == "xy" ) && integrator_type_ == "np"){
		group G_new = G;
		if ( G.get_group_type() == "mxy" || G.get_group_type() == "fmxy" ){
			G_new.fill_partition(); // Since G is const, it has to be copied here to update the partition
									// required for computing the interaction energy.
			if ( G.get_group_type() == "fmxy" ){
				G_new.generate_neighbor_list();
			}
		}
		H_0_ = G_new.calc_kinetic_energy() / (s_ * s_) + G_new.calc_interaction_energy() + .5*pi_*pi_/Q_ + g_ * par.kT() * std::log(s_);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
// SOLVERS
///////////////////////////////////////////////////////////////////////////////////////////
group integrator::integrate(const group& G, std::vector<double>& momdot) {
	if ( G.get_group_type() == "mxy" || G.get_group_type() == "fmxy" || G.get_group_type() == "xy" ) {
		group returngroup;
		if (integrator_type_ == "rk4"){
			returngroup = rk4(G);
		} else if (integrator_type_ == "lf"){
			returngroup = leapfrog(G,momdot);
// TODO
//		} else if (integrator_type_ == "mc"){
//			return mc(G);
		} else if (integrator_type_ == "nh"){
			returngroup = nh(G);
		} else if (integrator_type_ == "np"){
			returngroup = np(G);
		} else if (integrator_type_ == "langevin"){
			if (activity_ > 0) {
				returngroup = langevin_active(G,momdot,activity_);
			} else {
				returngroup = langevin(G,momdot);
			}
		} else {
			std::cout << "WARNING! Unknown integrator type " << integrator_type_ << " for group type " << G.get_group_type() << ".  Returning group" << std::endl;
			return G;
		}
		if ( (G.get_group_type() == "mxy") && (activity_ > 0) && (integrator_type_ != "langevin")) {
			add_activity(returngroup, activity_);
			returngroup.set_r_to_pbc();
			returngroup.fill_partition();
		}
		return returngroup;
	} else if ( G.get_group_type() == "vm" ||G.get_group_type() == "fvm") {
		if (integrator_type_ == "vm"){//<--
			return vm_rule(G);
		} else {
			std::cout << "WARNING! Unknown integrator type " << integrator_type_ << " for Vicsek model.  Returning group" << std::endl;
		}
	} else {
		std::cout << "WARNING! Unknown group type " << G.get_group_type() << ". Returning group" << std::endl;

	}
	return G;
}

group integrator::vm_rule(const group& G) {
	////////////////////////////////////////////////////////////
	// Three steps: Collision, noise and streaming
	////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	// Step A: Collision
	// Calculate average orientation around each particle.
	////////////////////////////////////////////////////////////
	group G_new = G;
	std::vector<double> distances; // Dummy vector, will not be filled.
	std::vector<int> nb;
	/* IMPLEMENTATION A: Only upper right neighbors.
	 * Every interaction is counted only once, but the process costs much
	 * more memory (a length N std::vector of topology::Vector2d has to
	 * be stored to keep track. For vectorial noise, one would also have to
	 * store a std::vector of int counting the number of neighbors.
	 */
	std::vector<topology::Vector2d> orientation_sum(G.get_N(),0);
	for (int i = 0; i < G.get_N(); i++){
		orientation_sum[i] += topology::vector_from_angle(G.get_theta(i));
		nb = G.get_neighbors(i,"ur",distances);
		for (size_t j = 0; j < nb.size(); j++){
			orientation_sum[i] += topology::vector_from_angle(G.get_theta(nb[j]));
			orientation_sum[nb[j]] += topology::vector_from_angle(G.get_theta(i));
		}

	}
	for (int i = 0; i < G.get_N(); i++){
		G_new.set_theta(atan2(orientation_sum[i].get_y(),orientation_sum[i].get_x()), i ); // Uses atan2 function to get the correct angle.
	}

	////////////////////////////////////////////////////////////
	// Step B: Noise
	// Adds noise to the angle
	////////////////////////////////////////////////////////////
	G_new.add_random_angle(G.get_vm_eta() / 2 * M_PI);
	G_new.set_theta_to_interval();

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Step C: Stream
	// Streams along the new orientation of the velocity. Only streams if model is mobile VM.
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( G.get_group_type() == "vm" ) {
		G_new.stream_along_spin(G.get_vm_v() * dt_);
		G_new.set_r_to_pbc();
		G_new.fill_partition();
	}
	return G_new;
}



group integrator::rk4(const group& G) {
	group G_new = G;
	if ( G.get_group_type() == "mxy" ){
		G_new.fill_partition();
	}
	group rkstep = G_new.time_derivative();
	group stepsum = rkstep * b_[0];
	for(int j=0; j<3; j++){
		G_new = G + rkstep*(dt_*a_[j]);
		if ( G.get_group_type() == "mxy" ){
			G_new.set_r_to_pbc();
			G_new.fill_partition();
		}
		rkstep = G_new.time_derivative();
		stepsum = stepsum + rkstep*b_[j+1];
	}
	G_new = G + stepsum * dt_;
	//////////////////////////////////////////////////////////////////////////////////////////
	// Refreshing of partition is not performed if activity is to be added in at
	// a later point
	//////////////////////////////////////////////////////////////////////////////////////////
	if ( G.get_group_type() == "mxy" && activity_ <= 0 ){
		G_new.set_r_to_pbc();
		G_new.fill_partition();
	}
	G_new.set_theta_to_interval();
	return G_new;
}



group integrator::leapfrog(const group& G, std::vector<double>& momdot) {
	group G_new = G;
	if (momdot.empty() || (G.get_group_type() == "mxy" && activity_ > 0) ){
		momdot = G_new.time_derivative_mom();
	}
	G_new.add_to_coord(G_new.time_derivative_coord(), dt_);
	G_new.add_to_coord_inertialscaling(momdot, .5 * dt_* dt_ );
	if ( G.get_group_type() == "mxy" ){
		G_new.set_r_to_pbc();
		G_new.fill_partition();
	}
	G_new.set_theta_to_interval();
	G_new.add_to_mom(momdot,.5 * dt_ );
	momdot = G_new.time_derivative_mom();
	G_new.add_to_mom(momdot,.5 * dt_);

	return G_new;
}




group integrator::leapfrog_active(const group& G, std::vector<double>& momdot, double activity) {
	group G_new = G;
	if (momdot.empty() || (G.get_group_type() == "mxy") ){
		momdot = G_new.time_derivative_mom();
	}
	G_new.add_to_coord(G_new.time_derivative_coord(), dt_);
	G_new.add_to_coord_inertialscaling(momdot, .5 * dt_* dt_ );
	add_activity(G_new,activity);
	if ( G.get_group_type() == "mxy" ){
		G_new.set_r_to_pbc();
		G_new.fill_partition();
	}
	G_new.set_theta_to_interval();
	G_new.add_to_mom(momdot,.5 * dt_ );
	momdot = G_new.time_derivative_mom();
	G_new.add_to_mom(momdot,.5 * dt_);
	return G_new;
}



group integrator::nh(const group& G) {
	// eq (6) to (16) of
	// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577381/
	// Kleinerman et al. J Chem Phy, 2008, Implementation of Nose-Hoover and Nose-Poincare [...]
	group G_new = G;
	G_new.fill_partition();
	double Ekin = G_new.calc_kinetic_energy();
	double s_hat;

	pi_ += .25 * dt_ * ( 2 * Ekin - g_ * kT_); // eq (6)
	s_hat = std::exp( - .5 * pi_ * dt_ / Q_); // eq (7)
	eta_ += .5 * pi_ * dt_ / Q_; // eq (8)
	pi_ += .25 * dt_ * ( 2 * Ekin * s_hat * s_hat - g_ * kT_); // eq (9). Correct?
	G_new.scale_mom(s_hat); // first part of (10)
	G_new.add_to_mom(vector_scale(G.time_derivative_mom(),.5*dt_)); // second part of (10)
	G_new.add_to_coord_inertialscaling(vector_scale(G_new.get_mom(), dt_ )); // eq (11)

	G_new.set_r_to_pbc();
	G_new.fill_partition();
	G_new.set_theta_to_interval();

	G_new.add_to_mom(vector_scale(G_new.time_derivative_mom(),.5*dt_)); // eq (12)
	Ekin = G_new.calc_kinetic_energy();
	pi_ += .25 * dt_ * ( 2 * Ekin - g_ * kT_); // eq (13)
	s_hat = std::exp( - .5 * pi_ * dt_ / Q_); // eq (14)
	eta_ += .5 * pi_ * dt_ / Q_; // eq (15)
	pi_ += .25 * dt_ * ( 2 * Ekin * s_hat * s_hat - g_ * kT_); // eq (16). cf eq (9)
	G_new.scale_mom(s_hat); // eq (17)

	return G_new;
}


group integrator::np(const group& G) {
	// eq (21) to (28) of
	// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577381/
	// Kleinerman et al. J Chem Phy, 2008, Implementation of Nose-Hoover and Nose-Poincare [...]

	group G_new = G;
	G_new.fill_partition();
	double interaction_energy = G.calc_interaction_energy();
	G_new.scale_mom(s_); // Transforming to tilde p = p * s.

	s_ *= std::pow(( 1 + .25 * pi_ * dt_ / Q_),2); // eq (21)
	pi_ /= ( 1 + .25 * pi_ * dt_ / Q_); // eq (22)
	G_new.fill_partition();
	G_new.add_to_mom(vector_scale(G_new.time_derivative_mom(),.5 * s_ * dt_)); // eq (23)
	G_new.add_to_coord_inertialscaling(vector_scale(G_new.get_mom(), dt_ / s_ )); // eq (24)
	G_new.set_r_to_pbc();
	G_new.set_theta_to_interval();
	G_new.fill_partition();
	G_new.add_to_mom(vector_scale(G_new.time_derivative_mom(),.5 * s_ * dt_)); // eq (25)
	pi_ += dt_ * ( G_new.calc_kinetic_energy() / ( s_ * s_ ) - .5 * (interaction_energy + G_new.calc_interaction_energy())
			- g_ * kT_ * (std::log(s_) + 1) - .5 * pi_ * pi_ / Q_ + H_0_ ); // eq (26)
	s_ *= std::pow(( 1 + .25 * pi_ * dt_ / Q_ ),2); // eq (27)
	pi_ /= ( 1 + .25 * pi_ * dt_ / Q_); // eq (28)

	G_new.scale_mom(1/s_); // Transforming back from \tilde p to p = \tilde p / s.

	return G_new;
}


// TODO mc for xy
/*
template<> xygroup integrator::mc(const xygroup& G){
	xygroup G_new = G;
	int index;
	for(int j=0; j < G.get_N(); j++){
		index = rand() % G.get_N();

	}
//	G_new = G + stepsum * dt_;
	G_new.set_r_to_pbc();
	G_new.set_theta_to_interval();
	G_new.fill_partition();
	return G_new;
}

// TODO mc for mxy
template<> mxygroup integrator::mc(const mxygroup& G){
	mxygroup G_new = G;
	G_new.fill_partition();
	int index;
	for(int j=0; j < G.get_N(); j++){
		index = rand() % G.get_N();

	}
//	G_new = G + stepsum * dt_;
	G_new.set_r_to_pbc();
	G_new.set_theta_to_interval();
	G_new.fill_partition();
	return G_new;
}
*/
//template<class GROUP>
//GROUP integrator::mc(const GROUP& G) {
//	GROUP G_new = G;
//	G_new.fill_partition();
//	GROUP rkstep = G_new.time_derivative();
//	GROUP stepsum = rkstep * b_[0];
//	for(int j=0; j<3; j++){
//		G_new = G + rkstep*(dt_*a_[j]);
//		G_new.set_r_to_pbc();
//		G_new.fill_partition();
//		rkstep = G_new.time_derivative();
//		stepsum = stepsum + rkstep*b_[j+1];
//	}
//	G_new = G + stepsum * dt_;
//	G_new.set_r_to_pbc();
//	G_new.fill_partition();
//	return G_new;
//}
//template xygroup integrator::mc(const xygroup& G);
//template mxygroup integrator::mc(const mxygroup& G);






group integrator::langevin(const group& G, std::vector<double>& momdot) {
	group G_new = G;
	if (momdot.empty() || (G.get_group_type() == "mxy" && activity_ > 0) ){
		momdot = G_new.time_derivative_mom();
	}

	// p(t+dt/2), (12.3a)
	G_new.add_to_mom(momdot,.5 * dt_ );

	// r(t+dt/2), (12.3b)
	G_new.add_to_coord(G_new.time_derivative_coord(), dt_);

	// p'(t+dt/2), (12.3c)
	std::default_random_engine generator;
	std::normal_distribution<double> nd(0.0,1.0);
	if ( gamma_ld_om_ > 0 || (G.get_group_type() == "mxy" && gamma_ld_p_ > 0)) {
		for (int i=0; i < G.get_N(); i++) {
			G_new.set_w(i,exp(-gamma_ld_om_ * dt_) * G_new.get_w(i) + sqrt(1 - exp(-2*gamma_ld_om_ * dt_)) * sqrt(G.get_I() * kT_te_ ) * nd(generator));
		}
		if ( G.get_group_type() == "mxy" ) {
			for (int i=0; i < G.get_N(); i++) {
				G_new.set_px(i,exp(-gamma_ld_p_ * dt_) * G_new.get_p(i).get_x() + sqrt(1 - exp(-2*gamma_ld_p_ * dt_)) * sqrt(G.get_m() * kT_r_ ) * nd(generator));
				G_new.set_py(i,exp(-gamma_ld_p_ * dt_) * G_new.get_p(i).get_y() + sqrt(1 - exp(-2*gamma_ld_p_ * dt_)) * sqrt(G.get_m() * kT_r_ ) * nd(generator));
			}
		}
	}

	// r(t+dt), (12.3d)
	G_new.add_to_coord(G_new.time_derivative_coord(), dt_);

	// handle coordinates and partition
	if ( G.get_group_type() == "mxy" ){
		G_new.set_r_to_pbc();
		G_new.fill_partition();
	}
	G_new.set_theta_to_interval();

	momdot = G_new.time_derivative_mom();

	// p(t+dt), (12.3e)
	G_new.add_to_mom(momdot,.5 * dt_);

	return G_new;
}



group integrator::langevin_active(const group& G, std::vector<double>& momdot, double activity) {
	group G_new = G;
	if (momdot.empty() || (G.get_group_type() == "mxy") ){
		momdot = G_new.time_derivative_mom();
	}

	// p(t+dt/2), (12.3a)
	G_new.add_to_mom(momdot,.5 * dt_ );

	// r(t+dt/2), (12.3b)
	G_new.add_to_coord(G_new.time_derivative_coord(), dt_);

	// Adding activity in (.5 * activity due to the half step)
	add_activity(G_new,.5*activity);


	// p'(t+dt/2), (12.3c)
	std::default_random_engine generator;
	std::normal_distribution<double> nd(0.0,1.0);
	if ( gamma_ld_om_ > 0 || (G.get_group_type() == "mxy" && gamma_ld_p_ > 0)) {
		for (int i=0; i < G.get_N(); i++) {
			G_new.set_w(i,exp(-gamma_ld_om_ * dt_) * G_new.get_w(i) + sqrt(1 - exp(-2*gamma_ld_om_ * dt_)) * sqrt(G.get_I() * kT_te_ ) * nd(generator));
		}
		if ( G.get_group_type() == "mxy" ) {
			for (int i=0; i < G.get_N(); i++) {
				G_new.set_px(i,exp(-gamma_ld_p_ * dt_) * G_new.get_p(i).get_x() + sqrt(1 - exp(-2*gamma_ld_p_ * dt_)) * sqrt(G.get_m() * kT_r_ ) * nd(generator));
				G_new.set_py(i,exp(-gamma_ld_p_ * dt_) * G_new.get_p(i).get_y() + sqrt(1 - exp(-2*gamma_ld_p_ * dt_)) * sqrt(G.get_m() * kT_r_ ) * nd(generator));
			}
		}
	}

	// r(t+dt), (12.3d)
	G_new.add_to_coord(G_new.time_derivative_coord(), dt_);

	// Adding activity in
	add_activity(G_new,.5*activity);

	// handle coordinates and partition
	if ( G.get_group_type() == "mxy" ){
		G_new.set_r_to_pbc();
		G_new.fill_partition();
	}
	G_new.set_theta_to_interval();

	momdot = G_new.time_derivative_mom();

	// p(t+dt), (12.3e)
	G_new.add_to_mom(momdot,.5 * dt_);

	return G_new;
}


///////////////////////////////////////////////////////////////////////////////////////////
// THERMOSTATS
///////////////////////////////////////////////////////////////////////////////////////////

double integrator::berendsen_thermostat(group& G, double T_desired, double tau) {
	if ( ( G.get_group_type() == "mxy" ) || ( G.get_group_type() == "xy" ) ){
		double T_cur = G.calc_temperature();
		G.scale_mom(sqrt(1 + dt_ / tau * (T_desired / T_cur - 1)));
		return T_cur;
	} else {
		return 0;
	}
}



void integrator::add_activity(group& G, double activity) const {
	G.stream_along_spin(activity * dt_);
}

