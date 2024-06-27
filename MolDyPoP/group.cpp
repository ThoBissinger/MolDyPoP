/*! \file group.cpp
 *
 *  \brief cpp-File to class declaration of group.
 *  Implements routines for the group.
 *
 *  Implements the functions declared in file group.h. Most code should be
 *  self-explanatory, see the documentation in group.h for an overview of
 *  what each function does.
 *
 *  \author Thomas Bissinger, additional contributions by Mathias Hoefler
 *
 *  \date Created: 2020-02-29 (full rewrite)
 *  \date Last Updated: 2023-07-23
 */

#include "group.h"
//#include "xyparticle.h"
//#include <iomanip>
//#include <math.h>


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// CLASS GROUP
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Constructor
group::group(const parameters& par) :
	group_type_(par.system()), N_(par.N()), sqrtN_(par.sqrtN()), L_(par.L()) {
	nb_rule_ = "ur";
	/*if ( group_type_ == "fvm" || group_type_ == "vm"){
		nb_rule_ = "all";
	} else {
		nb_rule_ = "ur";
	}*/
	nb_mult_factor_ = 1;

	nb_list_ = new neighbor_list;
	//neighbor_list nb_list;
	//nb_list_ = &nb_list;
	if ( group_type_ == "mxy" || group_type_ == "fmxy" || group_type_ == "xy") {
		J_ = par.J();
		I_ = par.I();

		if ( group_type_ == "mxy" ){
			U_ = par.U();
			m_ = par.m();

		}
	}
	if ( group_type_ == "mxy" || group_type_ == "fmxy" ||group_type_ == "vm" ) {
		cutoff_ = par.cutoff();
		partition_ = partition(N_,cutoff_,L_);
	} else if ( group_type_ == "xy" ||  group_type_ == "fvm" ) {
		lattice_type_ = par.lattice_type();
		if ( lattice_type_ != 's' && lattice_type_ != 't' && lattice_type_ != 'n' ){
			std::cerr << "  ERROR! Unknown lattice type!" << std::endl
								<< "  Unknown type:   " << lattice_type_ << std::endl
								<< "  Known types:    's'      square lattice" << std::endl
								<< "  Known types:    't'      trigonal lattice" << std::endl
								<< "  Known types:    'n'      no lattice (for mobile models)" << std::endl;
			exit(-1);
		}
	}
	if ( group_type_ == "vm" || group_type_ == "fvm" ){
		vm_eta_ = par.vm_eta();
		if ( group_type_ == "vm" ) {
			vm_v_ = par.vm_v();
		}
	}
	if ( group_type_ != "mxy" && group_type_ != "fmxy" && group_type_ != "xy" && group_type_ != "vm" && group_type_ != "fvm" ){
		std::cerr 	<< "  ERROR! Unknown group type!" << std::endl
					<< "         Unknown type:   " << group_type_ << std::endl
					<< "         Known types:    'mxy'    mobile xy model" << std::endl
					<< "                         'fmxy'   frozen mobile xy model" << std::endl
					<< "                         'xy'     static xy model" << std::endl
					<< "                         'vm'     Vicsek model" << std::endl
					<< "                         'fvm'    frozen (static) Vicsek model" << std::endl;
		exit(-1);
	}


}




group::group(const int N, const std::string group_type) : group_type_(group_type), N_(N), sqrtN_((int)sqrt(N_)) {
	// neighbor_list nb_list;
	// nb_list_ = &nb_list;
	// nb_mult_factor_ = 0;
};
//////////////////////////////////////////////////////////////////////////////
// Initialization functions
void group::clear() {
	theta_.clear();
	w_.clear();
	r_.clear();
	p_.clear();
	partition_.clear();
}



void group::initialize(const parameters& par){
	// Setting up the group according to parameters.
	if ( par.init_kT() >= 0){
		initialize_random(par.init_kT());
	} else {
		initialize_random(par.kT());
	}
	if ( par.init_mode() == "file" ){
		read_from_snapshot(par.init_file());
		mom_to_zero();
		// Lattice-based models still have to initialize their positions to the lattice
		// if ( par.system() == "mxy" || par.system() == "vm" || par.system() == "fmxy" ) {
		// 	fill_partition();
		//	if ( par.system() == "fmxy" ){
		//		generate_neighbor_list();
		//	}
		//}
		std::cout << "  ++ Reading initial coordinates from file " << par.init_file() << std::endl;
		if ( ( par.system() == "mxy" || par.system() == "xy" || par.system() == "fmxy" ) && par.init_kT() >= 0){
			std::cout << "  ++ Redrawing random temperatures. Target temperature " << par.init_kT() << std::endl;
			set_temperature(par.init_kT());
		}
	} else if ( par.init_mode() == "scale" ){
		// loads from subgroup
		std::cout << "  ++ Initial coordinates taken from scaling the file " << par.init_file() << std::endl;
		scale_from_subgroup(par.init_file());
		// Lattice-based models still have to initialize their positions to the lattice
		mom_to_zero();
		// adds random angle (if it was set)
		if ( par.init_random_angle() != 0){
			std::cout << "  ++ Adding noise to all angles. Magnitude " << par.init_random_angle() << std::endl;
			add_random_angle(par.init_random_angle());
		}
		// adds random displacement (if it was set)
		if ( par.init_random_displacement() != 0 && ( par.system() == "mxy" || par.system() == "vm" ) ){
			std::cout << "  ++ Adding noise to all positons. Magnitude " << par.init_random_displacement() << std::endl;
			add_random_displacement(par.init_random_displacement());
			set_r_to_pbc();
		}
		// randomizes velocities for new temperature (if it was set)
		if ( ( par.system() == "mxy" || par.system() == "xy" || par.system() == "fmxy" ) && par.init_kT() >= 0){
			std::cout << "  ++ Redrawing random temperatures. Target temperature " << par.init_kT() << std::endl;
			set_temperature(par.init_kT());
		}
	} else if ( par.init_mode() == "aligned" ){
		std::cout << "  ++ Aligned spin initialization" << std::endl;
		set_all_theta(random_angle());
	} else if ( par.init_mode() == "random" ){
		std::cout << "  ++ Random particle initialization" << std::endl;
	} else {
		std::cout << "!! WARNING: Initialization mode " << par.init_mode() << " unknown. Using random data." << std::endl;
		std::cerr << "!! WARNING: Initialization mode " << par.init_mode() << " unknown. Using random data." << std::endl;
	}
	if ( par.system() == "mxy" || par.system() == "vm" || par.system() == "fmxy" ) {
		fill_partition();
	}
	if ( par.system() == "xy" || par.system() == "fvm" || par.system() == "fmxy" ){
		generate_neighbor_list();
	}
}

void group::initialize_random(double kT) {
	randomize_particles(kT);
	mom_to_zero();
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {
		fill_partition();
	}
}


// Careful, no neighbors are prescribed.
void group::randomize_particles(double kT) {
	clear();
	if (group_type_ == "mxy") {
	  for (int i=0; i< N_; i++){
			theta_.push_back(random_angle());
			w_.push_back(random_boltzmann_double(I_ * kT));
			r_.push_back(topology::random_vector(L_));
			p_.push_back(topology::random_gaussian_vector(m_ * kT));
	  }
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ) {
		for (int i=0; i< N_; i++){
			theta_.push_back(random_angle());
			w_.push_back(random_boltzmann_double(I_ * kT));
		}
		if ( group_type_ == "xy") {
			r_to_lattice();
		}
	} else if (group_type_ == "vm") {
		for (int i=0; i< N_; i++){
			theta_.push_back(random_angle());
			r_.push_back(topology::random_vector(L_));
		}
	} else if (group_type_ == "fvm") {
		for (int i=0; i< N_; i++){
			theta_.push_back(random_angle());
		}
		r_to_lattice();
	}
}


void group::mom_to_zero() {
	// For spatial momenta (only required if defined)
	if ( p_.size() != 0 ) {
		topology::Vector2d p_sum = sum_p() / N_;
		for (int i = 0; i < N_; i++){
			p_[i] -= p_sum;
		}
	}
	// For spin momenta (only required if defined)
	if ( w_.size() != 0 ) {
		double w_sum = sum_w() / N_;
		for (int i = 0; i < N_; i++){
			w_[i] -= w_sum;
		}
	}
}

void group::r_to_squarelattice() {
	r_.clear();
	for (int i=0; i< N_; i++){
		r_.push_back(topology::vector_on_squarelattice(i,sqrtN_,sqrtN_,sqrtN_/L_.get_x()));
	}
}
void group::r_to_trigonallattice() {
	r_.clear();
	for (int i=0; i< N_; i++){
		r_.push_back(topology::vector_on_squarelattice(i,sqrtN_,sqrtN_,sqrtN_/L_.get_x()));
	}
}
void group::r_to_lattice() {
	if ( lattice_type_ == 't' ){
		r_to_trigonallattice();
	} else if ( lattice_type_ == 's' ){
		r_to_squarelattice();
	} else if ( lattice_type_ == 'n' ){
		;
	} else {
		std::cerr 	<< "! WARNING In function group::initialize::group_type_ '" << group_type_ << "' requires a valid lattice_type_." << std::endl
					<< "!         lattice_type_ ='" << lattice_type_ << "' unknown. " << std::endl
					<< "!         Known types:     't'     trigonal lattice" << std::endl
					<< "!                          's'     square lattice" << std::endl
					<< "!         Using trigonal lattice." << std::endl;
		r_to_trigonallattice();
	}
}

// Sets everything to zero.
void group::initialize_zero() {
	clear();
	if (group_type_ == "mxy") {
		for (int i=0; i< N_; i++){
			theta_.push_back(0);
			w_.push_back(0);
			r_.push_back(0);
			p_.push_back(0);
		}
	} else if (group_type_ == "xy" || group_type_ == "fmxy") {
		for (int i=0; i< N_; i++){
			theta_.push_back(0);
			w_.push_back(0);
		}
	} else if (group_type_ == "vm") {
		for (int i=0; i< N_; i++){
			theta_.push_back(0);
			r_.push_back(0);
		}
	}
}

void group::fill_partition() {
	if (group_type_ == "mxy" || group_type_ == "vm" || group_type_ == "fmxy") {
		partition_.fill(r_);
	}
}


void group::read_from_snapshot(std::string snapshotname){
	clear();
	std::ifstream snapfile(snapshotname, std::ios::in);
	std::string line;
	double theta, w, rx, ry, px, py;

	if (group_type_ == "mxy" || group_type_ == "fmxy") {
		// File format for mxy model:
		// t t t ...
		// w w w ...
		// rx ry rx ry rx ry ...
		// px py px py px py ...
		// t: theta, angles
		// w: omega, spin momenta
		// rx, ry: x and y positions of particle
		// px, py: x and y directions of momenta (ignored in case of fmxy model)
		getline(snapfile, line);
		std::istringstream iss(line);
		while(iss >> theta){
			theta_.push_back(theta);
		}

		getline(snapfile, line);
		iss=std::istringstream(line);
		while(iss >> w){
			w_.push_back(w);
		}

		getline(snapfile, line);
		iss=std::istringstream(line);
		while(iss >> rx >> ry){
			r_.push_back(topology::Vector2d(rx,ry));
		}

		if ( group_type_ == "mxy" ){
			getline(snapfile, line);
			iss=std::istringstream(line);
			while(iss >> px >> py){
				p_.push_back(topology::Vector2d(px,py));
			}
		}

	} else if (group_type_ == "xy") {
		// File format for xy model:
		// t t t ...
		// w w w ...
		// t: theta, angles
		// w: omega, spin momenta
		getline(snapfile, line);
		std::istringstream iss(line);
		while(iss >> theta){
			theta_.push_back(theta);
		}

		getline(snapfile, line);
		iss=std::istringstream(line);
		while(iss >> w){
			w_.push_back(w);
		}
	} else if (group_type_ == "vm") {
		// File format for Vicsek model:
		// t t t ...
		// rx ry rx ry rx ry ...
		// t: theta, angles/directions
		// rx, ry: x and y positions of particle
		getline(snapfile, line);
		std::istringstream iss(line);
		while(iss >> theta){
			theta_.push_back(theta);
		}

		getline(snapfile, line);
		iss=std::istringstream(line);
		while(iss >> rx >> ry){
			r_.push_back(topology::Vector2d(rx,ry));
		}

	}else if (group_type_ == "fvm") {
		// File format for frozen Vicsek model:
		// t t t ...
		// t: theta, angles
		// rx, ry: x and y positions of particle
		getline(snapfile, line);
		std::istringstream iss(line);
		while(iss >> theta){
			theta_.push_back(theta);
		}

	}
	snapfile.close();
	// Lattice-based groups need to set their r parameters (independent of snapshot)
	if (group_type_ == "xy" || group_type_ == "fvm") {
		r_to_lattice();
	}

}

void group::scale_from_subgroup(const group& G){
	double scalefac = 1;
	double scalepow = 0;
	while ( scalefac * G.get_sqrtN() < sqrtN_ ){
		scalefac *= 2;
		scalepow++;
	}
	if ( scalefac * G.get_sqrtN() != sqrtN_ ){
		std::cerr << "Error: subgroup dimension not compatible with group dimension" << std::endl;
	} else {
		if ( group_type_ == "mxy" || group_type_ == "fmxy" || group_type_ == "vm" ) {
			int oldgroup_sqrtN = G.get_sqrtN();
			int oldgroup_N = G.get_N();

			// Scale the positions to have uniform coverage of the box (some sizes may be incommensurable)
			double max_x = 0, max_y = 0;
			for (int i=0; i < oldgroup_N; i++){
				if ( G.get_r(i).get_x() > max_x ){
					max_x = G.get_r(i).get_x();
				}
				if ( G.get_r(i).get_y() > max_y ){
					max_y = G.get_r(i).get_y();
				}
			}
			// TODO: Scale from boxsize of old group
			double x_scale = L_.get_x() / max_x * oldgroup_sqrtN / sqrtN_;
			double y_scale = L_.get_y() / max_y * oldgroup_sqrtN / sqrtN_;
			for (int i=0; i < oldgroup_N; i++){
				r_[i] = topology::Vector2d(r_[i].get_x() * x_scale,r_[i].get_y() * y_scale);
			}

			int new_index, old_index;
			// More details on the idea behind the loop:
			// xygroup::scale_from_subgroup(std::string snapshotname) in xygroup.cpp
			for (int i=0; i < oldgroup_sqrtN; i++){
				for (int j = oldgroup_sqrtN - 1; j >= 0; j--){ // Top-Down is important to not overwrite data.
					new_index = i + j * sqrtN_; // ith column, jth row (from the bottom), at least in my way of viewing it. In the NEW group
					old_index = i + j * oldgroup_sqrtN; // ith column, jth row (from the bottom), at least in my way of viewing it. In the OLD group.
					theta_[new_index] = G.get_theta(old_index);
					r_[new_index] = G.get_r(old_index);
					if ( group_type_ == "mxy" ) {
						// Only defined for mxy group, vm doesn't have this
						w_[new_index] = G.get_w(old_index);
						p_[new_index] = G.get_p(old_index);
					}
					for (int k=0; k < scalefac; k++){
						for (int l=0; l < scalefac; l++){
							if (k + l > 0){ // to ignore redundancy
								// Copies particles
								theta_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = G.get_theta(new_index);
								r_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = G.get_r(new_index);
								if ( group_type_ == "mxy" ) {
									// Only defined for mxy group, vm doesn't have this
									w_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = G.get_w(new_index);
									p_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = G.get_p(new_index);
								}
								// particles_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = particles_[new_index];
								// Shifts positions
								r_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] += topology::Vector2d(k * L_.get_x(), l * L_.get_y()) / scalefac;
								//particles_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN].add_to_r(topology::Vector2d(k * L_.get_x(), l * L_.get_y()) / scalefac );
							}
						}
					}
				}
			}
		// In case the boundaries are not at the desired location, positions are set to pbc.
		// Should not be necessary when scaling the positions, but for safety...
		set_r_to_pbc();
		} else if ( group_type_ == "xy" || group_type_ == "fvm" ) {
			r_to_lattice();
			// TODO
			// Old code had it.
		} else {
			std::cerr 	<< "! WARNING Function group::scale_from_subgroup not defined for group_type_ = " << group_type_ << std::endl
						<< "!         Doing nothing." << std::endl;
		}
	}
}

void group::scale_from_subgroup(std::string snapshotname){
	std::ifstream snapfile(snapshotname, std::ios::in);
	std::string line;
	double theta, w, rx, ry, px, py;

	////////////////////////////////////////////////////////////////////////
	// Assigning variables from smaller group. Like read_from_snapshot.
	////////////////////////////////////////////////////////////////////////
	// theta always defined
	getline(snapfile, line);
	std::istringstream iss(line);
	int i = 0;
	while(iss >> theta){
		theta_[i] = theta;
		i++;
	}
	int oldgroup_N = i;
	int oldgroup_sqrtN = std::sqrt(oldgroup_N);

	if ( group_type_ == "mxy" || group_type_ == "xy" || group_type_ == "fmxy" ) {
		// w only defined in mxy and xy model
		getline(snapfile, line);
		iss=std::istringstream(line);
		i = 0;
		while(iss >> w){
			w_[i] = w;
			i++;
		}
	}
	double max_x = 0, max_y = 0;
	if ( group_type_ == "mxy" || group_type_ == "fmxy" || group_type_ == "vm" ) {
		// r only defined for mxy and Vicsek model
		getline(snapfile, line);
		iss=std::istringstream(line);
		i = 0;
		while(iss >> rx >> ry){
			r_[i] = topology::Vector2d(rx,ry);
			if ( rx > max_x ){
				max_x = rx;
			}
			if ( ry > max_y ){
				max_y = ry;
			}
			i++;
		}
	}

	if ( group_type_ == "mxy" ) {
		// p only defined for mxy model
		getline(snapfile, line);
		iss=std::istringstream(line);
		i = 0;
		while(iss >> px >> py){
			p_[i] = topology::Vector2d(px,py);
			i++;
		}
	}
	snapfile.close();


	////////////////////////////////////////////////////////////////////////
	// Assigning the rest. Like scale_from_subgroup.
	////////////////////////////////////////////////////////////////////////
	double scalefac = 1;
	double scalepow = 0;
	while ( scalefac * oldgroup_sqrtN < sqrtN_ ){
		scalefac *= 2;
		scalepow++;
	}
	if ( scalefac * oldgroup_sqrtN != sqrtN_ ){
		std::cerr << "Error: subgroup dimension not compatible with group dimension" << std::endl;
		std::cerr << "       old sqrtN = " << oldgroup_sqrtN << ", sqrtN = " << sqrtN_ << std::endl;
		std::cerr << "       scalefac  = " << scalefac << ", sclaepow = " << scalepow << std::endl;
	} else {
		if ( group_type_ == "mxy" || group_type_ == "vm" ) {
			// Scale the positions to have uniform coverage of the box (some sizes may be incommensurable)
			double x_scale = L_.get_x() / max_x * oldgroup_sqrtN / sqrtN_;
			double y_scale = L_.get_y() / max_y * oldgroup_sqrtN / sqrtN_;
			for (int i=0; i < oldgroup_N; i++){
				r_[i] = topology::Vector2d(r_[i].get_x() * x_scale,r_[i].get_y() * y_scale);
			}

			// More details on the idea behind the loop:
			// xygroup::scale_from_subgroup(std::string snapshotname) in xygroup.cpp
			int new_index, old_index;
			for (int i=0; i < oldgroup_sqrtN; i++){
				for (int j = oldgroup_sqrtN - 1; j >= 0; j--){ // Top-Down is important to not overwrite data.
					new_index = i + j * sqrtN_; // ith column, jth row (from the bottom), at least in my way of viewing it. In the NEW group
					old_index = i + j * oldgroup_sqrtN; // ith column, jth row (from the bottom), at least in my way of viewing it. In the OLD group.
					theta_[new_index] = theta_[old_index];
					r_[new_index] = r_[old_index];
					if ( group_type_ == "mxy" ) {
						// Only defined for mxy group, vm doesn't have this
						p_[new_index] = p_[old_index];
						w_[new_index] = w_[old_index];
					}
					for (int k=0; k < scalefac; k++){
						for (int l=0; l < scalefac; l++){
							if (k + l > 0){ // to ignore redundancy
								// Copies particles
								theta_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = theta_[new_index];
								r_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = r_[new_index];
								if ( group_type_ == "mxy" ) {
									// Only defined for mxy group, vm doesn't have this
									w_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = w_[new_index];
									p_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = p_[new_index];
								}
								// particles_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] = particles_[new_index];
								// Shifts positions
								r_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN] += topology::Vector2d(k * L_.get_x(), l * L_.get_y()) / scalefac;
								//particles_[new_index + (k + l  * sqrtN_) * oldgroup_sqrtN].add_to_r(topology::Vector2d(k * L_.get_x(), l * L_.get_y()) / scalefac );// Copies particles
							}
						}
					}
				}
			}
		// In case the boundaries are not at the desired location, positions are set to pbc.
		// Should not be necessary when scaling the positions, but for safety...
		set_r_to_pbc();
		} else if ( group_type_ == "xy" || group_type_ == "fvm" ) {
			r_to_lattice();
			// TODO
		} else {
			std::cerr 	<< "! WARNING Function group::scale_from_subgroup not defined for group_type_ = " << group_type_ << std::endl
						<< "!         Doing nothing." << std::endl;
		}
	}

}




//////////////////////////////////////////////////////////////////////////////
// Print functions
//////////////////////////////////////////////////////////////////////////////
void group::print_group(std::ofstream &outputfile) const {
	if (group_type_ == "mxy" || group_type_ == "fmxy" ) {
		for(int i = 0; i<N_; i++){
			outputfile << theta_[i] << ' ';
		}
		outputfile << std::endl;
		for(int i = 0; i<N_; i++){
			outputfile << w_[i] << ' ';
		}
		outputfile << std::endl;
		for(int i = 0; i<N_; i++){
			outputfile << r_[i].get_x() << ' ' << r_[i].get_y() << ' ';
		}
		outputfile << std::endl;
		if ( group_type_ == "mxy" ){
			for(int i = 0; i<N_; i++){
				outputfile << p_[i].get_x() << ' ' << p_[i].get_y() << ' ';
			}
			outputfile << std::endl;
		}
	} else if (group_type_ == "xy") {
		for(int i = 0; i<N_; i++){
			outputfile << theta_[i] << ' ';
		}
		outputfile << std::endl;
		for(int i = 0; i<N_; i++){
			outputfile << w_[i] << ' ';
		}
		outputfile << std::endl;
	} else if (group_type_ == "vm") {
		for(int i = 0; i<N_; i++){
			outputfile << theta_[i] << ' ';
		}
		outputfile << std::endl;
		for(int i = 0; i<N_; i++){
			outputfile << r_[i].get_x() << ' ' << r_[i].get_y() << ' ';
		}
		outputfile << std::endl;
	} else if (group_type_ == "fvm") {
		for(int i = 0; i<N_; i++){
			outputfile << theta_[i] << ' ';
		}
		outputfile << std::endl;
	}
}


void group::print_r(std::ofstream &outputfile) const {
	if ( ! r_.empty() ){
		for(int i = 0; i<N_; i++){
			outputfile << r_[i].get_x() << ' ' << r_[i].get_y() << ' ';
		}
		outputfile << std::endl;
	} else {
		std::cerr << "! WARNING: function group::print_r called but member variable r_ is empty !" << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// get functions
//////////////////////////////////////////////////////////////////////////////
// TODO: ERROR HANDLING
/* std::vector<double> group::get_theta() const{
	if ( ! theta_.empty() ){
		std::vector<double> theta(N_);
		for (int i=0; i < N_; i++){
			theta[i] = theta_[i];
		}
		return theta;
	} else {
		std::cerr << "! WARNING: function group::get_theta called but member variable theta_ is empty. Returning empty vector !" << std::endl;
		return {};
	}
}


std::vector<double> group::get_w() const{
	if ( ! theta_.empty() ){
		std::vector<double> w(N_);
		for (int i=0; i < N_; i++){
			w[i] = w_[i];
		}
		return w;
	} else {
		std::cerr << "! WARNING: function group::get_w called but member variable w_ is empty. Returning empty vector !" << std::endl;
		return {};
	}
}
std::vector<topology::Vector2d> group::get_r() const {
	if ( ! r_.empty() ){
		std::vector<double> r(N_);
		for (int i=0; i < N_; i++){
			r[i] = r_[i];
		}
		return r;
	} else {
		std::cerr << "! WARNING: function group::get_r called but member variable r_ is empty. Returning empty vector !" << std::endl;
		return {};
	}
}
std::vector<topology::Vector2d> group::get_p() const {
	if ( ! p_.empty() ){
		std::vector<double> p(N_);
		for (int i=0; i < N_; i++){
			p[i] = p_[i];
		}
		return p;
	} else {
		std::cerr << "! WARNING: function group::get_p called but member variable p_ is empty. Returning empty vector !" << std::endl;
		return {};
	}
}
*/

// TODO: Error handling.
std::vector<double> group::get_coord() const{
	std::vector<double> coord;
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {
		coord = std::vector<double>(3 * N_);
		for (int i=0; i < N_; i++){
			coord[i] = theta_[i];
			coord[N_ + i] = r_[i].get_x();
			coord[2 * N_ + i] = r_[i].get_y();
		}
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" || group_type_ == "fvm" ) {
		coord = std::vector<double>(N_);
		for (int i=0; i < N_; i++){
			coord[i] = theta_[i];
		}
	}
	return coord;
}



std::vector<double> group::get_mom() const{
	if ( group_type_ == "mxy"  ) {
		std::vector<double> mom(3 * N_);
		for (int i=0; i < N_; i++){
			mom[i] = w_[i];
			mom[N_ + i] = p_[i].get_x();
			mom[2 * N_ + i] = p_[i].get_y();
		}
		return mom;
	} else if ( group_type_ == "xy"  ) {
		std::vector<double> mom( N_);
		for (int i=0; i < N_; i++){
			mom[i] = w_[i];
		}
		return mom;
	} else if ( group_type_ == "vm"  ) {
		return {};
	} else {
		return {};
	}
}


std::vector<int> group::get_neighbors(int i, std::string cellselect, std::vector<double>& distances ) const{
	std::vector<int> nb;
	distances.clear();
	if ( ! nb_list_->is_empty() ){
		// distances are not assigned for lattice-based models, since all distances are unity there
		if ( group_type_ == "mxy" || group_type_ == "fmxy" ||group_type_ == "vm" ) {
			distances = nb_list_->get_dist(i);
		}
		return nb_list_->get_neighbors(i);
	} else if ( group_type_ == "mxy" || group_type_ == "fmxy" ||group_type_ == "vm" ) {
		std::vector<int> nb_cells, nb_in_cell;
		std::vector<double> distances_in_cell;
		int cell_index = partition_.find_cell(r_[i]);
		double cutoffsquared = cutoff_ * cutoff_;
		std::vector<topology::Vector2d> shifts;
		if ( cellselect == "ur" ) {
			///////////////////////////////////////////////////////////////
			// mxy or Vicsek model, above right neighbors
			///////////////////////////////////////////////////////////////
			// find neighboring cells (includes also the particle's cell). Excludes cells with no neighbors and unreachable cells.
			nb_cells = partition_.nb_cells_ur(r_[i],cutoffsquared,shifts);
	//		nb_cells = partition_.nb_cells_ur(partition_.find_cell(r_[i]),shifts);

			// checks for neighbor particles in the neighboring cells
			for (std::size_t j = 0; j < nb_cells.size(); j++ ){
				if (nb_cells[j] == cell_index){ // different rule for calculating neighbors in own cell.
					nb_in_cell = partition_.nb_in_cell_index_above(nb_cells[0], r_[i],
						cutoffsquared,r_,distances_in_cell,shifts[0]); // find neighbors in current cell
				} else { // calculating neighbors in adjacent cells
					nb_in_cell = partition_.nb_in_cell_index(nb_cells[j], r_[i],
						cutoffsquared,r_,distances_in_cell,shifts[j]); // find neighbors in cell
				}
				nb.insert(nb.end(),nb_in_cell.begin(),nb_in_cell.end());
				distances.insert(distances.end(),distances_in_cell.begin(),distances_in_cell.end());
			}
		} else if ( cellselect == "all" ){
			///////////////////////////////////////////////////////////////
			// mxy or Vicsek model, all neighbors
			///////////////////////////////////////////////////////////////
			nb_cells = partition_.nb_cells_all(partition_.find_cell(r_[i]),shifts);
			for (std::size_t j = 0; j < nb_cells.size(); j++ ){
				nb_in_cell = partition_.nb_in_cell_index(nb_cells[j], r_[i],
						cutoffsquared,r_,distances_in_cell,shifts[j]); // find neighbors in cell
				nb.insert(nb.end(),nb_in_cell.begin(),nb_in_cell.end());
				distances.insert(distances.end(),distances_in_cell.begin(),distances_in_cell.end());
			}
		} else if ( cellselect == "bruteforce" ){
			///////////////////////////////////////////////////////////////
			// mxy or Vicsek model, brute force neighbor calculation
			///////////////////////////////////////////////////////////////
			for (int j = 1; j < N_; j++){
				if ( (i != j) && (periodic_distance_squared(i,j) < cutoffsquared) ){
					nb.push_back(j);
					distances.push_back(periodic_distance(i,j));
				}
			}
		} else {
			std::cerr << "!! ERROR: In call to group::get_neighbors:\n"
								<< "!!        Incorrect nb_rule option " << cellselect << std::endl;
			exit(-1);
		}
	} else if ( group_type_ == "xy" || group_type_ == "fvm"  ) {
		if ( lattice_type_ == 's' ) {
			if ( cellselect == "ur" ){
				///////////////////////////////////////////////////////////////
				// xy model, square lattice, above right neighbors
				///////////////////////////////////////////////////////////////
				if ( (i + 1) > N_ - sqrtN_) { // upper neighbor
					nb.push_back(i - N_ + sqrtN_);
				} else {
					nb.push_back(i + sqrtN_);
				}
				if ((i + 1) % sqrtN_ == 0) { // right neighbor
					nb.push_back(i - sqrtN_ + 1);
				} else {
					nb.push_back(i + 1);
				}
			} else if ( cellselect == "all" ){
				///////////////////////////////////////////////////////////////
				// xy model, square lattice, all neighbors
				///////////////////////////////////////////////////////////////
				if ( (i + 1) > N_ - sqrtN_) { // upper neighbor
					nb.push_back(i - N_ + sqrtN_);
				} else {
					nb.push_back(i + sqrtN_);
				}
				if ((i + 1) % sqrtN_ == 0) { // right neighbor
					nb.push_back(i - sqrtN_ + 1);
				} else {
					nb.push_back(i + 1);
				}
				if (i < sqrtN_) { // lower neighbor
					nb.push_back(i + N_ - sqrtN_);
				} else {
					nb.push_back(i - sqrtN_);
				}
				if (i % sqrtN_ == 0) { // left neighbor
					nb.push_back(i + sqrtN_ - 1);
				} else {
					nb.push_back(i - 1);
				}
			}
		} else if ( lattice_type_ == 't' ) {
			int oddcheck = 0;
			// This checks whether the index is in an even or odd line. In a trigonal lattice, lines are
			// offset:
			//   x x x    2*sqrtN
			//    o o o   1*sqrtN
			//   x x x    0*sqrtN
			// etc. Therefore, boundary conditions depend on whether we are in an even or odd line.
			int x = i % sqrtN_;
			int y = (i - x) / sqrtN_;
			if ((y % 2) == 1) {
				oddcheck = 1;
			}
			int upshift = 0;
			int leftshift = 0;
			int rightshift = 0;
			int downshift = 0;
			if ( cellselect == "ur" ){
				///////////////////////////////////////////////////////////////
				// xy model, trigonal lattice, above right neighbors
				///////////////////////////////////////////////////////////////
				// Checks whether we are at the upper boundary, then right or left boundary
				// Lower boundary is not needed (half lattice)
				if (y == sqrtN_ - 1) { // upper boundary
					upshift = 1;
				}
				if (x == sqrtN_ - 1) { //right boundary
					rightshift = 1;
				} else if (x == 0) { // left boundary
					leftshift = 1;
				}
				nb.push_back(
						i - N_ * upshift + (1 - oddcheck) * (sqrtN_ * leftshift - 1)
								+ sqrtN_); // up left. Offset for even lines.
				nb.push_back(
						i - N_ * upshift + (oddcheck) * (-sqrtN_ * rightshift + 1)
								+ sqrtN_); // up right. Offset for odd lines.
				nb.push_back(i - sqrtN_ * rightshift + 1); // right. Same for even and odd lines

			} else if ( cellselect == "all" || cellselect == "bruteforce" ){
				///////////////////////////////////////////////////////////////
				// xy model, trigonal lattice, all neighbors
				///////////////////////////////////////////////////////////////
				// Checks whether we are at the upper or lower boundary, then right or left boundary
				if (y == sqrtN_ - 1) { // upper boundary
					upshift = 1;
				} else if (y == 0) { // lower boundary
					downshift = 1;
				}
				if (x == sqrtN_ - 1) { //right boundary
					rightshift = 1;
				} else if (x == 0) { // left boundary
					leftshift = 1;
				}

				nb.push_back(
						i - N_ * upshift + (1 - oddcheck) * (sqrtN_ * leftshift - 1)
								+ sqrtN_); // up left. Offset for even lines.
				nb.push_back(
						i - N_ * upshift + (oddcheck) * (-sqrtN_ * rightshift + 1)
								+ sqrtN_); // up right. Offset for odd lines.
				nb.push_back(i - sqrtN_ * rightshift + 1); // right. Same for even and odd lines
				nb.push_back(
						i + N_ * downshift + (oddcheck) * (-sqrtN_ * rightshift + 1)
								- sqrtN_); // down right. Offset for odd lines.
				nb.push_back(
						i + N_ * downshift + (1 - oddcheck) * (sqrtN_ * leftshift - 1)
								- sqrtN_); // down left. Offset for even lines.
				nb.push_back(i + sqrtN_ * leftshift - 1); // left. Same for even and odd lines
			}
		}
	}
	return nb;
}

void group::generate_neighbor_list(){
	free(nb_list_);
	nb_list_ = new neighbor_list(*this);
	//neighbor_list nb_list(*this);
	//nb_list_ = &nb_list;
}

//////////////////////////////////////////////////////////////////////////////
// Functions for setting values
//////////////////////////////////////////////////////////////////////////////
void  group::set_theta(double theta, int i) {
	theta_[i] = theta;
}

void  group::set_w(double w, int i) {
	w_[i] = w;
}
void  group::set_r(topology::Vector2d r, int i) {
	r_[i] = r;

}
void  group::set_rx(double x, int i) {
	r_[i].set_x(x);
}
void  group::set_ry(double y, int i) {
	r_[i].set_y(y);
}

void  group::set_p(topology::Vector2d p, int i) {
	p_[i] = p;
}
void  group::set_px(double px, int i) {
	p_[i].set_x(px);
}
void  group::set_py(double py, int i) {
	p_[i].set_y(py);
}


void group::set_all_theta(double theta){
	for(int i=0; i<N_; i++){
		theta_[i] = theta;
	}
}

void group::set_all_w(double w){
	for(int i=0; i<N_; i++){
		w_[i] = w;
	}
}

void group::set_all_p(topology::Vector2d p){
	for(int i=0; i<N_; i++){
		p_[i] = p;
	}
}


void  group::set_temperature(double kT, int i) {
	if ( ! w_.empty() ) {
		w_[i] = random_boltzmann_double(I_ * kT);
	}
	if ( ! p_.empty() ) {
		p_[i] = topology::random_gaussian_vector(m_ * kT);
	}
}

void group::set_temperature(double kT){
	if (kT > 0){
		for(int i=0; i<N_; i++){
			set_temperature(kT, i);
		}
	} else {
		set_all_w(0);
		set_all_p(0);
	}
}



void  group::set_temperature_p(double kT, int i) {
	if ( ! p_.empty() ) {
		p_[i] = topology::random_gaussian_vector(m_ * kT);
	}
}

void group::set_temperature_p(double kT){
	if (kT > 0){
		for(int i=0; i<N_; i++){
			set_temperature_p(kT, i);
		}
	} else {
		set_all_p(0);
	}
}


void  group::set_temperature_w(double kT, int i) {
	if ( ! w_.empty() ) {
		w_[i] = random_boltzmann_double(I_ * kT);
	}
}

void group::set_temperature_w(double kT){
	if (kT > 0){
		for(int i=0; i<N_; i++){
			set_temperature_w(kT, i);
		}
	} else {
		set_all_w(0);
	}
}

void group::set_theta_to_interval(){
	for(size_t i=0; i < theta_.size(); i++){
		theta_[i] = mod(theta_[i],-M_PI,M_PI);
	}
}

void group::set_r_to_pbc(){
	if ( group_type_ == "mxy" || group_type_ == "vm" ){
		for(size_t i=0; i < r_.size(); i++){
			r_[i].periodic_box(L_);
		}
	}
}

void group::scale_mom(double a){
	for(size_t i=0; i < w_.size(); i++){
		w_[i] *= a;
	}
	for(size_t i=0; i < p_.size(); i++){
		p_[i] *= a;
	}
}

void group::add_to_theta(const std::vector<double>& theta, double factor){
	for(size_t i=0; i < theta_.size(); i++){
		theta_[i] += factor * theta[i];
	}
}
void group::add_to_w(const std::vector<double>& w, double factor){
	for(size_t i=0; i < w_.size(); i++){
		w_[i] += factor * w[i];
	}
}
void group::add_to_r(const std::vector<topology::Vector2d>& r, double factor){
	for(size_t i=0; i < r_.size(); i++){
		r_[i] += factor * r[i];
	}
}

void group::add_to_r(const std::vector<double>& r, double factor){
	for(size_t i=0; i < r_.size(); i++){
		r_[i] += topology::Vector2d(factor * r[i],factor * r[N_ + i]);
	}
}

void group::add_to_p(const std::vector<topology::Vector2d>& p, double factor){
	for(size_t i=0; i < p_.size(); i++){
		p_[i] += factor * p[i];
	}
}
void group::add_to_p(const std::vector<double>& p, double factor){
	for(size_t i=0; i < p_.size(); i++){
		p_[i] += topology::Vector2d(factor * p[i],factor * p[N_ + i]);
	}
}

void group::add_to_coord(const std::vector<double>& coord, double factor){
	if ( group_type_ == "mxy" ) {
		for (int i=0; i < N_; i++){
			theta_[i] += factor * coord[i];
			r_[i] += (topology::Vector2d(factor * coord[N_ + i], factor * coord[2* N_ + i]));
		}
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ) {
		for (int i=0; i < N_; i++){
			theta_[i] += factor * coord[i];
		}
	} else {
		std::cerr << "! WARNING: function group::add_to_coord not defined for group type '" << group_type_ << "'. Doing nothing !" << std::endl;
	}
}
void group::add_to_coord_inertialscaling(const std::vector<double>& coord, double factor){
	if (group_type_ == "mxy") {
		double Iinv = 1.0 / I_;
		double minv = 1.0 / m_;
		for (int i=0; i < N_; i++){
			theta_[i] += factor * coord[i] * Iinv;
			r_[i] += (topology::Vector2d(factor * coord[N_ + i] * minv, factor * coord[2* N_ + i] * minv));
		}
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ) {
		double Iinv = 1.0 / I_;
		for (int i=0; i < N_; i++){
			theta_[i] += factor * coord[i] * Iinv;
		}
	} else {
		std::cerr << "! WARNING: function group::add_to_coord_inertialscaling not defined for group type '" << group_type_ << "'. Doing nothing !" << std::endl;
	}
}
void group::add_to_mom(const std::vector<double>& mom, double factor){
	if (group_type_ == "mxy") {
		for (int i=0; i < N_; i++){
			w_[i] += factor * mom[i];
			p_[i] += (topology::Vector2d(factor * mom[N_ + i], factor * mom[2* N_ + i]));
		}
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ) {
		for (int i=0; i < N_; i++){
			w_[i] += factor * mom[i];
		}
	} else {
		std::cerr << "! WARNING: function group::add_to_mom not defined for group type '" << group_type_ << "'. Doing nothing !" << std::endl;
	}
}


void group::add_random_angle(double angmax){
	for (int i=0; (unsigned)i < theta_.size(); i++){
		theta_[i] += random_angle(angmax);
	}
}
void group::add_random_displacement(double rmax){
	for (int i=0; (unsigned)i < r_.size(); i++){
		r_[i] += topology::random_vector(topology::Vector2d(-rmax,-rmax),topology::Vector2d(rmax,rmax));
	}
}


void group::stream_along_spin(double v) {
	for (int i=0; (unsigned)i < r_.size(); i++){
		r_[i] += topology::spin(theta_[i]) * v; // converts theta into vector and adds it to r.
	}
}
//////////////////////////////////////////////////////////////////////////////
// tools for energy and momentum checks
//////////////////////////////////////////////////////////////////////////////


double group::sum_w() const {
	double w=0;
	for (size_t i=0; i < w_.size(); i++){
		w += w_[i];
	}
	return w;
}

double group::sum_w_squared() const {
	double w=0;
	for (size_t i=0; i < w_.size(); i++){
		w += w_[i] * w_[i];
	}
	return w;
}

double group::sum_w_4() const {
	double w=0;
	for (size_t i=0; i < w_.size(); i++){
		w += w_[i] * w_[i] * w_[i] * w_[i];
	}
	return w;
}

double group::sum_theta() const {
	double te=0;
	for (size_t i=0; i < theta_.size(); i++){
		te += theta_[i];
	}
	return te;
}


topology::Vector2d group::sum_s() const {
	topology::Vector2d s = 0;
	for (size_t i=0; i < theta_.size(); i++){
		s += topology::spin(theta_[i]);
	}
	return s;
}

topology::Vector2d group::sum_p() const {
	topology::Vector2d p = 0;
	for (size_t i=0; i < p_.size(); i++){
		p += p_[i];
	}
	return p;
}
double group::sum_p_squared() const {
	double p_squared = 0;
	for (size_t i=0; i < p_.size(); i++){
		p_squared += topology::norm2(p_[i]);
	}
	return p_squared;
}
double group::sum_p_4() const {
	double p_4 = 0;
	for (size_t i=0; i < p_.size(); i++){
		p_4 += std::pow(topology::norm2(p_[i]),2);
	}
	return p_4;
}

double group::sum_e_squared() const {
	double e=0;
	for (int i=0; i < N_; i++){
		e += std::pow(calc_energy(i),2);
	}
	return e;
}


double group::sum_eint_squared() const {
	double eint=0;
	for (int i=0; i < N_; i++){
		eint += std::pow(calc_interaction_energy(i),2);
	}
	return eint;
}

double group::sum_ekin_squared() const {
	double ekin=0;
	for (int i=0; i < N_; i++){
		ekin += std::pow(calc_kinetic_energy(i),2);
	}
	return ekin;
}


// Attention! Uses a half lattice.
double group::calc_interaction_energy() const {
	double E = 0;
	// variable containing the neighbors of the current particle
	std::vector<int> nb;
	std::vector<double> distances;
	if ( group_type_ == "mxy" ){
		for (int i=0; i < N_; i++){
			nb = get_neighbors(i,nb_rule_,distances);
			for (std::size_t j = 0; j < nb.size(); j++ ){
				if (i != nb[j]){ // Particle at index i is contained in nb. Here to ignore wrong counting.
					E += - J_pot(distances[j]) * cos(theta_diff(i,nb[j])) + U_pot(distances[j]);
				}
			}
		}
	} else if ( group_type_ == "fmxy" ){
		for (int i=0; i < N_; i++){
			nb = get_neighbors(i,nb_rule_,distances);
			for (std::size_t j = 0; j < nb.size(); j++ ){
				if (i != nb[j]){ // Particle at index i is contained in nb. Here to ignore wrong counting.
					E += - J_pot(distances[j]) * cos(theta_diff(i,nb[j]));
				}
			}
		}
	} else if ( group_type_ == "xy" ) {
		for (int i=0; i < N_; i++){
			nb = get_neighbors(i,nb_rule_,distances);
			for (std::size_t j = 0; j < nb.size(); j++ ){
				if (i != nb[j]){ // Particle at index i is contained in nb. Here to ignore wrong counting
					E += cos(theta_diff(i,nb[j]));
				}
			}
			E *=  - J_ ;
		}
	}
	return nb_mult_factor_ * E;
}


double group::calc_interaction_energy(int index) const {
	double E = 0;
	// variable containing the neighbors of the current particle
	std::vector<int> nb;
	std::vector<double> distances;
	nb = get_neighbors(index,"all",distances);
	for (std::size_t j = 0; j < nb.size(); j++ ){
		if (index != nb[j]){ // Particle at index i is contained in nb. Here to ignore wrong counting.
			if ( group_type_ == "mxy" ){
				E += - J_pot(distances[j]) * cos(theta_diff(index,nb[j])) + U_pot(distances[j]);
			} else if ( group_type_ == "fmxy" ){
				E += - J_pot(distances[j]) * cos(theta_diff(index,nb[j]));
			} else if ( group_type_ == "xy" ){
				E += cos(theta_diff(index,nb[j]));
			}
		}
	}
	if ( group_type_ == "xy" ){
		E *= - J_;
	}
	return .5 * E;
}

double group::calc_kinetic_energy() const {
	if ( group_type_ == "mxy" ){
		return .5 * sum_w_squared()/ I_ + .5 * sum_p_squared()/ m_;
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ){
		return .5 * sum_w_squared()/ I_;
	}
	return 0; // To avoid warning message
}
double group::calc_kinetic_energy(int i) const {
	if ( group_type_ == "mxy" ){
		return .5 * std::pow(w_[i],2)/ I_ + .5 * p_[i].norm2()/ m_;
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ){
		return .5 * std::pow(w_[i],2)/ I_ ;
	}
	return 0; // To avoid warning message
}

double group::calc_temperature() const {
	if ( group_type_ == "mxy" ){
		return ( (sum_w_squared()- pow(sum_w(), 2.0) / N_)/ I_
				+ (sum_p_squared()   - topology::norm2(sum_p())/ N_ )/ m_)/(3 * N_ - 3); // 3 N_ - 3 degrees of freedom
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ){
		// return (sum_w_squared()- pow(sum_w(), 2.0) / N_)/ I_;
		return (sum_w_squared() - pow(sum_w(), 2.0)/ N_ )/(N_ * I_);
	}
	return 0;
}

std::vector<int> group::plaquette(int i) const {
	if ( group_type_ == "xy" || group_type_ == "fvm"){
		std::vector<double> distances;
		std::vector<int> plaq = get_neighbors(i, "ur", distances);
		std::vector<int> nb;
		// plaq now contains the right neighbor and the top neighbor.
		if (lattice_type_ == 't'){
			plaq.insert(plaq.begin(),i);
		}
		if (lattice_type_ == 's'){
			plaq.insert(plaq.begin(),i);
			// In a square lattice, we also need the top right neighbor. That is the top neighbor of the right neighbor.
			nb = get_neighbors(plaq[1], "ur", distances);
			plaq.insert(plaq.begin()+2,nb[1]);
		}
		// plaq goes around once, so the index i is in there twice.
		plaq.push_back(i);

		return plaq;
	} else {
		std::cerr 	<< " ! WARNING: in group::plaquette. Function ill-defined for group type '" << group_type_ << "' ! " << std::endl
					<< " !          returning empty vector." << std::endl;
		return std::vector<int>();
	}
}

double group::calc_vorticity(int index) const {
	double plaqsum = 0;
	std::vector<int> plaq = plaquette(index);
	if ( group_type_ == "xy" || group_type_ == "fvm"){
		for(size_t i = 0; i < plaq.size() - 1; i++){
			plaqsum += mod(theta_diff(plaq[i+1],plaq[i]),-M_PI,M_PI);
		}
	}
	return plaqsum / (2 * M_PI);
}

double group::calc_vortexdensity_unsigned() const {
	double rho_v = 0;
	if ( group_type_ == "xy" || group_type_ == "fvm"){
		for(int i = 0; i < N_; i++){
			rho_v += std::abs(calc_vorticity(i));
		}
	}
	return rho_v / (get_volume() * 2 * M_PI);
}

double group::calc_vortexdensity_signed() const {
	double rho_v = 0;
	if ( group_type_ == "xy" || group_type_ == "fvm"){
		for(int i = 0; i < N_; i++){
			rho_v += calc_vorticity(i);
		}
		return rho_v / (get_volume() * 2 * M_PI);
	} else {
		return 0;
	}

}


double group::calc_space_angular_mom() const {
	double L = 0;
	if ( r_.size() == (unsigned)N_ && p_.size() == (unsigned)N_ ) {
		for (int i=0; i < N_; i++){
			L += calc_space_angular_mom(i);
		}
	} else {
		std::cerr 	<< " ! WARNING: in group::calc_space_angular_mom. Size of r_ or p_ inappropriate ! " << std::endl
					<< " !          r_.size() = " << r_.size() << ", p_.size() = " << p_.size() << std::endl;
	}

	return L;
}

double group::calc_neighbor_mean(double te_pow, double r_pow, double cos_pow, double sin_pow, double J_pow, double Up_pow, double Upp_pow) const {
	if ( ( group_type_ == "xy" || group_type_ == "vm" || group_type_ == "fvm") && ( J_pow != 0 || Up_pow != 0 || Upp_pow != 0 )){
		std::cerr 	<< " ! WARNING: in group::calc_neighbor_mean. Problematic input. ! " << std::endl
					<< " !          values of J_pow, Up_pow and Upp_pow should be zero for group_type '" << group_type_ << "'." << std::endl
					<< " !          returning 0." << std::endl;
		return 0;
	}
	double returnval = 0;
	double auxval = 1;
	double theta_ij = 0;
	// variable containing the neighbors of the current particle
	std::vector<int> nb;
	std::vector<double> distances;
	for (int i=0; i < N_; i++){
		nb = get_neighbors(i,nb_rule_,distances);
		if ( group_type_ == "xy" || group_type_ == "fvm"){
			distances.clear();
			for (std::size_t j = 0; j < nb.size(); j++ ){
				distances.push_back(periodic_distance(i,nb[j]));
			}
		}
		for (std::size_t j = 0; j < nb.size(); j++ ){
			if (te_pow == 0 && r_pow == 0 && cos_pow == 0 && sin_pow == 0 && J_pow == 0 && Up_pow == 0 && Upp_pow == 0 ){ // For coordination number.
				returnval += nb.size();
			} else {
				if (i != nb[j]){ // Particle at index i is contained in nb. Here to ignore wrong counting.
					theta_ij = mod(theta_diff(i,nb[j]),-M_PI,M_PI);
					auxval = 1;
					if (te_pow != 0) {
						auxval *= std::pow(theta_ij,te_pow);
					}
					if (r_pow != 0) {
						auxval *= std::pow(distances[j],r_pow);
					}
					if (cos_pow != 0) {
						auxval *= std::pow(cos(theta_ij),cos_pow);
					}
					if (sin_pow != 0) {
						auxval *= std::pow(sin(theta_ij),sin_pow);
					}
					if (J_pow != 0) {
						auxval *= std::pow(J_pot(distances[j]),J_pow);
					}
					if (Up_pow != 0) {
						auxval *= std::pow(- J_pot_prime(distances[j]) * cos(theta_ij) + U_pot_prime(distances[j]) ,Up_pow);
					}
					if (Upp_pow != 0) {
						auxval *= std::pow(- J_pot_primeprime(distances[j]) * cos(theta_ij) + U_pot_primeprime(distances[j]) ,Upp_pow);
					}
					returnval += auxval;
				}
			}
		}

	}
	return nb_mult_factor_ * returnval;
}



std::vector<double> group::calc_helicity(double beta) const {
	double H_x = 0;
	double H_y = 0;
	double I_x = 0;
	double I_y = 0;

	double theta_ij;
	double J_ij;

	std::vector<int> nb;
	std::vector<double> distances;
	std::vector<topology::Vector2d> rij_vec;
	for (int i=0; i < N_; i++){
		nb = get_neighbors(i,nb_rule_,distances);
		distances.clear();
		rij_vec.clear();
		for (std::size_t j = 0; j < nb.size(); j++ ) {
			if ( i != nb[j] ) {
				distances.push_back(periodic_distance(i,nb[j]));
				rij_vec.push_back(periodic_distance_vector(i,nb[j]));
				theta_ij = mod(theta_diff(i,nb[j]),-M_PI,M_PI);
				if ( group_type_ == "mxy" || group_type_ == "fmxy" ) {
					J_ij = J_pot(distances[j]);
				} else {
					J_ij = J_;
				}
				H_x += J_ij * cos(theta_ij) * rij_vec[j].get_x() * rij_vec[j].get_x();
				H_y += J_ij * cos(theta_ij) * rij_vec[j].get_y() * rij_vec[j].get_y();
				I_x += J_ij * sin(theta_ij) * rij_vec[j].get_x();
				I_y += J_ij * sin(theta_ij) * rij_vec[j].get_y();
			}
		}
	}
	H_x *= nb_mult_factor_;
	H_y *= nb_mult_factor_;
	I_x *= nb_mult_factor_;
	I_y *= nb_mult_factor_;
	double helicity = 1/(2*L_.get_x()*L_.get_y()) * (H_x * H_y - beta * ( I_x * I_x + I_y * I_y) );
	std::vector<double> returnvec = {helicity, H_x, H_y, I_x, I_y};
	return returnvec;
}
/*	double returnval = 0;
	double auxval = 1;
	double theta_ij = 0;

	// variable containing the neighbors of the current particle
	std::vector<int> nb;
	std::vector<double> distances;
	std::vector<topology::Vector2d> rij_vec;
	for (int i=0; i < N_; i++){
		nb = get_neighbors(i,nb_rule_,distances);
		distances.clear();
		rij_vec.clear();
		for (std::size_t j = 0; j < nb.size(); j++ ) {
			if ( i != nb[j] ) {
				distances.push_back(periodic_distance(i,nb[j]));
				rij_vec.push_back(periodic_distance_vec(i,nb[j]));
				theta_ij = mod(theta_diff(i,nb[j]),-M_PI,M_PI);
				H_x += J_pot(distances[j]) * cos(theta_ij) * rij_vec[j].get_x() * rij_vec[j].get_x();
				H_y += J_pot(distances[j]) * cos(theta_ij) * rij_vec[j].get_y() * rij_vec[j].get_y();
				I_x += J_pot(distances[j]) * sin(theta_ij) * rij_vec[j].get_x();
				I_y += J_pot(distances[j]) * sin(theta_ij) * rij_vec[j].get_y();
			}
		}
	}
	H_x *= nb_mult_factor_;
	H_y *= nb_mult_factor_;
	I_x *= nb_mult_factor_;
	I_y *= nb_mult_factor_;
	double helicity = 1/(2*L_.get_x()*L_.get_y()) * (H_x * H_y - beta * ( I_x * I_x + I_y * I_y) );
	std::vector<double> returnvec = {helicity, H_x, H_y, I_x, I_y};
	return returnvec;
}
*/
//
//double group::binder_cumulant() const {
//double binder_cumulant = 0;
//topology::Vector2d magnetization_4=0;
//topology::Vector2d magnetization_2=0;
//	for (int i=0; i < N_; i++){
//		magnetization_4 += topology::spin(theta_[i])*topology::spin(theta_[i])*topology::spin(theta_[i])*topology::spin(boretheta_[i]);
//	}
//return binder_cumulant;
//}



//double group::ang_MSD(xygroup& initial) const{
//	double ang_MSD=0;
//	double angle=0;
//	for(int i=0;i<N_;i++){
//		if(theta_[i] > 0){
//		 angle = acos(cos(theta_[i]));}
//		else{
//			 angle = -acos(cos(theta_[i]));
//		}
//		if(initial[i].get_theta() > 0){  //this if bracket for initial is because of later starting MSD calculations
//			 initial[i].set_theta(acos(cos(initial[i].get_theta())));}
//			else{
//				 initial[i].set_theta( -acos(cos(initial[i].get_theta())));
//			}
////		std::cout << angle << "   " << initial[i].get_theta() << std::endl;
//		ang_MSD+= (angle-initial[i].get_theta())*(angle-initial[i].get_theta());
//	}
//
//	return ang_MSD/N_;
//}
//
//
//double group::ang_MSD_nonsat(const group& initial) const{
//	double ang_MSD=0;
//	for(int i=0;i<N_;i++){
////		std::cout << angle << "   " << initial[i].get_theta() << std::endl;
//		ang_MSD+= (theta_[i]-initial[i].get_theta())*(theta_[i]-initial[i].get_theta());
//	}
//
//	return ang_MSD/N_;
//}
////////////////////////////////////////////////////////////////////////////////
//
//double group::spin_autocorrelation(const group& corr_G) const{
//	double avg_acr = 0;
//	for(int i=0;i<N_;i++){
//            avg_acr += topology::spin(corr_G[i].get_theta())*topology::spin(theta_[i]);
//	}
//	return avg_acr/N_;
//}

std::complex<double> group::calc_mxq(const topology::Vector2d q, double Mx_0) const{
	std::complex<double> mq = 0;
	for (int j=0; j < N_; j++){
		mq += (cos(theta_[j]) - Mx_0) * calc_eiqr(q, j);
	}
	return mq / (double)sqrtN_;
}

std::complex<double> group::calc_myq(const topology::Vector2d q, double My_0) const{
	std::complex<double> mq = 0;
	for (int j=0; j < N_; j++){
		mq += (sin(theta_[j]) - My_0)* calc_eiqr(q, j);
	}
	return mq / (double)sqrtN_;
}


std::complex<double> group::calc_wq(const topology::Vector2d q, double W_0) const{
	std::complex<double> wq = 0;
	for (int j=0; j < N_; j++){
		wq += (w_[j] - W_0) * calc_eiqr(q, j);
	}
	return wq / (double)sqrtN_;
}

std::complex<double> group::calc_eq(const topology::Vector2d q, double E_0) const{
	std::complex<double> eq = 0;
	std::vector<double> distances;
	std::vector<int> nb;
	for (int j=0; j < N_; j++){
		nb = get_neighbors(j,nb_rule_,distances); // Take only upper right cells.;
		for (std::size_t k = 0; k < nb.size(); k++ ){
			if ( group_type_ == "mxy" ){
				eq += ( - J_pot(distances[k]) * cos(theta_diff(j,nb[k]))
						+ U_pot(distances[k]) ) * (calc_eiqr(q,j) + calc_eiqr(q,nb[k]));
			} else if ( group_type_ == "fmxy" ){
				eq += ( - J_pot(distances[k]) * cos(theta_diff(j,nb[k])) )
						* (calc_eiqr(q,j) + calc_eiqr(q,nb[k]));
			} else if ( group_type_ == "xy" ){
				eq += ( cos(theta_diff(j,nb[k])) ) * (calc_eiqr(q,j) + calc_eiqr(q,nb[k]));
			}
		}
		if ( group_type_ == "mxy" ){
			eq += (w_[j] * w_[j]/(2 * I_) + p_[j].norm2() / (2 * m_) - E_0) * calc_eiqr(q, j);
		} else if ( group_type_ == "xy" || group_type_ == "fmxy" ){
			eq += (w_[j] * w_[j]/(2 * I_) - E_0) * calc_eiqr(q, j);
		}
	}
	return nb_mult_factor_ * eq / (double)sqrtN_;
}


std::complex<double> group::calc_teq(const topology::Vector2d q, double Te_0) const {
	std::complex<double> teq = 0;
	for (int j=0; j < N_; j++){
		teq += (mod(theta_[j], - M_PI, M_PI) - Te_0) * calc_eiqr(q, j);
	}
	return teq / (double)sqrtN_;
}


std::complex<double> group::calc_rq(const topology::Vector2d q) const{
	std::complex<double> rq = 0;
	for (int j=0; j < N_; j++){
		rq += calc_eiqr(q, j);
	}
	return rq / (double)sqrtN_;
}


std::vector<std::complex<double> > group::calc_jq(const topology::Vector2d q, topology::Vector2d J_0) const {
	std::vector<std::complex<double> > jq = {0,0};
	J_0 = J_0 * m_; // Turning it into a momentum
	for (int j=0; j < N_; j++){
		for (int k=0; k < 2; k++){
			jq[k] += (p_[j][k] - J_0[k]) * calc_eiqr(q, j);
		}
	}
	return vector_scale(jq, 1.0/ m_ / sqrtN_);
}

std::complex<double> group::calc_jqpar(const topology::Vector2d q, topology::Vector2d J_0) const{
	std::complex<double> jqpar = 0;
	J_0 = J_0 * m_; // Turning it into a momentum
	if ( q.norm2() > 0 ){
		for (int j=0; j < N_; j++){
			jqpar += topology::parallel_projection(p_[j] - J_0,q) * calc_eiqr(q, j);
		}
	} else {
		jqpar = .5 * (sum_p().get_x() + sum_p().get_y());
		std::cerr << "+  WARNING: in calc_jqpar\n" <<
				"+           q = 0" << std::endl;
	}
	return jqpar / m_ / (double)sqrtN_;
}

std::complex<double> group::calc_jqperp(const topology::Vector2d q, topology::Vector2d J_0) const{
	std::complex<double> jqperp = 0;
	J_0 = J_0 * m_; // Turning it into a momentum
	if ( q.norm2() > 0 ){
		for (int j=0; j < N_; j++){
			jqperp += topology::orthogonal_projection(p_[j] - J_0,q) * calc_eiqr(q, j);
		}
	} else {
		jqperp = 0;
		std::cerr << "+  WARNING: in calc_jqperp\n" <<
				"+           q = 0" << std::endl;
	}
	return jqperp / m_ / (double)sqrtN_;
}


std::complex<double> group::calc_lq(const topology::Vector2d q, double L_0) const{
	std::complex<double> lq = 0;
	for (int j=0; j < N_; j++){
		lq += (calc_space_angular_mom(j)  - L_0) * calc_eiqr(q, j);
	}
	return lq / (double)sqrtN_;
}


double group::calc_fieldfluct_average(std::string fluctname, topology::Vector2d q) const {
	double sum_aux = 0;
	if (fluctname == "mxq"){
		return sum_s().get_x();
	} else if (fluctname == "myq"){
		return sum_s().get_y();
	} else if (fluctname == "wq"){
		return sum_w();
	} else if (fluctname == "eq"){
		return calc_energy();
	} else if (fluctname == "teq"){
		for (int i=0; i < N_; i++){
			sum_aux += theta_[i];
		}
		return sum_aux;
	} else if (fluctname == "rq"){
		return N_;
	} else if (fluctname == "lq"){
		for (int i=0; i < N_; i++){
			sum_aux += calc_space_angular_mom(i);
		}
		return sum_aux;
	} else if (fluctname == "jparq"){
		if ( q.get_x() != 0 || q.get_y() != 0 ){
			return topology::parallel_projection(sum_p(),q);
		}
	} else if (fluctname == "jperpq"){
		if ( q.get_x() != 0 || q.get_y() != 0 ){
			return topology::orthogonal_projection(sum_p(),q);
		}
	}
	std::cerr << "ERROR: Wrong input in group::calc_one_particle_density.\n"
			<<   "       Either wrong q value or wrong fluctname.\n"
			<<   "       fluctname: " << fluctname << "\n"
			<<   "       q = (" << q.get_x() << ", " << q.get_y() << ")" << std::endl;
	exit(-1);

	return 0; // Just to keep the compiler from complaining. Error handling should be improved.
}


double group::calc_one_particle_density(int index, std::string fluctname, topology::Vector2d q) const {
	if (fluctname == "mxq"){
		return std::cos(theta_[index]);
	} else if (fluctname == "myq"){
		return std::sin(theta_[index]);
	} else if (fluctname == "wq"){
		return w_[index];
	} else if (fluctname == "eq"){
		return calc_energy(index);
	} else if (fluctname == "teq"){
		return theta_[index];
	} else if (fluctname == "rq"){
		return 1;
	} else if (fluctname == "lq"){
		return calc_space_angular_mom(index);
	} else if (fluctname == "jparq"){
		if ( q.get_x() != 0 || q.get_y() != 0 ){
			return topology::parallel_projection(p_[index],q);
		}
	} else if (fluctname == "jperpq"){
		if ( q.get_x() != 0 || q.get_y() != 0 ){
			return topology::orthogonal_projection(p_[index],q);
		}
	}
	std::cerr << "ERROR: Wrong input in group::calc_one_particle_density.\n"
			<<   "       Either wrong q value or wrong fluctname.\n"
			<<   "       fluctname: " << fluctname << "\n"
			<<   "       q = (" << q.get_x() << ", " << q.get_y() << ")" << std::endl;
	exit(-1);

	return 0; // Just to keep the compiler from complaining. Error handling should be improved.
}

std::vector<std::complex<double> > group::calc_fieldfluct(const std::vector<topology::Vector2d> qvals, std::string fluctname) const{
	std::vector<std::complex<double> > returnvec;
	if (fluctname == "mxq"){
		double Mx_0 = sum_s().get_x() / N_;
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_mxq(qvals[k], Mx_0));
		}
	} else if (fluctname == "myq"){
		double My_0 = sum_s().get_y() / N_;
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_myq(qvals[k], My_0));
		}
	} else if (fluctname == "wq"){
		double W_0 = sum_w() / N_;
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_wq(qvals[k], W_0));
		}
	} else if (fluctname == "eq"){
		double E_0 = calc_energy() / N_;
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_eq(qvals[k], E_0));
		}
	} else if (fluctname == "teq"){
		double Te_0 = sum_theta() / N_;
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_teq(qvals[k], Te_0));
		}
	} else if (fluctname == "rq"){
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_rq(qvals[k]));
		}
	} else if (fluctname == "lq"){
		double L_0 = calc_space_angular_mom() / N_;
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_lq(qvals[k], L_0));
		}
	} else if (fluctname == "jparq"){
		topology::Vector2d J_0 = sum_p() / (N_ * m_);
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_jqpar(qvals[k], J_0));
		}
	} else if (fluctname == "jperpq"){
		topology::Vector2d J_0 = sum_p() / (N_ * m_);
		for(size_t k = 0; k < qvals.size(); k++){
			returnvec.push_back(calc_jqpar(qvals[k], J_0));
		}
	} else {
		std::cerr << "WARNING: Fluctname \"" << fluctname << "\" unknown in call to group::calc_fieldfluct. Returning empty vector." << std::endl;
	}
	return returnvec;
}

std::vector<std::complex<double> > group::calc_fieldfluct_convolution(const std::vector<topology::Vector2d> qvals,
		std::string fluctname_1, std::string fluctname_2) const {
	std::vector<std::complex<double> > returnvec;
	if ( fluctname_1 == "rq" ) {
		return calc_fieldfluct(qvals, fluctname_2);
	} else if ( fluctname_2 == "rq" ) {
		return calc_fieldfluct(qvals, fluctname_1);
	}
	std::complex<double> fluct_aux;
	for(size_t k = 0; k < qvals.size(); k++){
//		double average_1 = calc_fieldfluct_average(fluctname_1,qvals[k]);
//		double average_2 = calc_fieldfluct_average(fluctname_2,qvals[k]);
		fluct_aux = 0;
		for (int j=0; j < N_; j++){
			fluct_aux += calc_one_particle_density(j,fluctname_1,qvals[k]) * // - average_1)
					calc_one_particle_density(j,fluctname_2,qvals[k]) * //- average_2)*
					calc_eiqr(qvals[k], j);
		}
		returnvec.push_back(fluct_aux);
	}
	return returnvec;
}


topology::Vector2d group::calc_tau() const {
	topology::Vector2d tau = 0;
//	std::vector<double> wderi=time_derivative_w();
//	for (int j=0; j < N_; j++){
//		tau += wderi[j] *  get_r(j);
//	}
	return tau ;
}


topology::Vector2d group::calc_je() const {
	topology::Vector2d je = 0;
//	std::vector<int> cur_nei;
//	std::vector<double> distances;
//	std::vector<int> nb;
//	for (int j=0; j < N_; j++){
//		nb = get_neighbors(j,r_,"ur",distances); // Take only upper right cells.;
//		for (std::size_t k = 0; k < nb.size(); k++ ){
//			je += J_pot(distances[k]) * (get_w(j) + get_w(nb[k]) )
//				* std::sin(theta_diff(nb[k],j)) *
//				(get_r(j) + get_r(nb[k]));
//		}
//	}
	return je * .5 / I_ ;
}



topology::Vector2d group::calc_current(std::string currentname) const {
	if (currentname == "tau") {
		return calc_tau();
	} else { // if (currentname == "je") {
		return calc_je();
	}
}



//////////////////////////////////////////////////////////////////////////////
// SCF functions
//////////////////////////////////////////////////////////////////////////////
// SCF Spin
std::vector<double> group::calc_SCF_S_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_S(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_S[binindex]+= std::cos(theta_diff(index, i));
			counts[binindex]++;
		}
	}
	return SCF_S;
}

std::vector<double> group::calc_SCF_S_oriented_individual( const int index, const std::vector<double> rbin, const double& orientation_angle, std::vector<int>& counts) const {
	std::vector<double> SCF_S(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_S[binindex]+= std::cos(theta_[i] - orientation_angle);
			counts[binindex]++;
		}
	}
	double s_local = std::cos(theta_[index] - orientation_angle);
	for (int i=0; i < (int)rbin.size(); i++){
		SCF_S[i] *= s_local;
	}
	return SCF_S;
}
/*
std::vector<double> group::calc_SCF_S( const std::vector<double> rbin, const double av_M_squared) const {
	std::vector<double> SCF_S(rbin.size(), 0.0);
	for (int i = 0; i < N_; i++){
		SCF_S = vector_sum(SCF_S,calc_SCF_S_individual(i,rbin,av_M_squared));
	}
	return vector_scale(SCF_S,1./N_);
}


std::vector<double> group::calc_SCF_S( const std::vector<double> rbin, const double av_M_squared, int number_of_points) const {
	if ( number_of_points > .5* N_){
		return	calc_SCF_S(rbin,av_M_squared);
	} else {
		std::vector<double> SCF_S(rbin.size(), 0.0);
		std::vector<int> already_included;
		int index;
		while ((int)already_included.size() < number_of_points){ // (int) just to avoid warnings
			index = random_int(0,N_ - 1);
			if ( std::find(already_included.begin(), already_included.end(), index) == already_included.end()){
				already_included.push_back(index);
				SCF_S = vector_sum(SCF_S,calc_SCF_S_individual(index,rbin,av_M_squared));
			}
		}
		return vector_scale(SCF_S,1./already_included.size());
	}
}
*/

// Direct correlation function
std::vector<double> group::calc_SCF_g_individual( const int index, const std::vector<double> rbin) const {
	std::vector<double> g(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if (dist_squared < maxdist_squared ){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			g[binindex]++;
		}
	}
	g[0]--; // ignores the particle itself, which appears at index 0.
	g[0] /= (2 * M_PI * rbin[0] * rbin[0] * get_density());
	for (size_t k = 1; k < rbin.size(); k++){
		if (g[k] > 0){ // just a small speedup
			g[k] /= (2 * M_PI * rbin[k] * (rbin[k] - rbin[k - 1] ) * get_density());
		}
	}
	return g;
}


std::vector<double> group::calc_SCF_g( const std::vector<double> rbin, int number_of_points) const {
	std::vector<double> g(rbin.size(), 0.0);
	int index;

	std::vector<int> particle_selection(N_);
	std::iota(particle_selection.begin(),particle_selection.end(),0);
	std::random_shuffle(particle_selection.begin(),particle_selection.end());


	for (int i = 0; i < std::min(number_of_points,N_);i++) {
		index = particle_selection[i];
		g = vector_sum(g,calc_SCF_g_individual(index,rbin));
	}
	return vector_scale(g,1./std::min(number_of_points,N_));
}


// Anglediff
std::vector<double> group::calc_SCF_anglediff_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_anglediff(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_anglediff[binindex]+= theta_diff(index, i);
			counts[binindex]++;
		}
	}
	return SCF_anglediff;
}


std::vector<double> group::calc_SCF_E_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_E(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_E[binindex]+= calc_energy(i);
			counts[binindex]++;
		}
	}
	double E_self = calc_energy(index);
	for (size_t k = 0; k < rbin.size(); k++){
		if (counts[k] > 0){
			SCF_E[k] *= E_self ;
		}
	}
	return SCF_E;
}

std::vector<double> group::calc_SCF_Ekin_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_Ekin(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_Ekin[binindex]+= calc_kinetic_energy(i);
			counts[binindex]++;
		}
	}
	double E_self = calc_kinetic_energy(index);
	for (size_t k = 0; k < rbin.size(); k++){
		if (counts[k] > 0){
			SCF_Ekin[k] *= E_self;
		}
	}
	return SCF_Ekin;
}

std::vector<double> group::calc_SCF_Eint_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_Eint(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_Eint[binindex]+= calc_interaction_energy(i);
			counts[binindex]++;
		}
	}
	double E_self = calc_interaction_energy(index);
	for (size_t k = 0; k < rbin.size(); k++){
		if (counts[k] > 0){
			SCF_Eint[k] *= E_self;
		}
	}
	return SCF_Eint;
}

// Momentum
std::vector<double> group::calc_SCF_P_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_P(rbin.size(), 0.0);
	topology::Vector2d p_self = p_[index];
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_P[binindex]+= topology::innerproduct(p_self,p_[i]);
			counts[binindex]++;
		}
	}
	return SCF_P;
}

std::vector<double> group::calc_SCF_W_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const {
	std::vector<double> SCF_w(rbin.size(), 0.0);
	double dist_squared;
	int binindex;
	double maxdist_squared = std::pow(rbin.back(),2);
	for (int i = 0; i < N_; i++){
		dist_squared = periodic_distance_squared(index,i);
		if ((index != i) && (dist_squared < maxdist_squared )){
			binindex = find_bin(rbin, std::sqrt(dist_squared));
			SCF_w[binindex]+= w_[i];
			counts[binindex]++;
		}
	}
	double w_self = w_[index];
	for (size_t k = 0; k < rbin.size(); k++){
		if (counts[k] > 0){
			SCF_w[k] *= w_self;
		}
	}
	return SCF_w;
}


std::vector<double> group::calc_SCF_averaged( const std::vector<double> rbin,
		int number_of_points, std::string name) const{
	std::vector<double> return_vec(rbin.size(), 0.0);
	std::vector<int> counts(rbin.size(), 0);

	std::vector<int> particle_selection(N_);
	std::iota(particle_selection.begin(),particle_selection.end(),0);
	std::random_shuffle(particle_selection.begin(),particle_selection.end());


	double orientation_par, orientation_perp;
	int index;
	if (name == "S_par"){
		topology::Vector2d m_cur = sum_s() / N_;
		orientation_par = topology::angle_from_vector(m_cur);
	} else if (name == "S_perp"){
		topology::Vector2d m_cur = sum_s() / N_;
		orientation_perp = topology::angle_from_vector(m_cur) + .5 * M_PI;
	}
	if ( (name != "g") && (name != "anglediff") && (name != "S") && (name != "S_par") && (name != "S_perp")
				&& (name != "P") && (name != "W") && (name != "E") && (name != "Ekin") && (name != "Eint")){
		std::cerr << "WARNING: Name \"" << name << "\" unknown in call to group::calc_SCF_averaged. Returning 0 vector." << std::endl;
		return std::vector<double>(rbin.size());
	}
	for (int i = 0; i < std::min(number_of_points,N_);i++) {
		index = particle_selection[i];
		if ( name == "g" ){
			return_vec = vector_sum(return_vec,calc_SCF_g_individual(index,rbin));
		} else if (name == "anglediff"){
			return_vec = vector_sum(return_vec,calc_SCF_anglediff_individual(index, rbin, counts));
		} else if (name == "S"){
			return_vec = vector_sum(return_vec,calc_SCF_S_individual(index, rbin, counts));
		} else if (name == "S_par"){
			return_vec = vector_sum(return_vec,calc_SCF_S_oriented_individual(index, rbin, orientation_par, counts));
		} else if (name == "S_perp"){
			return_vec = vector_sum(return_vec,calc_SCF_S_oriented_individual(index, rbin, orientation_perp, counts));
		} else if (name == "P"){
			return_vec = vector_sum(return_vec,calc_SCF_P_individual(index, rbin, counts));
		} else if (name == "W"){
			return_vec = vector_sum(return_vec,calc_SCF_W_individual(index, rbin, counts));
		} else if (name == "E"){
			return_vec = vector_sum(return_vec,calc_SCF_E_individual(index, rbin, counts));
		} else if (name == "Ekin"){
			return_vec = vector_sum(return_vec,calc_SCF_Ekin_individual(index, rbin, counts));
		} else if (name == "Eint"){
			return_vec = vector_sum(return_vec,calc_SCF_Eint_individual(index, rbin, counts));
		}
	}
	return vector_scale(return_vec,1./std::min(number_of_points,N_));
}







//
//
//double group::calc_vorticity(int index) const {
//	std::vector<int> plaq = plaquette(index);
//	double plaqsum = 0;
//	for(size_t i = 0; i < plaq.size() - 1; i++){
//		plaqsum += mod(theta_diff(plaq[i+1],plaq[i]),-M_PI,M_PI);
//	}
//	return plaqsum;
//}
//
//double group::calc_vortexdensity_unsigned() const {
//	double rho_v = 0;
//	for(int i = 0; i < N_; i++){
//		rho_v += std::abs(calc_vorticity(i));
//	}
//	return rho_v / (get_box_size() * get_box_size() * 2 * M_PI);
//}
//
//double group::calc_vortexdensity_signed() const {
//	double rho_v = 0;
//	for(int i = 0; i < N_; i++){
//		rho_v += calc_vorticity(i);
//	}
//	return rho_v / (get_box_size() * get_box_size() * 2 * M_PI);
//}

//void group::calc_field_fluct(const topology::Vector2d q, std::complex<double>& mxq,
//		std::complex<double>& wq, std::complex<double>& eq, std::complex<double>& teq) const {
//	mxq = calc_mxq(q);
//	wq = calc_wq(q);
//	eq = calc_eq(q);
//	teq = calc_teq(q);
//}
//
//void group::calc_field_fluct(const std::vector<topology::Vector2d > q, std::vector<std::complex<double> >& mxq,
//			std::vector<std::complex<double> >& wq, std::vector<std::complex<double> >& eq,
//			std::vector<std::complex<double> >& teq) const {
//	mxq.clear();
//	wq.clear();
//	eq.clear();
//	teq.clear();
//	for (size_t i=0; i < q.size(); i++){
//		mxq.push_back(calc_mxq(q[i]));
//		wq.push_back(calc_wq(q[i]));
//		eq.push_back(calc_eq(q[i]));
//		teq.push_back(calc_teq(q[i]));
//	}
//}
//
//void group::calc_susceptibilities(const topology::Vector2d q, std::complex<double>& chimxq,
//		std::complex<double>& chiwq, std::complex<double>& chieq, std::complex<double>& chiteq) const {
//	chimxq = calc_chimxq(q);
//	chiwq = calc_chiwq(q);
//	chieq = calc_chieq(q);
//	chiteq = calc_chiteq(q);
//}
//
//void group::calc_susceptibilities(const std::vector<topology::Vector2d > q, std::vector<std::complex<double> >& chimxq,
//		std::vector<std::complex<double> >& chiwq, std::vector<std::complex<double> >& chieq,
//		std::vector<std::complex<double> >& chiteq) const {
//	chimxq.clear();
//	chiwq.clear();
//	chieq.clear();
//	chiteq.clear();
//	for (size_t i=0; i < q.size(); i++){
//		chimxq.push_back(calc_chimxq(q[i]));
//		chiwq.push_back(calc_chiwq(q[i]));
//		chieq.push_back(calc_chieq(q[i]));
//		chiteq.push_back(calc_chiteq(q[i]));
//	}
//}



double group::calc_ACF_S(const group& G_initial) const{
	double sum_ACF_S = 0;
	for (int i = 0; i < N_; i++){
		sum_ACF_S += std::cos(theta_[i] - G_initial.get_theta(i));
	}
	return sum_ACF_S / (double)N_ ;
	// Does not subtract total sums -- averages can be subtracted manually in data analysis.
	// Since they are no conserved quantities, subtraction here could lead to problems
	// because of differences to the actual average, I think.
//	return (sum_ACF_S / (double)N_ - topology::innerproduct(G_initial.sum_s(),sum_s())/ ((double)N_ * (double)N_));
}

double group::calc_ACF_anglediff(const group& G_initial) const{
	double sum_ACF_anglediff = 0;
	for (int i = 0; i < N_; i++){
		sum_ACF_anglediff += std::pow(theta_[i] - G_initial.get_theta(i),2);
	}
	return sum_ACF_anglediff / (double)N_;
	// Does not subtract total sums -- averages can be subtracted manually in data analysis.
	// Since they are no conserved quantities, subtraction here could lead to problems
	// because of differences to the actual average, I think.
//	return (sum_ACF_anglediff / (double)N_ - G_initial.sum_theta() * sum_theta()/((double)N_ * (double)N_));
}

double group::calc_ACF_sp(const group& G_initial, const std::string name) const{
	double sum_ACF = 0;
	topology::Vector2d M_old = 0, M_new = 0; // average magnetization, useful for some of the quantities. Avoids repeated calculation.
	double phi_old = 0, phi_new = 0;         // average magnetization angle, useful for some of the quantities. Avoids repeated calculation.
	if ( (name == "Ppar") || (name == "Pperp") || (name == "Spar") || (name == "Sperp")) {
		M_old = G_initial.sum_s();
		M_new = sum_s();
		if (topology::norm2(M_old) < 1e-12 || topology::norm2(M_new) < 1e-12){
			return 0; // to avoid numerical errors. Error bar maybe even a bit too low, but we'll see.
		} else {
			M_old.normalized();
			M_new.normalized();
			phi_old = M_old.angle();
			phi_new = M_new.angle();
		}
	}
	for (int i = 0; i < N_; i++){
		if (name == "S"){
			sum_ACF += std::cos(theta_[i] - G_initial.get_theta(i));
		} else if (name == "Sx"){
			sum_ACF += std::cos(theta_[i]) * std::cos(G_initial.get_theta(i));
		} else if (name == "Sy"){
			sum_ACF += std::sin(theta_[i]) * std::sin(G_initial.get_theta(i));
		} else if (name == "Spar"){
			sum_ACF += std::cos(theta_[i] - phi_new) * std::cos(G_initial.get_theta(i) - phi_old);
		} else if (name == "Sperp"){
			sum_ACF += std::sin(theta_[i] - phi_new) * std::sin(G_initial.get_theta(i) - phi_old);
		} else if (name == "anglediff"){
			sum_ACF += std::pow(theta_[i] - G_initial.get_theta(i),2);
		} else if (name == "P"){
			sum_ACF += topology::innerproduct(p_[i], G_initial.get_p(i));
		} else if (name == "Px"){
			sum_ACF += p_[i].get_x() * G_initial.get_p(i).get_x();
		} else if (name == "Py"){
			sum_ACF += p_[i].get_y() * G_initial.get_p(i).get_y();
		} else if (name == "Ppar"){
			sum_ACF += topology::parallel_projection(p_[i],M_new)
					* topology::parallel_projection(G_initial.get_p(i),M_old);
		} else if (name == "Pperp"){
			sum_ACF += topology::orthogonal_projection(p_[i],M_new)
					* topology::orthogonal_projection(G_initial.get_p(i),M_old);
		} else if (name == "W"){
			sum_ACF += w_[i] * G_initial.get_w(i);
		} else if (name == "E"){
			sum_ACF += calc_energy(i) * G_initial.calc_energy(i);
		} else if (name == "Eint"){
			sum_ACF += calc_interaction_energy(i) * G_initial.calc_interaction_energy(i);
		} else if (name == "Ekin"){
			sum_ACF += calc_kinetic_energy(i) * G_initial.calc_kinetic_energy(i);
		} else if (name == "MSD"){ // does not work with boundary condition updates.
			sum_ACF += topology::periodic_distance_squared(r_[i], G_initial.get_r(i),get_L());
		} else {
			std::cerr << "WARNING: Name \"" << name << "\" unknown in call to group::calc_ACF_sp. Returning 0." << std::endl;
			return 0;
		}
	}
	return sum_ACF / (double)N_;
}

double group::calc_ACF_q0(const group& G_initial, const std::string name) const{
	double quantity_old = 0;		// the old quantity
	topology::Vector2d aux_old; 	// sometimes necessary
	double quantity_new = 0;		// the new quantity
	topology::Vector2d aux_new; 	// sometimes necessary
	if (name == "M" || name == "S"){
		aux_old = G_initial.sum_s() / (double)N_;
		aux_new = sum_s() / (double)N_;
		return topology::innerproduct(aux_old,aux_new);
	} else if (name == "absM" || name == "mpar"){
		quantity_old = std::sqrt(G_initial.sum_s().norm2()) / (double)N_;
		quantity_new = std::sqrt(sum_s().norm2()) / (double)N_;
	} else if (name == "Ppar"){
		aux_old = G_initial.sum_s();
		quantity_old = topology::parallel_projection(sum_p(),aux_old) / (double)N_;
		quantity_new = topology::parallel_projection(sum_p(),aux_old) / (double)N_;
	} else if (name == "Pperp"){
		aux_old = G_initial.sum_s();
		quantity_old = topology::orthogonal_projection(sum_p(),aux_old) / (double)N_;
		quantity_new = topology::orthogonal_projection(sum_p(),aux_old) / (double)N_;
	} else {
		std::cerr << "WARNING: Name \"" << name << "\" unknown in call to group::calc_ACF_sp. Returning 0." << std::endl;
		return 0;
	}
	return quantity_new * quantity_old;
}
std::vector<std::complex<double> > group::calc_TCF(const group& G_initial, std::vector<topology::Vector2d> qvals,
		std::string fluctname_initial, std::string fluctname_current) const {
	std::vector<std::complex<double> > fieldfluct_initial = G_initial.calc_fieldfluct(qvals,fluctname_initial);
	std::vector<std::complex<double> > fieldfluct_current = calc_fieldfluct(qvals,fluctname_current);
	for (size_t i = 0; i < qvals.size(); i++){
		fieldfluct_current[i] *= std::conj(fieldfluct_initial[i]);
	}
	return fieldfluct_current;
}
/*====================================================================
double group::ss_correlation(const group& corr_G) const{
	double somedouble = 0;

	std::vector <std::vector<double> > dist_spin(2, std::vector<double>(N_, 1.23));


	for(int i=0;i<N_;i++){
		for(int j=0;j<N_;j++){
dist_spin[0][j]=sqrt(topology::norm2(r_[i]-particles_[j].get_r()));
dist_spin[1][j]=topology::spin(theta_[i])*topology::spin(theta_[j]);

		}
	}
//	std::vector <std::vector<double> > vec2D(4, std::vector<double>(4, 1.23));  //matrix example
//	vec2D[1][1]=0.346;
//	std::vector<double> size(vec2D.size(),vec2D[0].size());
//	 std::cout << "vec2Dsize = " << size[0] <<"x" << size[1] << std::endl;
//
//	for(auto vec : vec2D){
//		for(auto x : vec)
//			std::cout<<x << " ";
//		std::cout << std::endl;
//	}


//still some work in progress here
	return somedouble/N_;
}


=================================================================   */


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// functions adding and removing particles

/* void group::add_particle(mxyparticle part) {
   particles_.push_back(part);
   N_++;
}


void group::remove_particle(int index) {
   particles_.erase(particles_.begin() + index);
   N_--;
}
*/

std::vector<double> group::time_derivative_theta() const {
	if ( group_type_ == "mxy" || group_type_ == "xy" || group_type_ == "fmxy" ) {
		std::vector<double> theta(N_);
		for (int i = 0; i < N_; i++){
			// Time derivative of theta
			theta[i] = w_[i]/I_;
		}
		return theta;
	} else {
		std::cerr 	<< "! WARNING Function group::time_derivative_theta not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning empty vector." << std::endl;
		return std::vector<double>();
	}
}



std::vector<double> group::time_derivative_w() const {
	if ( group_type_ == "mxy" || group_type_ == "xy" || group_type_ == "fmxy" ) {
		std::vector<double> w(N_,0.0), distances;
		std::vector<int> nb;
		double w_ia = 0;

		for (int i = 0; i < N_; i++) {
			nb = get_neighbors(i,nb_rule_,distances);
			for (std::size_t j = 0; j < nb.size(); j++ ){
				if ( i != nb[j] ){ // Particle at index i is contained in nb as well. Here to ignore wrong counting.
					if ( group_type_ == "mxy" || group_type_ == "fmxy" ) {
						w_ia = - J_pot(distances[j]) * sin(theta_diff(nb[j],i));	// Time derivative of w. Interaction sum.
					} else if ( group_type_ == "xy" ) {
						w_ia = - J_ * sin(theta_diff(nb[j],i));
					}
					w[i] += w_ia;
					w[nb[j]] -= w_ia;
				}
			}
		}
		if (nb_rule_ == "ur") {
			return w;
		} else {
			return vector_scale(w,nb_mult_factor_);
		}
	} else {
		std::cerr 	<< "! WARNING Function group::time_derivative_w not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning empty vector." << std::endl;
	}
	return std::vector<double>();
}


std::vector<double> group::time_derivative_r() const {
	if ( group_type_ == "mxy" ) {
		std::vector<double> r(2*N_);
		for (int i = 0; i < N_; i++){
			// Time derivative of r
			r[i] = p_[i].get_x()/m_;
			r[N_+i] = p_[i].get_y()/m_;
		}
		return r;
	} else {
		std::cerr 	<< "! WARNING Function group::time_derivative_r not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning empty vector." << std::endl;
		return std::vector<double>();
	}
}


std::vector<double> group::time_derivative_p() const {
	if ( group_type_ == "mxy" ) {
		std::vector<double> p(2 * N_, 0.0), distances;
		topology::Vector2d unit_vec;
		topology::Vector2d p_ia;
		std::vector<int> nb;
		for (int i = 0; i < N_; i++) {
			nb = get_neighbors(i,nb_rule_,distances); // Take only upper right cells.
				// Avoids double counting, shortens computation time.
			for (std::size_t j = 0; j < nb.size(); j++ ){
				if (distances[j] > 0){ // Particle at index i is contained in nb as well. Here to ignore wrong counting.
					unit_vec = periodic_distance_vector(nb[j], i) / distances[j];
					p_ia = - U_pot_prime(distances[j]) * unit_vec; 								// Time derivative of p. Spatial part.
					p_ia += J_pot_prime(distances[j]) * cos(theta_diff(nb[j],i)) * unit_vec; 	// Time derivative of p. Spin part.
					p[i] += p_ia.get_x();
					p[N_ + i] += p_ia.get_y();
					p[nb[j]] -= p_ia.get_x();
					p[N_ + nb[j]] -= p_ia.get_y();
				}
			}
		}
		if (nb_rule_ == "ur") {
			return p;
		} else {
			return vector_scale(p,nb_mult_factor_);
		}
	} else {
		std::cerr 	<< "! WARNING Function group::time_derivative_p not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning zero vector." << std::endl;
		return std::vector<double>();
	}
}



std::vector<double> group::time_derivative_coord() const {
	if ( group_type_ == "mxy" ) {
		std::vector<double> coord_dot(3 * N_);
		double minv = 1.0/m_; // division is more expensive than multiplication. But this may be overdoing it ;)
		double Iinv = 1.0/I_;
		for (int i = 0 ; i < N_; i++){
			coord_dot[i] = w_[i] * Iinv;
			coord_dot[N_ + i] = p_[i].get_x() * minv;
			coord_dot[2 * N_ + i] = p_[i].get_y() * minv;
		}
		return coord_dot;
	} else if ( group_type_ == "xy" || group_type_ == "fmxy"  ) {
		return time_derivative_theta();
	} else {
		std::cerr 	<< "! WARNING Function group::time_derivative_coord not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning empty vector." << std::endl;
		return std::vector<double>();
	}
}

std::vector<double> group::time_derivative_mom() const {
	if ( group_type_ == "mxy" ) {
		std::vector<double> mom_dot(3 * N_,0.0);
	//	mom_dot.reserve(3 * N_);
		std::vector<double> distances;
	//	std::fill(mom_dot.begin(), mom_dot.end(), 0.0);
		topology::Vector2d unit_vec;
		topology::Vector2d p_ia; // interaction of p.
		double w_ia; // interaction of w.
		std::vector<int> nb;
		for (int i = 0; i < N_; i++) {
			nb = get_neighbors(i,nb_rule_,distances); // Take only upper right cells.
				// Avoids double counting, shortens computation time.
			for (std::size_t j = 0; j < nb.size(); j++ ){
				if (distances[j] > 1e-10){ // If particles are too close, computation of unit_vec may disturb the computation
					w_ia = - J_pot(distances[j]) * sin(theta_diff(nb[j],i));	// Time derivative of w. Interaction sum.
					mom_dot[i] += w_ia;
					mom_dot[nb[j]] -= w_ia;
					unit_vec = periodic_distance_vector(nb[j], i) / distances[j];
					p_ia = (- U_pot_prime(distances[j]) + J_pot_prime(distances[j]) * cos(theta_diff(nb[j],i))) * unit_vec;
					mom_dot[N_ + i] += p_ia.get_x();
					mom_dot[2*N_ + i] += p_ia.get_y();
					mom_dot[N_ + nb[j]] -= p_ia.get_x();
					mom_dot[2*N_ + nb[j]] -= p_ia.get_y();
				}
			}
		}
		if (nb_rule_ == "ur") {
			return mom_dot;
		} else {
			return vector_scale(mom_dot,nb_mult_factor_);
		}
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" ) {
		return time_derivative_w();
	} else {
		std::cerr 	<< "! WARNING Function group::time_derivative_mom not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning empty vector." << std::endl;
		return std::vector<double>();
	}
}


group group::time_derivative() const {
	group G(N_,group_type_);
	G.initialize_zero();
	if ( ! theta_.empty() ) {
		G.add_to_theta(time_derivative_theta());
	}
	if ( ! w_.empty() ) {
		G.add_to_w(time_derivative_w());
	}
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {// Here, group type is important, r can be non-empty for static groups
		G.add_to_r(time_derivative_r());
	}
	if ( ! p_.empty() ) {
		G.add_to_p(time_derivative_p());
	}
	return G;
}




std::vector<double> group::coord_diff(const group& G) const {
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {
		std::vector<double> coord_diff(3 * N_);
		topology::Vector2d aux_distance;
		for (int i = 0 ; i < N_; i++){
			aux_distance = topology::periodic_distance_vector(r_[i], G.get_r(i), L_);
			coord_diff[i] = mod(theta_[i] - G.get_theta(i), - M_PI, M_PI);
			coord_diff[N_ + i] = aux_distance.get_x();
			coord_diff[2 * N_ + i] = aux_distance.get_y();
		}
		return coord_diff;
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" || group_type_ == "fvm" ) {
		std::vector<double> coord_diff(N_);
		for (int i = 0 ; i < N_; i++){
			coord_diff[i] = mod(theta_[i] - G.get_theta(i), - M_PI, M_PI);
		}
		return coord_diff;
	} else {
		std::cerr 	<< "! WARNING Function group::coord_diff not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Returning empty vector." << std::endl;
		return std::vector<double>();
	}
}


void group::accumulative_MSD(std::vector<double>& MSD, const group& last_G) const {
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {
		topology::Vector2d aux_distance;
		for (int i = 0 ; i < N_; i++){
			aux_distance = topology::periodic_distance_vector(r_[i], last_G.get_r(i), L_);
			MSD[i] += mod(theta_[i] - last_G.get_theta(i), - M_PI, M_PI);
			MSD[N_ + i] += aux_distance.get_x();
			MSD[2 * N_ + i] += aux_distance.get_y();
		}
	} else if ( group_type_ == "xy" || group_type_ == "fmxy" || group_type_ == "fvm") {
		for (int i = 0 ; i < N_; i++){
			MSD[i] += mod(theta_[i] - last_G.get_theta(i), - M_PI, M_PI);
		}
	} else {
		std::cerr 	<< "! WARNING Function group::accumulative_MSD not defined for group_type_ = " << group_type_ << std::endl
					<< "!         Doing nothing." << std::endl;
	}
}

/*
void group::lf_integrate(double dt){
	std::vector<double> momdot = time_derivative_mom();
	double minvdt = 1.0 / m_ * dt;
	double Iinvdt = 1.0 / I_ * dt;
	for (int i = 0; i < N_; i ++){
		theta_[i] += theta_[i] * Iinvdt;
		r_[i] +=  p_[i] * minvdt;

		theta_[i] += .5 * dt* dt * momdot[i];
		r_[i] += topology::Vector2d(.5 * dt* dt * momdot[N_ + i],
				.5 * dt* dt * momdot[2 * N_ + i]);
	}
	set_r_to_pbc();
	fill_partition();
	add_to_mom(momdot,.5 * dt );
	add_to_mom(time_derivative_mom(),.5 * dt);

//	add_to_coord(time_derivative_coord(), dt);
//	add_to_coord_inertialscaling(momdot, .5 * dt* dt );
//
//	set_r_to_pbc();
//	fill_partition();
//	add_to_mom(momdot,.5 * dt );
//	add_to_mom(time_derivative_mom(),.5 * dt);
}
*/
//////////////////////////////////////////////////////////////////////////////
// Returns pointer to element.

//mxyparticle& group::operator[](int index) {
//	  return particles_[index];
//}
//////////////////////////////////////////////////////////////////////////////
// Returns a pointer to const element.

//const mxyparticle& group::operator[](int index) const {
//  return particles_[index];
//}

group& group::operator+=(const group& G){
	if ( ! theta_.empty() ) {
		for (int i=0; i < N_; i++){
			theta_[i] += G.get_theta(i);
		}
	}
	if ( ! w_.empty() ) {
		for (int i=0; i < N_; i++){
			w_[i] += G.get_w(i);
		}
	}
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {// Here, group type is important, r can be non-empty for static groups
		for (int i=0; i < N_; i++){
			r_[i] += G.get_r(i);
		}
	}
	if ( ! p_.empty() ) {
		for (int i=0; i < N_; i++){
			p_[i] += G.get_p(i);
		}
	}
	return *this;
}


group& group::operator*=(const double a){
	if ( ! theta_.empty() ) {
		for (int i=0; i < N_; i++){
			theta_[i] *= a;
		}
	}
	if ( ! w_.empty() ) {
		for (int i=0; i < N_; i++){
			w_[i] *= a;
		}
	}
	if ( group_type_ == "mxy" || group_type_ == "vm" ) {// Here, group type is important, r can be non-empty for static groups
		for (int i=0; i < N_; i++){
			r_[i] *= a;
		}
	}
	if ( ! p_.empty() ) {
		for (int i=0; i < N_; i++){
			p_[i] *= a;
		}
	}
	return *this;
}
