/*! \file group.h
 *
 *  \brief Header-File to class declaration of group.
 *  Introduces the group, the central data structure.
 *
 *  \date Created: 2020-02-29 (full rewrite)
 *  \date Last Updated: 2023-07-23
 *
 *  \class group.
 *  \brief A group of polar particles. Stores vectors with particle positions, velocities,
 *  spin orientations and spin rotation velocity, as well as further group properties.
 *
 *  Contains the main data to be manipulated in a simulation of the MXY model and the
 *  other models.
 *
 *  Functionalities
 *
 *  - Constructors
 *
 *  - Functions for clearing and initialization as well as handling the partition member variable
 *
 *  - Operations for reading and copying from other groups
 *
 *  - Printing operations
 *
 *  - Simple information extraction
 *
 *  - Simple arithmetic operations on individual particles and their properties (differences,
 *    scaling, setting to new values etc.)
 *
 *  - Calculation of physical properties (kinetic temperature, energz, momentum, helicity etc.)
 *
 *  - Calculation of field fluctuations
 *
 *  - Calculation for spatial and temporal correlation functions as well as correlations in
 *    reciprocal space
 *
 *  - Calculation of time derivatives
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: 2020-02-29 (full rewrite)
 *  \date Last Updated: 2023-07-23
 */
#ifndef GROUP_H_
#define GROUP_H_
#include "computations.h"
#include "topology.h"
#include "partition.h"
#include "neighbor_list.h"
#include "parameters.h"

#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <complex>
#include <random>


using namespace std::complex_literals; // makes 1i the complex unit.

// Forward declaration because of interdependencies between classes.
class neighbor_list;

class group {
public:
	group() {};
	//   ===================================================================================================
	//   Constructors:
	//   ===================================================================================================
	group(const parameters& par); ///< Constructor from values stored in parameters. Only sets simulation parameters, does not initialize particle data.

	group(const int N, const std::string group_type); ///< Reduced constructor, useful for time derivative group.
	//   ===================================================================================================

	//   ===================================================================================================
	//   Clearing, initialization, partition handling
	//   ===================================================================================================
	/// Clears particles and partition.
	void clear();

	/// Initializes particle data for the group based on parameters given
	void initialize(const parameters& par);
	/// Initializes the mobile group with random particle positions and fills the partition.
	void initialize_random(double kbT = 0);
	/// Sets particles to random values.
	void randomize_particles(double kbT = 0);
	/// Sets momenta to zero by shifts
	void mom_to_zero();
	/// Sets positions to square lattice
	void r_to_squarelattice();
	/// Sets positions to trigonal lattice. CAREFUL! Trigonal lattice does not fit well into square box.
	void r_to_trigonallattice();
	/// Sets positions to lattice. Decides which lattice depending on lattice_type_ member variable.
	void r_to_lattice();
	/// Sets all particles to zero
	void initialize_zero();
	/// Fills, i.e. computes the partition.
	void fill_partition();
	//   ===================================================================================================

	//   ===================================================================================================
	//   Reading and copying from other groups
	//   ===================================================================================================
	/// Reads coordinates and momenta from file snapshotname
	void read_from_snapshot(std::string snapshotname);

	/// Takes a subgroup (smaller group) and scales it up to the correct size of the group by copying.
	/*! Subgroups must be smaller by powers of 4. For proper results, this-object must be initialized.
	 */
	// TODO Works for MXY model, still problems with lattice models.
	void scale_from_subgroup(const group& G);



	/// Takes a subgroup (smaller group) and scales it up to the correct size of the group by copying.
	/// Subgroups must be smaller by powers of 4. For proper results, this-object must be initialized.
	/// Unlike the same function that takes a group as input, this function needs only a coordinate file.
	/// It is thus simpler to use.
	// TODO Works for MXY model, still problems with lattice models.
	void scale_from_subgroup(std::string snapshotname);
	///< Reads a subgroup (smaller group) from a file and scales it up to the correct size of the group by copying.
	//   ===================================================================================================

	//   ===================================================================================================
	//   Printing functions
	//   ===================================================================================================
	///  Prints the entire group to the outputfile
	void print_group(std::ofstream &outputfile) const;
	///  Prints only position coordinates of group to the outputfile
	void print_r(std::ofstream &outputfile) const;

	//   ===================================================================================================
	//   Extracting information
	//   ===================================================================================================
	/// Returns number of particles
	inline int get_N() const { return N_; };
	/// Same as get_N()
	inline int size() const { return N_; };
	/// Returns sqrt of number of particles
	inline int get_sqrtN() const { return sqrtN_; };
	/// Returns simulation box size
	inline topology::Vector2d get_L() const { return L_; };
	/// Returns smallest box length.
	inline double get_boxsize() const { return std::min(L_.get_x(),L_.get_y()); };
	/// Returns volume
	inline double get_volume() const { return L_.get_x() * L_.get_y() ; };
	/// Returns density.
	inline double get_density() const { return N_ / get_volume(); };
	/// Returns member variable I_ (spin inertia)
	inline double get_I() const { return I_; };
	/// Returns member variable J_ (spin coupling strength)
	inline double get_J() const { return J_; };
	/// Returns member variable m_ (particle mass)
	inline double get_m() const { return m_; };
	/// Returns member variable cutoff_ (interaction cutoff length)
	inline double get_cutoff() const { return cutoff_; };
	/// Returns member variable vm_v_ (Vicsek model velocity)
	inline double get_vm_v() const { return vm_v_; };
	/// Returns member variable vm_eta_ (Vicsek model noise strength)
	inline double get_vm_eta() const { return vm_eta_; };
	/// Returns member variable group_type_ (type of group)
	inline std::string get_group_type() const { return group_type_; };

	/// Returns member vector theta_ (spin angles). Length N.
	inline std::vector<double> get_theta() const { return theta_; };
	/// Returns member vector w_ (spin momenta). Length N.
	inline std::vector<double> get_w() const { return w_; };
	/// Returns member vector r_ (positions). Length N.
	inline std::vector<topology::Vector2d> get_r() const { return r_; };
	/// Returns member vector p_ (linear momenta). Length N.
	inline std::vector<topology::Vector2d> get_p() const { return p_; };

	/// Returns vector of all coordinates (angles theta_ and poitions r_, length 3N)
	std::vector<double> get_coord() const;
	/// Returns vector of all momenta (spin momenta w_ and linear momenta p_, length 3N)
	std::vector<double> get_mom() const;

	/// Returns spin angle theta_[i] of particle i
	inline double get_theta(int i) const { return theta_[i]; };
	/// Returns spin momentum w_[i] of particle i
	inline double get_w(int i) const { return w_[i]; };
	/// Returns position r_[i] of particle i
	inline topology::Vector2d get_r(int i) const { return r_[i]; };
	/// Returns linear momentum p_[i] of particle i
	inline topology::Vector2d get_p(int i) const { return p_[i]; };


	/// Returns spin interaction potential (distance-dependence)
	inline double J_pot(double dist) const {  return J_ * std::pow(1 - dist, 2); };
	/// Returns spatial interaction potential
	inline double U_pot(double dist) const {  return 4 * U_ * std::pow(1 - dist, 2); };
	/// Returns derivative of spin interaction potential (distance-dependence)
	inline double J_pot_prime(double dist) const {  return 2 * J_ * (dist - 1); };
	/// Returns derivative of spatial interaction potential
	inline double U_pot_prime(double dist) const {  return 8 * U_ * (dist - 1);  };
	/// Returns second derivative of spin interaction potential (distance-dependence)
	inline double J_pot_primeprime(double dist) const {  return 2 * J_ ; };
	/// Returns second derivative of spatial interaction potential
	inline double U_pot_primeprime(double dist) const {  return 8 * U_ ;  };


	//   ===================================================================================================
	//   Neighborhood determination
	//   ===================================================================================================
	/// Returns indices of neighbors of the particle. Selection of cells possible.
    /*! @param[in] i Particle index.
     *  @param[in] cellselect Cell selection command. "all" is for all neighboring cells,
     *  			"ur" is for the cell of the particle, the three cells above and the cell to the right
     *  			"single" is just for the cell the particle is in.
     *  @param[out] distances Stores the distances to all neighbors. Saves computation time. Only filled in case of mobile particles.
     *
     */
	std::vector<int> get_neighbors(int i, std::string cellselect, std::vector<double>& distances ) const;
	/// Returns indices of neighbors of the particle. Uses the member variable nb_rule_ to determine which cells to select.
	inline std::vector<int> get_neighbors(int i, std::vector<double>& distances ) const {
		return get_neighbors(i, nb_rule_, distances );
	};

	/// Fills the variables nb_index_, nb_first_, nb_dist_ according to the current neighborhood situation. Strongly recommended for fmxy model, recommended for xy and fvm model with small system sizes.
	void generate_neighbor_list();
	//   ===================================================================================================

	//   ===================================================================================================
	//   Calculating differences
	//   ===================================================================================================
	/// Difference in angles of two different particles. $\theta_{ij}$ in Bore paper.
	inline double theta_diff(int i, int j) const { return get_theta(j) - get_theta(i); };
	/// Returns squared distance between particle i and j considering periodic boundaries (square box). Squared function faster to calculate.
	inline double periodic_distance_squared(int i, int j) const {
		return topology::periodic_distance_squared(get_r(j), get_r(i), L_);
	};
	/// Returns distance between particle i and j considering periodic boundaries (square box). Taking sqrt takes more time than returning the squared quantity by dist_periodic_squared.
	inline double periodic_distance(int i, int j) const { return std::sqrt(periodic_distance_squared(i,j)); };
	/// Returns distance vector between particle i and j considering periodic boundaries (square box).
	inline topology::Vector2d periodic_distance_vector(int i, int j) const {
		return topology::periodic_distance_vector(get_r(i), get_r(j), L_);
	};

	//   ===================================================================================================
	//   Setting and scaling values
	//   ===================================================================================================
	void set_theta(double theta, int i) ; 		///< Gives theta_ of particle i a specified value.
	void set_w(double w, int i) ;				///< Gives w_ of particle i a specified value.
	void set_r(topology::Vector2d r, int i) ;	///< Gives r_ of particle i a specified value.
	void set_rx(double x, int i) ;				///< Gives x-component of r_ of particle i a specified value.
	void set_ry(double y, int i) ;				///< Gives y-component of r_ of particle i a specified value.
	void set_p(topology::Vector2d p, int i) ;	///< Gives p_ of particle i a specified value.
	void set_px(double px, int i) ;				///< Gives x-component of p_ of particle i a specified value.
	void set_py(double py, int i) ;				///< Gives y-component of p_ of particle i a specified value.
	void set_particle(double theta, double w, topology::Vector2d r, topology::Vector2d p,int i) ;
		///< Sets all values theta_, w_, r_, p_ of particle i to the designated values
	void set_all_w(double w);					///< Sets all w to given value (useful for setting T = 0)
	void set_all_p(topology::Vector2d p);		///< Sets all p to given value (useful for setting T = 0)
	void set_all_theta(double theta);			///< Sets all theta to given value (useful for perfect spin alignment)
	void set_temperature(double kT, int i);		///< Randomizes momenta to be in agreement with given kT of particle i
	void set_temperature(double kT);			///< Randomizes momenta to be in agreement with given kT of all particles.
	void set_temperature_p(double kT, int i);	///< Randomizes linear momenta to be in agreement with given kT of particle i
	void set_temperature_p(double kT);			///< Randomizes linear momenta to be in agreement with given kT of all particles.
	void set_temperature_w(double kT, int i);	///< Randomizes spin momenta to be in agreement with given kT of particle i
	void set_temperature_w(double kT);			///< Randomizes spin momenta to be in agreement with given kT of all particles.

	void scale_mom(double a);					///< Scales all momenta (w_, p_) by a factor a
	//   ===================================================================================================

	//   ===================================================================================================
	//   Adding to variables
	//   ===================================================================================================
	void add_to_theta(const std::vector<double>& theta, double factor = 1); ///< Adds vector of theta to theta_, scales by factor. Vector must have length of at least N_. Vector of doubles
	void add_to_r(const std::vector<topology::Vector2d>& r, double factor = 1); ///< Adds vector of r to r_, scales by factor. Vector must have length of at least 2 * N_. Vector of doubles
	void add_to_r(const std::vector<double>& r, double factor = 1); ///< Adds vector of r to r_, scales by factor. Vector must have length of at least N_. Vecot of topology::Vector2d
	void add_to_coord(const std::vector<double>& coord, double factor = 1); ///< Adds vector of coord to all coordinates (r and theta, if available). Vector must have length of at least N_. Vector of doubles. Does nothing for Vicsek type models.
	void add_to_coord_inertialscaling(const std::vector<double>& coord, double factor = 1); ///< Adds vector of coord to all coordinates (r and theta, if available). Scales by factor and inverse intertia (1/m or 1/I, respectively). Vector must have length of at least N_. Vector of doubles. Does nothing for Vicsek type models.
	void add_to_w(const std::vector<double>& w, double factor = 1); ///< Adds vector of w to w_, scales by factor. Vector must have length of at least N_. Vector of doubles
	void add_to_p(const std::vector<topology::Vector2d>& p, double factor = 1); ///< Adds vector of p to p_, scales by factor. Vector must have length of at least N_. Vector of doubles
	void add_to_p(const std::vector<double>& p, double factor = 1); ///< Adds vector of p to p_, scales by factor. Vector must have length of at least N_. Vector of topology::Vector2d
	void add_to_mom(const std::vector<double>& mom, double factor = 1); ///< Adds vector of mom to momenta (w and p, if available), scales by factor. Vector must have length of at least N_. Vector of doubles. Does nothing for Vicsek type models.
	void add_random_angle(double angmax);			///< Adds a uniformly distributed angle in (-angmax, angmax) to each particle's theta_.
	void add_random_displacement(double rmax);		///< Adds a uniformly distributed displacement (-rmax, rmax)^2 to each particle's r_.
	void stream_along_spin(double v); 				///< Streams along spin, r_new = r + v * spin(theta)

	//   ===================================================================================================
	//	Peridic boundary handling
	//   ===================================================================================================
	void set_theta_to_interval(); ///< Sets theta_ values to interval (-pi, pi)
	void set_r_to_pbc(); ///< Sets particle positions according to boundary conditions
	//   ===================================================================================================

	//   ===================================================================================================
	//	Summing over particles
	//   ===================================================================================================
	/// Returns sum over omega, basically N_\<w\>. Extensive.
	double sum_w() const;
	/// Omega squared, basically N_<w^2>. Proportional to kinetic energy. Extensive.
	double sum_w_squared() const ;
	/// Omega to the fourth power, basically N_<w^4>. Proportional to kinetic energy. Extensive.
	double sum_w_4() const ;
	/// Returns sum over all theta, basically N_<theta>. Extensive. Probably pointless.
	double sum_theta() const;
	/// Magnetization. Basically N_\<s\>. Extensive.
	topology::Vector2d sum_s() const;
	/// Magnetization squared. Basically N_\<s\>^2. Extensive.
	inline double sum_s_squared() const { return topology::norm2(sum_s()); };
	/// Magnetization to the fourth power. Basically N_\<s\>^4. Extensive.
	inline double sum_s_4() const { return std::pow(sum_s_squared(),2); };
	/// Total momentum. Basically N_\<p\>. Extensive.
	topology::Vector2d sum_p() const;
	/// Sum over momentum squared. Basically N_<p^2>. Extensive.
	double sum_p_squared() const ;
	/// Sum over momentum to the fourth power. Basically N_<p^4>. Extensive.
	double sum_p_4() const ;
	/// Total energy squared, basically N_<e_i^2>. Extensive.
	double sum_e_squared() const ;
	/// Kinetic energy squared, basically N_<e_{i,kin}^2>. Extensive.
	double sum_ekin_squared() const ;
	/// Interaction energy squared, basically N_<e_{i,int}^2>. Extensive.
	double sum_eint_squared() const ;
	//   ===================================================================================================

	//   ===================================================================================================
	//  Calculation of physical properties
	//   ===================================================================================================
	/// Binder cumulant. 1 - <s^4> / (3 <s^2>). Intensive.
	inline double binder_cumulant() const { return 1 - N_ * sum_s_4() / ( 3.0 * std::pow(sum_s_squared(), 2.0) ) ; } ;
	/// System interaction energy. Extensive.
	double calc_interaction_energy() const;
	/// Interaction energy of particle i.
	double calc_interaction_energy(int i) const;
	/// Returns system energy. Extensive.
	double calc_kinetic_energy() const;
	/// Returns kinetic energy of particle i.
	double calc_kinetic_energy(int i) const;
	/// Returns system energy. Extensive.
	inline double calc_energy() const { return calc_interaction_energy() + calc_kinetic_energy(); };
	/// Energy of particle i. Extensive.
	inline double calc_energy(int i) const { return calc_kinetic_energy(i) + calc_interaction_energy(i); };
	/// Returns temperature. Careful, this function returns ((<p^2>-<p>^2)/m + (<w^2>-<w>^2)/I)/3, not the kinetic energy. Intensive.
	double calc_temperature() const;
	/// Returns spin angular momentum temperature.
	inline double calc_temperature_w() const {
		return (sum_w_squared() - pow(sum_w(), 2.0)/(N_ * I_) )/(N_);
	};
	/// Returns linear momentum temperature.
	inline double calc_temperature_p() const {
		return (sum_p_squared() - topology::norm2(sum_p())/(2 * N_ * m_))/(2 * N_);
	};

	/// Return the plaquette the particle i belongs to. i is in the lower left corner. Only works for lattice-based models.
	std::vector<int> plaquette(int i) const;
	/// Returns vorticity along the plaquette at index.
	double calc_vorticity(int index) const;
	/// Returns the unsigned vortex density (i.e. number of vortices divided by box area).
	double calc_vortexdensity_unsigned() const;
	/// Returns the signed vortex density (i.e. number of positive vortices minus number of negative vortices divided by box area).
	double calc_vortexdensity_signed() const;
	/// Returns total spatial angular momentum of particles
	double calc_space_angular_mom() const ;
	/// Returns spatial angular momentum of the particle with index i
	inline double calc_space_angular_mom(int i) const { return r_[i].get_x() * p_[i].get_y() - r_[i].get_y() * p_[i].get_x(); };
	/// Calculates the mean over nearest neighbors
	/*! For each pair of particles, the function determines
	 * "theta^te_pow * r^r_pow * cos(theta)^cos_pow * sin(theta)^sin_pow * J(r)^J_pow * U'(r)^Up_pow  * U''(r)^Upp_pow"
	 * and adds all the values up. Here, theta is the spin angle difference between the two particles,
	 * r is the distance between the particles.
	 * Can be used to calculate transport coefficients or other diagnostics.
	 */
	double calc_neighbor_mean(double te_pow, double r_pow, double cos_pow, double sin_pow, double J_pow, double Up_pow, double Upp_pow) const;

	/// Calculates the helicity modulus and auxiliary quantities. Output is a vector with entries (Upsilon,H_x,H_y,I_x,I_y)
	std::vector<double> calc_helicity(double beta) const;

	//  ===================================================================================================
	//	Functions implemented by Mathias Hoefler. Redundant and with old group structure
	//  ===================================================================================================
	//	double ang_MSD(xygroup &initial) const;
	//	double ang_MSD_nonsat(const xygroup &initial) const;
	//	double spin_autocorrelation(const xygroup &corr_G) const;


	// ======================================================================================================
	// Calculation of field fluctuations
	// ======================================================================================================
	/// Calculates e^(i q r_i) for particle i.
	inline std::complex<double> calc_eiqr(const topology::Vector2d q, int i) const {
		return std::complex<double>(std::cos(topology::innerproduct(q,get_r(i))),
				std::sin(topology::innerproduct(q,get_r(i))));
	};
	/// Calculates $m_{x,q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_mxq(const topology::Vector2d q, double Mx_0) const;
	/// Calculates $m_{y,q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_myq(const topology::Vector2d q, double My_0) const;
	/// Calculates $w_{q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_wq(const topology::Vector2d q, double W_0) const;
	/// Calculates $e_{q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_eq(const topology::Vector2d q, double E_0) const;
	/// Calculates $theta_{q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_teq(const topology::Vector2d q, double Te_0) const;
	/// Calculates $rho_{q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_rq(const topology::Vector2d q) const;
	/// Calculates $j_{q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::vector<std::complex<double> > calc_jq(const topology::Vector2d q, topology::Vector2d J_0) const;
	/// Calculates $j_{q,L}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_jqpar(const topology::Vector2d q, topology::Vector2d J_0) const;
	/// Calculates $j_{q,T}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_jqperp(const topology::Vector2d q, topology::Vector2d J_0) const;
	/// Calculates $l_{q}$ (see Bissinger PhD thesis) and \link calc_fieldfluct calc_fieldfluct \endlink
	std::complex<double> calc_lq(const topology::Vector2d q, double L_0) const;

	/// Calculates the average of the field fluctuation fluctname. (e.g. for "wq" this returns sum omega_i.)
	/*! For further details on the variable fluctname, see \link calc_fieldfluct calc_fieldfluct \endlink
	 */
	double calc_fieldfluct_average(std::string fluctname, topology::Vector2d q = 0) const;
	/// Calculates the one-particle density associated with the field fluctuation fluctname. (e.g. for "wq" this returns omega_index.)
	/*! For further details on the variable fluctname, see \link calc_fieldfluct calc_fieldfluct \endlink
	 */
	double calc_one_particle_density(int index, std::string fluctname, topology::Vector2d q = 0) const;
	/// Calculates the field fluctuation for the quantity specified in fluctname.
	/*! \fn calc_fieldfluct
	 *  A field fluctuation is a quantity
	 *  \f[
	 *  	a_{\mathbf{q}} = \frac{1}{\sqrt{N}} \sum_j a_j e^{-i \mathbf{q}\cdot\mathbf{r}_j}
	 *  \f]
	 *  for some property \f$a_j\f$ carried by each particle.
	 *
	 *  This function obtains the value of the field fluctuation at all values of \f$\mathbf{q}\f$ stored in the
	 *  vector qvals. The type of the field fluctuation is defined by the value of fluctname.
	 *
	 *  @param[in]    	qvals Vector containing the wavevectors \f$\mathbf{q}\f$
	 *  @param[in]    	fluctname Name of the fluctuation, specifies which one to compute
	 *  			  \parblock
	 *  			  <table>
     *                <caption id="multi_row">Values of fluctname</caption>
     *                  <tr><th> name value <th>\f$a_\mathbf{q}\f$				<th>\f$a_i\f$				<th>Meaning
	 *                	<tr><td>"mxq"		<td>\f$m_{x,\mathbf{q}}\f$			<td>\f$s_{i,x}\f$			<td>Magnetization in x-direction
	 *                	<tr><td>"myq"		<td>\f$m_{y,\mathbf{q}}\f$			<td>\f$s_{i,y}\f$			<td>Magnetization in y-direction
	 *                	<tr><td>"wq"		<td>\f$w_{\mathbf{q}}\f$			<td>\f$\omega_i\f$			<td>Spin angular momentum
	 *                	<tr><td>"eq"		<td>\f$e_{\mathbf{q}}\f$			<td>\f$e_i\f$				<td>Energy density
	 *                	<tr><td>"teq"		<td>\f$\theta_{\mathbf{q}}\f$		<td>\f$\theta_i\f$			<td>Spin angle (not recommended to use)
	 *                	<tr><td>"rq"		<td>\f$\rho_{\mathbf{q}}\f$			<td>\f$1\f$					<td>Density
	 *                	<tr><td>"lq"		<td>\f$l_{\mathbf{q}}\f$			<td>\f$l_i\f$				<td>Spatial angular momentum (not recommended to use)
	 *                	<tr><td>"jparq"		<td>\f$j_{\parallel,\mathbf{q}}\f$	<td>\f$v_{i}^{\parallel}\f$	<td>Longitudinal velocity fluctuation (along \f$\mathbf{q}\f$)
	 *                	<tr><td>"jparq"		<td>\f$j_{\perp,\mathbf{q}}\f$		<td>\f$v_{i}^{\perp}\f$		<td>Transversal velocity fluctuation (perpendicular to \f$\mathbf{q}\f$)
	 *                  <tr><td>Other		<td> -- 							<td> -- 					<td> For any other entry, the return value is set to 0. A warning is printed to std::cerr.
	 *              	</table>
	 *                \endparblock
	 */
	std::vector<std::complex<double> > calc_fieldfluct(const std::vector<topology::Vector2d> qvals, std::string fluctname) const;
	/// Calculates the field fluctuation for the quantity specified in fluctname.
	std::vector<std::complex<double> > calc_fieldfluct_convolution(const std::vector<topology::Vector2d> qvals, std::string fluctname_1, std::string fluctname_2) const;

	// ======================================================================================================
	// Functions to calculate transport coefficients. These do not work properly
	// ======================================================================================================
	/// Calculates \f$\tau\f$, as defined for the xy model. NOT CORRECT FOR THE MOBILE CASE.
	topology::Vector2d calc_tau() const;
	/// Calculates \f$j^e\f$, as defined for the xy model. NOT CORRECT FOR THE MOBILE CASE.
	topology::Vector2d calc_je() const;
	/// Calculates the current for the quantity specified in currentname. currentname = {"tau","je"}. NOT CORRECT FOR THE MOBILE CASE.
	topology::Vector2d calc_current(std::string currentname) const;

	//   ===================================================================================================
	//   Static correlation functions
	//   ===================================================================================================
	/// Calculates the static spin correlation function for a specific particle at index
	/*! More precisely, obtains \<S_index * S_j\>=\<cos(theta_[index]-theta_[j])\> for all j in the sample,
	 *  sums the results and sorts them into bins.
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     *  @param[inout] counts Counts how many particle are found within a bin. As the function is typically
     *                repeatedly to sample many particles, this parameter can be updated for each individual
     *                call to the function. Must be initialized to zero before the first call.
	 */
	std::vector<double> calc_SCF_S_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;

	/// Calculates the static oriented spin correlation function for a specific particle at index
	/*! More precisely, obtains \<cos(theta_[index])*cos(theta_[j])\> for all j in the sample,
	 *  sums the results and sorts them into bins.
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     *  @param[inout] counts Counts how many particle are found within a bin. As the function is typically
     *                repeatedly to sample many particles, this parameter can be updated for each individual
     *                call to the function. Must be initialized to zero before the first call.
	 */
	std::vector<double> calc_SCF_S_oriented_individual( const int index, const std::vector<double> rbin, const double& orientation_angle, std::vector<int>& counts) const;

	/// Calculates the pair distribution function g(r) for a specific particle at index
	/*! More precisely, counts particles j within the interval (rbin[k],rbin[k+1]].
	 *  The result is multiplied by 2 * pi * rbin[k] * dr[k] * rho, with the density rho and the
	 *  bin width dr[k] = rbin[k+1] - rbin[k].
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     */
	std::vector<double> calc_SCF_g_individual( const int index, const std::vector<double> rbin) const;

	/// Calculates the overall pair distribution function g(r)
		/*! Uses calc_SCF_g_individual for number_of_points many randomly chosen particles and
		 *  averages over the result
		 *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
		 *                for an equidistant bin. Particles further from the focal particle than rmax
		 *                are not considered.
		 *  @param[in]    number_of_points Determines how often the function calls calc_SCF_g_individual.
		 */
	std::vector<double> calc_SCF_g( const std::vector<double> rbin, int number_of_points) const;

	/// Calculates the static angle difference correlation function for a specific particle at index
	/*! More precisely, obtains \<theta_[index]-theta_[j]\> for all j in the sample,
	 *  sums the results and sorts them into bins.
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     *  @param[inout] counts Counts how many particle are found within a bin. As the function is typically
     *                repeatedly to sample many particles, this parameter can be updated for each individual
     *                call to the function. Must be initialized to zero before the first call.
     */
	std::vector<double> calc_SCF_anglediff_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;

	/// Calculates the static total energy correlation function for a specific particle at index
	/*! More precisely, obtains \<e(index)*e(j)\> for all j in the sample, with e the total energy of a particle,
	 *  sums the results and sorts them into bins.
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     *  @param[inout] counts Counts how many particle are found within a bin. As the function is typically
     *                repeatedly to sample many particles, this parameter can be updated for each individual
     *                call to the function. Must be initialized to zero before the first call.
	 */
	std::vector<double> calc_SCF_E_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;

	/// Calculates the static kinetic energy correlation function for a specific particle at index
	/*! More precisely, obtains \<e_kin(index)*e_kin(j)\> for all j in the sample, with e_kin the kinetic energy of a particle,
	 *  sums the results and sorts them into bins.
	 */
	std::vector<double> calc_SCF_Ekin_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;

	/// Calculates the static interaction energy correlation function for a specific particle at index
	/*! More precisely, obtains \<e_int(index)*e_int(j)\> for all j in the sample, with e_int the interaction energy of a particle,
	 *  sums the results and sorts them into bins.
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     *  @param[inout] counts Counts how many particle are found within a bin. As the function is typically
     *                repeatedly to sample many particles, this parameter can be updated for each individual
     *                call to the function. Must be initialized to zero before the first call.
	 */
	std::vector<double> calc_SCF_Eint_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;

	/// Calculates the static momentum correlation function for a specific particle at index
	/*! More precisely, obtains \<p_[index]*p[j]\> for all j in the sample (meaning the inner product in this case),
	 *  sums the results and sorts them into bins.
	 */
	std::vector<double> calc_SCF_P_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;
	/// Calculates the static spin momentum correlation function for a specific particle at index
	/*! More precisely, obtains \<w_[index]*w[j]\> for all j in the sample,
	 *  sums the results and sorts them into bins.
	 *  @param[in]    index Particle index.
     *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
     *                for an equidistant bin. Particles further from the focal particle than rmax
     *                are not considered.
     *  @param[inout] counts Counts how many particle are found within a bin. As the function is typically
     *                repeatedly to sample many particles, this parameter can be updated for each individual
     *                call to the function. Must be initialized to zero before the first call.
	 */
	std::vector<double> calc_SCF_W_individual( const int index, const std::vector<double> rbin, std::vector<int>& counts) const;

	/// Calculates the static correlation function specified by name for number_of_points many random particles
	/*! Uses one of the individual SCF calculation functions for number_of_points many randomly chosen particles and
	 *  averages over the result.
	 *
	 *  **Improvement possibilities**
	 *
	 *  - *Case handling*. Case handling for different names follows syntactic simplicity. One could
	 *    rewrite the code to drastically reduce calls to if-cases.
	 *
	 *  - *Information efficiency*. In this function, points within a specific bin are determined. If one calls
	 *    for this function repeatedly, these points are always calculated anew. This is a great loss of
	 *    efficiency and could be mended by more careful code.
	 *
	 *  @param[in]    rbin Bin prescription. Contains bin edges, e.g. rbin = (0,dr,2*dr,...,rmax-dr,rmax)
	 *                for an equidistant bin. Particles further from the focal particle than rmax
	 *                are not considered.
	 *  @param[in]    number_of_points Determines how often the function calls an SCF_individual function.
	 *  @param[in]    name Used to choose a specific correlation function to be calcualted. Options are
	 *  			  \parblock
	 *  			  <table>
     *                <caption id="multi_row">Values of name</caption>
     *                  <tr><th> name value <th>Operation
	 *                	<tr><td>"g"			<td>uses calc_SCF_g_individual, calculates g(r)
	 *                	<tr><td>"anglediff"	<td>uses calc_SCF_anglediff_individual, calculates mean angle difference
	 *                	<tr><td>"S"			<td>uses calc_SCF_S_individual, calculates mean spin alignment
	 *                	<tr><td>"S_par"		<td>uses calc_SCF_S_oriented_individual with the orientation along
	 *                	                        the total magnetization angle
	 *                	<tr><td>"S_perp"	<td>uses calc_SCF_S_oriented_individual with the orientation
	 *                	                        perpendicular to the total magnetization angle
	 *                	<tr><td>"P"			<td>uses calc_SCF_P_individual, calculates mean momentum alignment
	 *                	<tr><td>"W"			<td>uses calc_SCF_W_individual, calculates mean spin momentum correlation
	 *                	<tr><td>"E"			<td>uses calc_SCF_E_individual, calculates mean energy correlation
	 *                	<tr><td>"Ekin"		<td>uses calc_SCF_Ekin_individual, calculates mean kinetic energy correlation
	 *                	<tr><td>"Eint"		<td>uses calc_SCF_Eint_individual, calculates mean interaction energy correlation
	 *                  <tr><td>Other		<td>For any other entry, the return value is set to 0. A warning is printed to std::cerr.
	 *              	</table>
	 *                \endparblock
	 */
	std::vector<double> calc_SCF_averaged( const std::vector<double> rbin, int number_of_points, const std::string name) const;


	//   ===================================================================================================
	//   Time correlation functions
	//   ===================================================================================================
	/// Calculates the spin autocorrelation-function averaged over all indices.
	/*! Computes \f$ \frac{1}{N}\sum_{i=1}^N \mathbf{s}_i^{\texttt{G}} \cdot \mathbf{s}^{\texttt{G\_initial}}\f$,
	 *  where \f$\mathbf{s}_i^{\texttt{G\_initial}}\f$ is the i-th spin in the group G_initial
	 *  and \f$\mathbf{s}_i^{\texttt{G}}\f$ is the i-th spin in the current instance of group
	 *  for which calc_ACF_S is called.
	 */
	double calc_ACF_S(const group& G_initial) const;

	/// Calculates the angle difference autocorrelation-function averaged over all indices.
	/*! Computes \f$ \frac{1}{N}\sum_{i=1}^N (\theta_i^{\texttt{G\_initial}}
	 *  - \theta_i^{\texttt{G}})^2\f$,
	 *  where \f$\theta_i^{\texttt{G\_initial}}\f$ is the i-th spin angle in group G_initial
	 *  and \f$\theta_i^{\texttt{G}}\f$ is the i-th spin angle in the current instance of group
	 *  for which calc_ACF_anglediff is called.
	 */
	double calc_ACF_anglediff(const group& G_initial) const;

	/// Calculates the single-particle autocorrelation-function for some quantity specified by name.
	/*! Computes \f$ \frac{1}{N}\sum_{i=1}^N \big\langle a_i^{\texttt{G\_initial}} \cdot a_i^{\texttt{G}} \big\rangle\f$,
	 *  where \f$a_i\f$ is a quantity defined for each particle individually.
	 *  \f$a_i^{\texttt{G\_initial}}\f$ is then the quantity associated to the i-th particle in the group G_initial,
	 *  while \f$a_i^{\texttt{G}}\f$ the the quantity associated to the i-th particle in the current instance of group
	 *  for which calc_ACF_sp is called.
	 *
	 *  **Improvement possibilities**
	 *
	 *  - *Case handling*. Case handling for different names follows syntactic simplicity. One could
	 *    rewrite the code to drastically reduce calls to if-cases.
	 *
	 *  @param[in]    	G_initial group with which the correlation is to be compared. In most cases, this is the
	 *  				simulated group at a previous time.
	 *  @param[in]    	name Used to choose a specific correlation function to be calcualted. Options are
	 *  			  \parblock
	 *  			  <table>
     *                <caption id="multi_row">Values of name</caption>
     *                  <tr><th> name value <th>Operation
	 *                	<tr><td>"S"			<td>\f$a_i = \mathbf{s}_i\f$, same as calc_ACF_S
	 *                	<tr><td>"Sx"		<td>\f$a_i = s_{i,x}\f$, x-component of spin
	 *                	<tr><td>"anglediff"	<td>\f$a_i = \theta_{i}\f$, averges over \f$(\theta_i^{\texttt{G\_initial}}
	 *                                      	- \theta_i^{\texttt{G}})^2\f$. See calc_ACF_anglediff
	 *                	<tr><td>"Spar"		<td>\f$a_i = s_{i,\parallel}\f$, that is spins oriented along
	 *                	                    	the total magnetization angle
	 *                	<tr><td>"Sperp"		<td>\f$a_i = s_{i,\perp}\f$, that is spins oriented
	 *                	                    	perpendicular to the total magnetization angle
	 *                	<tr><td>"P"			<td>\f$a_i = \mathbf{p}_i\f$, linear momentum
	 *                	<tr><td>"Px"		<td>\f$a_i = \mathbf{p}_{i,x}\f$, x-component of linear momentum
	 *                	<tr><td>"Py"		<td>\f$a_i = \mathbf{p}_{i,y}\f$, y-component of linear momentum
	 *                	<tr><td>"Ppar"		<td>\f$a_i = \mathbf{p}_{i,\parallel}\f$, component of linear momentum
	 *                	                    	parallel to the total magnetization
	 *                	<tr><td>"Pperp"		<td>\f$a_i = \mathbf{p}_{i,\perp}\f$, component of linear momentum
	 *                	                    	perpendicular to the total magnetization
	 *                	<tr><td>"W"			<td>\f$a_i = \omega_{i}\f$, spin momentum
	 *                	<tr><td>"E"			<td>\f$a_i = e_{i}\f$, energy per particle
	 *                	<tr><td>"Ekin"		<td>\f$a_i = e_{\textrm{kin},i}\f$, kinetic energy per particle
	 *                	<tr><td>"Eint"		<td>\f$a_i = e_{\textrm{int},i}\f$, interaction energy per particle
	 *                	<tr><td>"MSD"		<td>Averages over \f$\big(\mathbf{r}_i^{\texttt{G}} - \mathbf{r}_i^{\texttt{G\_initial}}\big)^2\f$.
	 *                                      	Careful, does not take periodic boundary into consideration
	 *                  <tr><td>Other		<td>For any other entry, the return value is set to 0. A warning is printed to std::cerr.
	 *              	</table>
	 *                \endparblock
	 *
	 */
	double calc_ACF_sp(const group& G_initial, const std::string name) const;

	/// Calculates the (q=0)-autocorrelation-function for some quantity specified by name.
	/*! Computes \f$ \frac{1}{N}\sum_{i=1}^N \big\langle a_i^{\texttt{G\_initial}} \cdot a_i^{\texttt{G}} \big\rangle\f$,
	 *  where \f$a_i\f$ is a quantity defined for each particle individually.
	 *  \f$a_i^{\texttt{G\_initial}}\f$ is then the quantity associated to the i-th particle in the group G_initial,
	 *  while \f$a_i^{\texttt{G}}\f$ the the quantity associated to the i-th particle in the current instance of group
	 *  for which calc_ACF_sp is called.
	 *
	 *  **Improvement possibilities**
	 *
	 *  - *Case handling*. Case handling for different names follows syntactic simplicity. One could
	 *    rewrite the code to drastically reduce calls to if-cases.
	 *
	 *  @param[in]    	G_initial group with which the correlation is to be compared. In most cases, this is the
	 *  				simulated group at a previous time.
	 *  @param[in]    	name Used to choose a specific correlation function to be calcualted. Options are
	 *  			  \parblock
	 *  			  <table>
	 *                <caption id="multi_row">Values of name</caption>
	 *                  <tr><th> name value <th>Operation
	 *                	<tr><td>"S"			<td>\f$a_i = \mathbf{s}_i\f$, same as calc_ACF_S
	 *                	<tr><td>"Sx"		<td>\f$a_i = s_{i,x}\f$, x-component of spin
	 *                	<tr><td>"anglediff"	<td>\f$a_i = \theta_{i}\f$, averges over \f$(\theta_i^{\texttt{G\_initial}}
	 *                                      	- \theta_i^{\texttt{G}})^2\f$. See calc_ACF_anglediff
	 *                	<tr><td>"Spar"		<td>\f$a_i = s_{i,\parallel}\f$, that is spins oriented along
	 *                	                    	the total magnetization angle
	 *                	<tr><td>"Sperp"		<td>\f$a_i = s_{i,\perp}\f$, that is spins oriented
	 *                	                    	perpendicular to the total magnetization angle
	 *                	<tr><td>"P"			<td>\f$a_i = \mathbf{p}_i\f$, linear momentum
	 *                	<tr><td>"Px"		<td>\f$a_i = \mathbf{p}_{i,x}\f$, x-component of linear momentum
	 *                	<tr><td>"Py"		<td>\f$a_i = \mathbf{p}_{i,y}\f$, y-component of linear momentum
	 *                	<tr><td>"Ppar"		<td>\f$a_i = \mathbf{p}_{i,\parallel}\f$, component of linear momentum
	 *                	                    	parallel to the total magnetization
	 *                	<tr><td>"Pperp"		<td>\f$a_i = \mathbf{p}_{i,\perp}\f$, component of linear momentum
	 *                	                    	perpendicular to the total magnetization
	 *                	<tr><td>"W"			<td>\f$a_i = \omega_{i}\f$, spin momentum
	 *                	<tr><td>"E"			<td>\f$a_i = e_{i}\f$, energy per particle
	 *                	<tr><td>"Ekin"		<td>\f$a_i = e_{\textrm{kin},i}\f$, kinetic energy per particle
	 *                	<tr><td>"Eint"		<td>\f$a_i = e_{\textrm{int},i}\f$, interaction energy per particle
	 *                	<tr><td>"MSD"		<td>Averages over \f$\big(\mathbf{r}_i^{\texttt{G}} - \mathbf{r}_i^{\texttt{G\_initial}}\big)^2\f$.
	 *                                      	Careful, does not take periodic boundary into consideration
	 *                	<tr><td>Other		<td>For any other entry, the return value is set to 0. A warning is printed to std::cerr.
	 *              	</table>
	 *                \endparblock
	 *
	 */
	double calc_ACF_q0(const group& G_initial, const std::string name) const;


	/// Calculates time-correlation function between two different groups. Fluctuation names must be specified.
	/*! \fn calc_TCF
	 *  Calculates \f$a_{\mathbf{q}}^* b_{\mathbf{\qq}}(t)\f$, or more accurately
	 *  \f$(a_{\mathbf{q}}^{\texttt{G\_initial}})^* b_{\mathbf{\qq}}^{\texttt{G}}\f$,
	 *  where \f$\mathbf{q}\f$ is a wave vector and \f$a_{\mathbf{q}}^{\texttt{G\_initial}}\f$ and
	 *  \f$b_{\mathbf{q}}^{\texttt{G}}\f$ are field fluctuations associated with the group
	 *  \f$\texttt{G\_initial}\f$ and \f$\texttt{G}\f$ (the current instance of group for which
	 *  calc_TCF is called), respectively.
	 *
	 *  Returns a vector whose entries correspond to the wavevectors in qvals.
	 *
	 *  @param[in]    	G_initial group with which the correlation is computed. In most cases, this is the
	 *  				simulated group at a previous time.
	 *  @param[in]    	qvals vector of \f$\mathbf{q}\f$-values for which the product of field fluctuations
	 *                  is calculated
	 *  @param[in]    	fluctname_initial specifies \f$a_{\mathbf{q}}^{\textttt{G_initial}}\f$,
	 *  				the field fluctuation of G_initial. For details, see
	 *  				\link calc_fieldfluct calc_fieldfluct \endlink
	 *
	 *  @param[in]    	fluctname_current specifies \f$b_{\mathbf{q}}^{\textttt{G}}\f$,
	 *  				the field fluctuation in G. For details, see
	 *  				\link calc_fieldfluct calc_fieldfluct \endlink
	 */
	std::vector<std::complex<double> > calc_TCF(const group& G_initial, std::vector<topology::Vector2d> qvals,
			std::string fluctname_initial, std::string fluctname_current) const;


	// ======================================================================================================
	// Calculation of time derivatives and coordinate differences
	// ======================================================================================================
	/// Returns theta (spin angle) time derivative (splitting useful for leapfrog)
	std::vector<double> time_derivative_theta() const;
	/// Returns omega (spin momentum) time derivative (splitting useful for leapfrog)
	std::vector<double> time_derivative_w() const;
	/// Returns r (particle position) time derivative. First N_ entries are x direction, N_+1 to 2N_ is y direction (splitting useful for leapfrog)
	std::vector<double> time_derivative_r() const;
	/// Returns p (linear momentum) time derivative. First N_ entries are x direction, N_+1 to 2N_ is y direction (splitting useful for leapfrog)
	std::vector<double> time_derivative_p() const;

	/// Returns coordinate time derivative (first N_ entries are theta, then r_x, then r_y).
	std::vector<double> time_derivative_coord() const;
	/// Returns momenta time derivative (first N_ entries are omega, then p_x, then p_y).
	std::vector<double> time_derivative_mom() const;
	/// Returns time derivative of the entire group.
	group time_derivative() const;

	/// Returns coordinate difference between this group and another one, with proper care of boundaries. First N_ entries are theta, then r_x, then r_y.
	std::vector<double> coord_diff(const group& G) const;
	/// Same as coord_diff, but adding the difference to an MSD vector.
	void accumulative_MSD(std::vector<double>& MSD, const group& last_G) const;


	//   ===================================================================================================
	//   Addition operator of groups and multiplication operator of a group by a double value
	//   ===================================================================================================
	/// Adds particle entries (used for adding time derivatives and such).
	group& operator+=(const group& G);
	/// Multiplies particles by constant (used for adding time derivatives and such).
	group& operator*=(const double a);

protected:

	/// Type of the group. Can be "xy" for the XY model, "mxy" for the mobile XY model, "fmxy" for a mobile XY model frozen in place, "vm" for the Vicsek model and "fvm" for the frozen (static) Vicsek model
	std::string group_type_;

	/// Size of the group.
	int N_;
	/// Square root of the group size (often useful).
	int sqrtN_;

	/// Size of the box (some functions only defined for square boxes yet).
	topology::Vector2d L_;

	/// Spin inertia
	double I_ = 1;
	/// Mass
	double m_ = 1;

	/// Nearest neighbor interaction strength
	double J_ = 1;
	/// Spatial repulsion interaction strength
	double U_ = 1;

	/// Interaction cutoff radius.
	double cutoff_ = 1;

	std::vector<topology::Vector2d> r_;			///< Particle positions in the group
	std::vector<topology::Vector2d> p_;			///< Particle momenta/velocities in the group
	std::vector<double> theta_;					///< Particle spin angles in the group
	std::vector<double> w_;						///< Particle spin momenta in the group

	/// Partition (cell list) for neighborhood interaction
	partition partition_;


	std::string nb_rule_;				///< Neighbor calculation rule. Possible values: "bruteforce", "all", "ur" for full, (partition with) all and (partition with) upper right neighbors.
	double nb_mult_factor_;				///< If neighbor rule leads to double counting, this factor has to be .5, otherwise 1.

	neighbor_list* nb_list_;			///< neighbor-list

	/// lattice type. (type 's': square, type 't': trigonal, type 'n': none (mxy model etc))
	char lattice_type_;

	/// Vicsek model parameter eta: Angle for random noise.
	double vm_eta_;

	/// Vicsek model parameter v: Streaming velocity
	double vm_v_;



};

/// Addition operator. Adds all particle entries of two groups
inline group operator+(const group& G, const group& G2){
	return  group(G) += G2;
};


/// Right multiplication operator, multiplies all group elements by a scalar
inline group operator*(const group& G, const double a){
	return group(G) *= a;
};

/// Left multiplication operator, multiplies all group elements by a scalar
inline group operator*(const double a, const group& G){
	return group(G) *= a;
};

#endif /* GROUP_H_ */
