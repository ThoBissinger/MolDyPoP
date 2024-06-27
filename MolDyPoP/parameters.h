/*! \file parameters.h
 *
 *  \brief Header-File to class declaration of parameters. Introduces parameters,
 *  the data structure associated with input data that governs the simulation
 *  run.
 *
 *  Contains functionalities for reading input files, setting and sharing parameter values.
 *
 *  \date Created: 2019-04-12
 *  \date Last Updated: 2023-08-02
 *
 *
 *  \class parameters.
 *
 *  \brief Contains the run parameters of a simulation.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: 2019-04-12
 *  \date Last Updated: 2023-08-02
 *
 */


#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "computations.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>


class parameters {
protected:
	/*! \brief Type of system.
	 *
	 *   Possible values:
	 *   \parblock
	 *   <table>
	 *   <caption id="multi_row">Values of system </caption>
	 *   <tr><th> system value	<th>System chosen
	 *   <tr><td>"xy"			<td>XY model
	 *   <tr><td>"mxy"			<td>mobile XY (MXY) model
	 *   <tr><td>"fmxy"			<td>disordered XY (DXY) model / frozen mobile XY model
	 *   <tr><td>"vm"			<td>Vicsek model (not fully tested)
	 *   <tr><td>"fvm"			<td>frozen Vicsek model (not fully tested)
	 *   </table>
	 *   \endparblock
	 */
	std::string system_ = "xy";
	/*! \brief Run mode for the simulation.
	 *
	 *   Possible values:
	 *   <table> <caption id="multi_row">Values of mode</caption>
	 *   <tr><th> mode value 		<th>What happens during run?
	 *   <tr><td>"none"				<td>does nothing
	 *   <tr><td>"test"				<td>like non, unless user specified. See main.cpp
	 *   <tr><td>"integrate"		<td>performs time integration, possibly with storing
	 *   <tr><td>"integrate_cont"	<td>picks up an aborted integration and continues (relevant in older version, not recommended)
	 *   <tr><td>"samp"				<td>samples a list of stored snapshots (has to be provided)
	 *   <tr><td>"equilibrate"		<td>similar to integrate, but performs equilibration run only
	 *   </table>
	 */
	std::string mode_ = "none";
	/*! \brief Initialization mode.
	 *
	 *   Possible values:
	 *   <table> <caption id="multi_row">Values of init_mode</caption>
	 *   <tr><th> init_mode value 	<th>initialization
	 *   <tr><td>"random"			<td>all particles placed and oriented randomly
	 *   <tr><td>"aligned"			<td>particles placed randomly, but oriented along some random direction
	 *   <tr><td>"file"				<td>reads initial configuration from file
	 *   </table>
	 *
	 *   For initialization from file, the path has to be specified in the parameter init_file_
	 */
	std::string init_mode_ = "random";
	std::string init_file_ = "input/init_snap.in"; 		///< Path to initialization file. Only relevant if init_mode is set to "file"
	double init_kT_ = -1; 								///< Initialization temperature, temperature at which simulation should start. Has no effect when chosen negative.
	double init_random_displacement_ = 0; 				///< Random particle displacement at initialization (for scaled initialization only)
	double init_random_angle_ = 0; 						///< Random angle displacement at initialization (for scaled initialization only)


	std::string job_id_ = ""; 							///< Job identifier. If specified, data will be stored in data_'job_id'.out etc. Otherwise, it is just the default data.out
	std::string outfilename_ = "output/data.out";		///< Name of outfile.
	int N_ = 256;										///< Number of particles
	int sqrtN_ = 16;									///< Square root of the number of particles
	double dof_;										///< Degrees of freedom. Calculated automatically from N_ and the system_.
	double L_ = 16;										///< Length of simulation box
	double dt_ = 0.01;									///< Integration time step
	double Tmax_ = 100;									///< Runtime for sampling
	double samplestart_ = 0;							///< Time at which sampling starts in the integration run
	double samplestep_ = 1;								///< Time separation of two samples
	double av_time_spacing_ = 1e2; 						///< Minimum time spacing for average calculations in sampling
	int Nsamp_ = 30;									///< Number of samples computed in sampling. Overwrites samplestep.
	int randomseed_ = 1;								///< Random seed
	double kT_ = .89;									///< Temperature (units of energy)
	double I_ = 1;										///< Spin inertia
	double m_ = 1;										///< Mass
	double J_ = 1;										///< Spin interaction strength
	double U_ = 1;										///< Spatial interaction strength
	double cutoff_ = 1; 								///< Cutoff length for interaction potential (default value 1)
	/*! \brief In static system: Type of lattice.
	 *
	 *  Possible values:
	 *   <table> <caption id="multi_row">Values of lattice_type</caption>
	 *   <tr><th> lattice_typ value <th>lattice
	 *   <tr><td> 's'				<td>square lattice (default)
	 *   <tr><td> 't'				<td>trigonal lattice
	 *   </table>
	 */
	char lattice_type_ = 's';
	double activity_ = 0;								///< Activity in the mobile XY model (stream velocity)
	double vm_v_ = 0;									///< Vicsek model streaming velocity.
	double vm_eta_ = 0;									///< Vicsek model noise strength (between 0 and 1, actual noise interval will be [-eta/2*Pi,eta/2 * Pi]).

	/*!  \brief Integrator used for equilibration.
	 *
	 *   Possible values:
	 *   <table> <caption id="multi_row">Values of eq_integrator_type_</caption>
	 *   <tr><th> eq_integrator_type_ value	<th>integrator used
	 *   <tr><td> "lf"						<td>leap-frog scheme
	 *   <tr><td> "rk4"						<td>4th order Runge-Kutta scheme
	 *   <tr><td> "langevin"				<td>Langevin-thermostatted time-evolution (stochastic, not tested)
	 *   <tr><td> "nh"						<td>Nosé-Hoover thermostatted integrator (not extensively tested)
	 *   <tr><td> "np"						<td>Nosé-Poincaré thermostatted integrator (possibly bugged, not extensively tested)
	 *   <tr><td> "mc"						<td>Monte-Carlo integration: Has not been implemented yet
	 *   </table>
	 *
	 *   TODO: Langevin dynamics uses same parameters as sampling integrator. Potentially
	 *         different parameters for equilibration may be implemented
	 *
	 *   TODO: Monte Carlo does not work yet.
	 */
	std::string eq_integrator_type_;
	/*! \brief Equilibration mode.
	 *
	 *  Possible values:
	 *   <table> <caption id="multi_row">Values of eq_mode_</caption>
	 *   <tr><th> eq_mode_ value	<th>mode of equilibration
	 *   <tr><td> "anneal"			<td>reaches the desired temperature with simulated annealing (default)
	 *   <tr><td> "berendsen"		<td>reaches the desired temperature with Berendsen thermostat
	 *   <tr><td> "brownian"		<td>reaches the desired temperature by repeatedly setting it (hard Brownian thermostat)
	 *   </table>
	 *
	 *   TODO: Add options for Langevin thermostat.
	 *
	 *   TODO: Equilibrate without temperature, keep energy fixed.
	 */
	std::string eq_mode_ = "anneal";
	/*! \brief Break condition for equilibration.
	 *
	 *  Possible values:
	 *   <table> <caption id="multi_row">Values of eq_breakcond_</caption>
	 *   <tr><th> eq_breakcond_ value	<th>break condition
	 *   <tr><td> "temperature"			<td>breaks when the desired temperature has been maintained for a certain time (default)
	 *   <tr><td> "time"				<td>breaks when an equilibration time has passed and checks temperature,
	 *   									continues equilibration if temperature is not reached
	 *   <tr><td> "time_hard"			<td>breaks when an equilibration time has passed, regardless of temperature
	 *   <tr><td> "any"					<td>breaks if the time or the temperature condition is met
	 *   </table>
	 */
	std::string eq_breakcond_ = "temperature";
	double eq_Tmax_ = 100;								///< Maximum equilibration time
	double eq_agreement_threshold_ = 1e-2; 				///< Temperature agreement threshold to check successful equilibration
	double eq_av_time_ = 0; 							///< Averaging time in equilibration (set to 0 for no averaging)
	double tau_berendsen_;								///< tau coefficient of the Berendsen thermostat
	double eq_anneal_rate_ = .999; 						///< Annealing rate
	double eq_anneal_step_ = 1e1; 						///< Annealing time step
	double eq_Tprintstep_ =  std::numeric_limits<double>::quiet_NaN();							///< Temperature print interval. Prints intermediate data when equilibration reaches T
	double eq_brownian_kT_omega_ = -1;					///< Temperature for spin momentum thermostat during equilibration in a Brownian dynamics simulation (applied when non-negative)
	double eq_brownian_kT_p_ = -1;						///< Temperature for momentum thermostat during equilibration in a Brownian dynamics simulation (applied when non-negative)
	double eq_brownian_timestep_ = 1.00;				///< Step length for Brownian timestep during equilibration (in units of system time)
	bool eq_sampswitch_ = 0;							///< switch controlling whether sampling data is stored in equilibration
	int eq_Nsamp_ = 100;								///< Number of time samplings performed during equilibration. Only has an effect if eq_sampswitch = 1
	std::string eq_samp_time_sequence_ = "lin";			///< Sampling time sequence. 'lin' or 'log'


	/*! \brief Integrator used for sampling.
	 *
	 *   Possible values:
	 *   <table> <caption id="multi_row">Values of eq_integrator_type_</caption>
	 *   <tr><th> eq_integrator_type_ value	<th>integrator used
	 *   <tr><td> "lf"						<td>leap-frog scheme
	 *   <tr><td> "rk4"						<td>4th order Runge-Kutta scheme
	 *   <tr><td> "langevin"				<td>Langevin-thermostatted time-evolution (stochastic, not tested)
	 *   <tr><td> "nh"						<td>Nosé-Hoover thermostatted integrator (not extensively tested)
	 *   <tr><td> "np"						<td>Nosé-Poincaré thermostatted integrator (possibly bugged, not extensively tested)
	 *   <tr><td> "mc"						<td>Monte-Carlo integration: Has not been implemented yet
	 *   </table>
	 *
	 *   TODO: Monte Carlo does not work yet.
	 */
	std::string sample_integrator_type_;
	std::string ensemble_;								///< Desired ensemble. "nvt" only works with certain integrators. "nve" is good standard.
	double nhnp_pi_ = 0; 								///< Nose-Hoover or Nose-Poincare pi value
	double nhnp_Q_ = -1; 								///< Nose-Hoover or Nose-Poincare Q value. Default value negative, only used if positive. Method using nhnp_tau is more recommended
	double nhnp_tau_ = 0.01;							///< Nose-Hoover or Nose-Poincare tau value (time scale that determines Q by Q = g * kT * tau^2). Order of a microscopic collision time.
	double nh_eta_ = 0; 								///< Nose-Hoover eta value
	double np_s_ = 1; 									///< Nose-Poincare s value (initial time scaling factor)
	double mc_steplength_theta_ = 0.1;					///< Maximal trial step length of MC algorithm for angle theta
	double brownian_kT_omega_ = -1;						///< Temperature for spin momentum thermostat in a Brownian dynamics simulation (applied when non-negative)
	double brownian_kT_p_ = -1;							///< Temperature for momentum thermostat in a Brownian dynamics simulation (applied when non-negative)
	double brownian_timestep_ = 1.00;					///< Step length for Brownian timestep (in units of system time)
	double gamma_ld_p_ = 0;								///< Damping rate gamma of the Langevin dynamics. Associated to the linear angular momentum p.
	double gamma_ld_om_ = 0;							///< Damping rate gamma of the Langevin dynamics. Associated to the spin angular velocity omega.
	double mc_steplength_r_ = 0.1;						///< Maximal trial step length of MC algorithm for distance r

	std::string sampling_time_sequence_ = "lin";		///< Specifies how sampling should be performed. Values are lin (linear sequence of sampling times) and log (logarithmic sequence of sampling times).
	int N_rbin_;										///< Number of bins for r values
	double min_binwidth_r_;								///< Minimal width of r bins (should be around 1)
	std::string qbin_type_ = "mult";					///< Type of qbin. "all" uses all possible q values on the grid, "mult" uses only integer multiples of \f$2\pi/L\f$ (recommended)
	double qmax_ = 2 * M_PI;							///< Maximum value for q in bin. Default is \f$2\pi\f$ (used when qbin_type = 'all')
	double min_binwidth_q_ = .015;						///< Minimal width of q bins, in fractions of 2pi / L
	int N_qbin_;										///< Number of bins for q values (used when qbin_type = 'mult')
	std::vector<double> rbin_;							///< Binning values for the vector r. Set internally.
	std::vector<double> qbin_;							///< Binning values for the vector q. Set internally.
	double qfullmax_;									///< q value for full coverage (everything beyond is randomly sampled)
	int qsamps_per_bin_;								///< Number of q samples on a circle.
	int n_rsamps_ = 30; 								///< Number of (randomly chosen) r sampling points for the calculation of correlation functions in position space (like spin correlation functions or mean squared displacements etc.)
	bool print_snapshots_ = false;						///< Specifies if snapshots should be printed
	bool on_fly_sampling_ = true;						///< Specifies if sampling should be performed on-the-fly.
	std::string output_folder_ = "output";							///< Path to output folder. Added by David Stadler.
	std::string snap_overview_file_ = "output/snapshot_overview.out"; ///< File where the names of all snapshots are stored. Only relevant when run in sampling mode.


public:
// constructors
	/// Empty constructor
	parameters(){};

	////////////////////////////////////////////////////////////////////////////////////
	// Filling possibilities
	////////////////////////////////////////////////////////////////////////////////////
	/// Reads parameters from an infile.
	int read_from_file(std::ifstream& infile);
	/// Corrects potentially wrong input
	int correct_values(std::ofstream& stdoutfile);
	/// Initializes both bins
	void initialize_bins();
	/*! \brief Initializes the qbin.
	 *
	 *  Also initializes qfullmax_ so that the first bin that is randomly
	 *  selected contains more than 1.5 * qsamps_per_bin_ q-values.
	 */
	void initialize_qbin();
	/// Initializes the rbin
	void initialize_rbin(int N_rbin);
	/// tau_berendsen can be scaled with this.
	void scale_tau(double scale_factor);



	////////////////////////////////////////////////////////////////////////////////////
	// get functions
	////////////////////////////////////////////////////////////////////////////////////
	inline std::string system() const { return system_; }; ///< Returns system_
	inline std::string mode() const { return mode_; }; ///< Returns mode_
	inline std::string job_id() const { return job_id_; }; ///< Returns job_id_
	inline std::string outfilename() const { return outfilename_; }; ///< Returns outfilename_

	inline int N() const { return N_; }; ///< Returns N_
	inline int sqrtN() const { return sqrtN_; }; ///< Returns sqrtN_
	inline double dof() const { return dof_; }; ///< Returns dof_
	inline double L() const { return L_; }; ///< Returns L_
	inline double dt() const { return dt_; }; ///< Returns dt_
	inline double Tmax() const { return Tmax_; }; ///< Returns Tmax_

	inline double samplestart() const { return samplestart_; }; ///< Returns samplestart_
	inline double samplestep() const { return samplestep_; }; ///< Returns samplestep_
	inline double av_time_spacing() const { return av_time_spacing_; }; ///< Returns av_time_spacing_
	inline int Nsamp() const { return Nsamp_; }; ///< Returns Nsamp_
	inline int randomseed() const { return randomseed_; }; ///< Returns randomseed_
	inline double kT() const { return kT_; }; ///< Returns kT_
	inline double I() const { return I_; }; ///< Returns I_
	inline double m() const { return m_; }; ///< Returns m_
	inline double J() const { return J_; }; ///< Returns J_
	inline double U() const { return U_; }; ///< Returns U_
	inline double cutoff() const { return cutoff_; }; ///< Returns cutoff_
	inline char lattice_type() const { return lattice_type_; }; ///< Returns lattice_type_
	inline double activity() const { return activity_; };///< Returns activity_
	inline double vm_v() const { return vm_v_; }; ///< Returns vm_v_
	inline double vm_eta() const { return vm_eta_; }; ///< Returns vm_eta_

	inline std::string init_mode() const { return init_mode_; }; ///< Returns init_mode_
	inline std::string init_file() const { return init_file_; }; ///< Returns init_file_
	inline double init_kT() const { return init_kT_; }; ///< Returns init_kT_
	inline double init_random_displacement() const { return init_random_displacement_; }; ///< Returns init_random_displacement_
	inline double init_random_angle() const { return init_random_angle_; }; ///< Returns init_random_angle_


	inline std::string eq_mode() const { return eq_mode_; }; ///< Returns eq_mode_
	inline std::string eq_integrator_type() const { return eq_integrator_type_; }; ///< Returns eq_integrator_type_
	inline double eq_Tmax() const { return eq_Tmax_; }; ///< Returns eq_Tmax_
	inline std::string eq_breakcond() const { return eq_breakcond_; }; ///< Returns eq_breakcond_
	inline double eq_agreement_threshold() const { return eq_agreement_threshold_; }; ///< Returns eq_agreement_threshold_
	inline double eq_av_time() const { return eq_av_time_; }; ///< Returns eq_av_time_
	inline double tau_berendsen() const { return tau_berendsen_; }; ///< Returns tau_berendsen_
	inline double eq_anneal_rate() const { return eq_anneal_rate_; }; ///< Returns eq_anneal_rate_
	inline double eq_anneal_step() const { return eq_anneal_step_; }; ///< Returns eq_anneal_step_
	inline double eq_Tprintstep() const { return eq_Tprintstep_; }; ///< Returns eq_Tprintstep_
	inline double eq_brownian_kT_p() const { return eq_brownian_kT_p_; }; ///< Returns eq_brownian_kT_p_
	inline double eq_brownian_timestep() const { return eq_brownian_timestep_; }; ///< Returns eq_brownian_timestep_
	inline double eq_brownian_kT_omega() const { return eq_brownian_kT_omega_; }; ///< Returns eq_brownian_kT_omega_
	inline bool eq_sampswitch() const { return eq_sampswitch_; }; ///< Returns eq_sampswitch_
	inline int eq_Nsamp() const { return eq_Nsamp_; }; ///< Returns eq_Nsamp_
	inline std::string eq_samp_time_sequence() const { return eq_samp_time_sequence_; }; ///< Returns eq_samp_time_sequence_


	inline std::string sample_integrator_type() const { return sample_integrator_type_; }; ///< Returns sample_integrator_type_
	inline std::string ensemble() const { return ensemble_; }; ///< Returns ensemble_
	inline double nhnp_pi() const { return nhnp_pi_; }; ///< Returns nhnp_pi_
	inline double nhnp_Q() const { return nhnp_Q_; }; ///< Returns nhnp_Q_
	inline double nhnp_tau() const { return nhnp_tau_; }; ///< Returns nhnp_tau_
	inline double nh_eta() const { return nh_eta_; }; ///< Returns nh_eta_
	inline double np_s() const { return np_s_; }; ///< Returns np_s_
	inline double mc_steplength_theta() const { return mc_steplength_theta_; }; ///< Returns mc_steplength_theta_
	inline double brownian_kT_omega() const { return brownian_kT_omega_; }; ///< Returns brownian_kT_omega_
	inline double brownian_kT_p() const { return brownian_kT_p_; }; ///< Returns brownian_kT_p_
	inline double brownian_timestep() const { return brownian_timestep_; }; ///< Returns brownian_timestep_
	inline double gamma_ld_om() const { return gamma_ld_om_; }; ///< Returns gamma_ld_om_
	inline double gamma_ld_p() const { return gamma_ld_p_; }; ///< Returns gamma_ld_p_
	inline double mc_steplength_r() const { return mc_steplength_r_; }; ///< Returns mc_steplength_r_


	inline std::string sampling_time_sequence() const { return sampling_time_sequence_; }; ///< Returns sampling_time_sequence_
	inline int N_rbin() const { return N_rbin_; }; ///< Returns N_rbin_
	inline int N_qbin() const { return N_qbin_; }; ///< Returns N_qbin_
	inline double min_binwidth_r() const { return min_binwidth_r_; }; ///< Returns min_binwidth_r_
	inline std::string qbin_type() const { return qbin_type_; }; ///< Returns qbin_type_
	inline double qmax() const { return qmax_; }; ///< Returns qmax_
	inline double min_binwidth_q() const { return min_binwidth_q_; }; ///< Returns min_binwidth_q_
	inline std::vector<double> rbin() const { return rbin_; }; ///< Returns rbin_
	inline std::vector<double> qbin() const { return qbin_; }; ///< Returns qbin_
	inline int qsamps_per_bin() const { return qsamps_per_bin_; }; ///< Returns qsamps_per_bin_
	inline int n_rsamps() const { return n_rsamps_; }; ///< Returns n_rsamps_
	inline double qfullmax() const { return qfullmax_; }; ///< Returns qfullmax_
	inline bool print_snapshots() const { return print_snapshots_; }; ///< Returns print_snapshots_
	inline bool on_fly_sampling() const { return on_fly_sampling_; }; ///< Returns on_fly_sampling_
	inline std::string snap_overview_file() const { return snap_overview_file_; }; ///< Returns snap_overview_file_
	inline std::string output_folder() const { return output_folder_; }; ///< Returns output_folder_ (added by David Stadler)



};
#endif /* PARAMETERS_H_ */
