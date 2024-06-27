/*! \file integrator.h
 *
 *  \brief Header-file to class declaration of integrator. Introduces the integrator,
 *  the data structure associated with discrete time evolution. Incorporates multiple
 *  ODE solvers, deterministic as well as stochastic.
 *
 *  \author Thomas Bissinger, with contributions by Mathias Höfler

 *  \date Created: 2019-04-12
 *  \date Last Updated: 2023-08-06
 *
 *  \class integrator
 *  \brief Defines various integration methods for groups. Also includes thermostats.
 *  Integrators include: fourth order Runge-Kutta (rk4) and Leapfrog (lf)
 *
 *  \author Thomas Bissinger, RK4 implementation by Mathias Hoefler
 *
 *  \date Created: 2019-04-12
 *  \date Last Updated: 2023-08-06
 *
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H


#include <iostream>
#include <fstream>
#include <vector>
#include "topology.h"
#include "computations.h"
#include "parameters.h"
#include "group.h"


#include <chrono>

class integrator {
private:
	/*! \brief Type of the integrator.
	 *
	 *   Possible values:
	 *   <table> <caption id="multi_row">Values of integrator_type_</caption>
	 *   <tr><th> integrator_type_ value	<th>integrator used
	 *   <tr><td> "lf"						<td>leap-frog scheme
	 *   <tr><td> "rk4"						<td>4th order Runge-Kutta scheme
	 *   <tr><td> "langevin"				<td>Langevin-thermostatted time-evolution (stochastic, not tested)
	 *   <tr><td> "nh"						<td>Nosé-Hoover thermostatted integrator (not extensively tested)
	 *   <tr><td> "np"						<td>Nosé-Poincaré thermostatted integrator (possibly bugged, not extensively tested)
	 *   <tr><td> "mc"						<td>Monte-Carlo integration: Has not been implemented yet
	 *   </table>
	 *
	 *   TODO: NH and NP are not tested and have had numerical problems in some instances already.
	 *         Have to be tested and potentially corrected.
	 *
	 *   TODO: Langevin dynamics has been tested and should work without problems, but has not
	 *         been tested extensively.
	 *
	 *   TODO: Monte Carlo does not work yet.
	 */
	std::string integrator_type_;

	double dt_; 						///< Integration time step
	const std::vector<double> a_ {1.0/2.0, 1.0/2.0, 1.0}; ///< RK4 vector components
	const std::vector<double> b_ {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}; ///< RK4 matrix components


	double kT_;							///< Required for any thermostat
	double kT_te_;						///< Optional: Temperature associated to spin degrees of freedom
	double kT_r_;						///< Optional: Temperature associated to positional degrees of freedom
	double activity_;					///< Activity in the mobile XY model

	double g_;							///< Required for any thermostat
	double Q_;							///< Required for any thermostat
	double H_0_;						///< Required in Nose-Poincare
	double pi_;							///< Required in Nose-Poincare and Nose-Hoover
	double eta_;						///< Required in Nose-Hoover
	double s_;							///< Required in Nose-Poincare (time scaling factor)
	double tau_;						///< Not required, but typically used for setting Q_

	double gamma_ld_om_;				///< Langevin gamma associated to spin degrees of freedom
	double gamma_ld_p_;					///< Langevin gamma associated to spatial degrees of freedom


	double mc_steplength_theta_;		///< Maximal trial step length of MC algorithm for angle theta
	double mc_steplength_r_;			///< Maximal trial step length of MC algorithm for distance r

public:
//   ===================================================================================================
//   Constructors and intialization
//   ===================================================================================================//	/////////////////////////////

/// Empty constructor
integrator(){};
/// Constructor with given time steps
integrator(double dtin, std::string type);
/// Sets integrator for equilibration run
void integrator_eq(const parameters& par);
/// Sets integrator for sampling run
void integrator_sample(const parameters& par);

/*! \brief Initializes integrator from given parameters and a given group.
 *
 *  The group is needed to initialize some further values like the reference energy H_0
 *  in a Nosé-Poincaré scheme. Makes a call to initialize_parameters for the
 *  parameter readout.
 */
void initialize(const parameters& par, const group& G);
/// Initializes integrator from given parameters.
void initialize_parameters(const parameters& par);



//   ===================================================================================================
//   Simple get-functions
//   ===================================================================================================//	/////////////////////////////
/// Returns dt
inline double get_dt() const { return dt_; } ;
/// Returns H_0_
inline double get_H_0() const { return H_0_; } ;
/// Returns pi_
inline double get_pi() const { return pi_; } ;
/// Returns eta_
inline double get_eta() const { return eta_; } ;
/// Returns s_
inline double get_s() const { return s_; } ;

//   ===================================================================================================
//   Thermostats
//   ===================================================================================================//	/////////////////////////////
/// Implements the Berendsen thermostat.
double berendsen_thermostat(group& G, double T_desired, double berendsen_tau) ;

//   ===================================================================================================
//   Solvers
//   ===================================================================================================//	/////////////////////////////
/*! \brief Performs an integration time step.
 *
 *  Depending on integrator_type_, different solvers are called applied.
 *  For details, see also definition of integrator_type_.
 *
 *  Not every method requires a variable momdot, yet some make use
 *  of a staggered time grid or the velocity at a previous time step,
 *  these need momdot.
 */
group integrate(const group& G, std::vector<double>& momdot) ;

/// Implements Vicsek model time-evolution rules.
group vm_rule(const group& G) ;

/// Implements the fourth-order Runge-Kutta method
group rk4(const group& G) ;

/// Implements leapfrog solver. Can get initial momdot for reduced computation time.
group leapfrog(const group& G, std::vector<double>& momdot) ;
/// Implements leapfrog solver with added activity. Can get initial momdot for reduced computation time.
group leapfrog_active(const group& G, std::vector<double>& momdot, double activity) ;



/*! \brief Nosé-Poincare integration solved via velocity Verlet.
 *
 *  For some details, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577381/, jchemphys 2008 128(24), Kleinerman et al., eq (21-28)
 */
group np(const group& G, double kT, double& s, double& pi, double Q, double H_0) ;

/*! \brief Nosé-Hoover integration solved via velocity Verlet.
 *
 *  For some details, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577381/, jchemphys 2008 128(24), Kleinerman et al., eq (6-17)
 */
group nh(const group& G) ;

/*! \brief Nosé-Poincare integration solved via velocity Verlet.
 *
 *  For some details, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577381/, jchemphys 2008 128(24), Kleinerman et al., eq (21-28)
 */
inline group np(const group& G) ;


/*! \brief Langevin dynamics solver.
 *
 *  Can get initial momdot for reduced computation time. Details in Allen and Tildesley, "Computer Simulation of Liquids", Chapter 12
 */
group langevin(const group& G, std::vector<double>& momdot) ;

/*! \brief Langevin dynamics solver including active motion of particles along spin.
 *
 *  Can get initial momdot for reduced computation time. Details in Allen and Tildesley, "Computer Simulation of Liquids", Chapter 12
 */
group langevin_active(const group& G, std::vector<double>& momdot, double activityc ) ;

// TODO mc
// template<class GROUP>
// GROUP mc(const GROUP& G) ; ///< Monte-Carlo integration solver (not working yet)


/// The Hamiltonian for Nosé-Hoover and Nosé-Poincaré schemes
inline double NoseHamiltonian(group& G, double kT, double s , double pi, double Q) const {
	return G.calc_kinetic_energy() / (s * s) + G.calc_interaction_energy() + .5 * pi * pi / Q + 2 * G.get_N() * kT * std::log(s);
};

/// Adds activity by moving the particles along spin direction times activity times dt
void add_activity(group& G, double activity) const ;



};
#endif /* INTEGRATOR_H_ */
