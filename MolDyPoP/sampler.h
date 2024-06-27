/*! \file sampler.h
 *
 *  \brief Header-file to class declaration of sampler. Introduces a data structure
 *  calculating and storing different properties of the system during the runtime.
 *
 *  This data structure is central to sampling. See description of class sampler
 *
 *
 *  \author Thomas Bissinger

 *  \date Created: 2020-01-20
 *  \date Last Updated: 2023-08-06
 *
 *
 *  \class sampler
 *  \brief Stores and handles all data sampling performed during a run or in a later diagnostic.
 *
 *  Has a wide variety of quantities that may be of physical interest during the run.
 *  Can compute these quantities when called upon and stores them in vecotrs. Most of
 *  these vectors can get time averaged at the end of the run.
 *
 *  All results are printed in matlab-readable form and can be further processed with matlab code.
 *
 *  TODO: Alternate output than matlab, ideally binary or CSV.
 *
 * 	\date Created on 2020-01-20
 *
 *  \author Thomas Bissinger
 *
 *
 */


#ifndef SAMPLER_H
#define SAMPLER_H



#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "topology.h"
#include "parameters.h"
#include "inputoutput.h"
#include "group.h"

#include <chrono> // Not necessary in code, just for testing.

class sampler {
protected:
	parameters par_;												///< Simulation parameters

	int NumberOfSwitches_ = 15;										///< Total number of switches. Could be useful.

	int nsamp_ = 0;													///< Number of samples (for time averages later on)
	int nsnap_ = 0;													///< Number of snapshots

	std::vector<double> averaging_times_;							///< Stores sampling times (for averaging)
	std::vector<double> TCF_times_;									///< Stores sampling times (for time correlation function)

	std::vector<double> H_;											///< Stores energies
	std::vector<double> H_2_;										///< Stores individual energies squared (<e_i^2>)
	std::vector<double> Hkin_2_;									///< Stores individual kinetic energies squared
	std::vector<double> Hint_2_;									///< Stores individual interaction energies squared
	std::vector<double> W_;											///< Stores momenta
	std::vector<double> W_2_;										///< Stores momenta squared
	std::vector<topology::Vector2d> M_;								///< Stores magnetizations
	std::vector<double> M_2_;										///< Stores magnetizations squared
	std::vector<double> M_4_;										///< Stores magnetizations to the fourth power
	std::vector<double> absM_;										///< Stores absolute magnetizations
	std::vector<double> M_angle_;									///< Stores angle of the magnetization (with respect to the spin x-axis)
	std::vector<double> Theta_;										///< Stores mean theta (with respect to the spin x-axis)
	std::vector<double> Theta_2_;									///< Stores mean theta^2 (with respect to the spin x-axis)
	std::vector<double> Theta_4_;									///< Stores mean theta^4 (with respect to the spin x-axis)
	std::vector<double> Theta_rel_to_M_;							///< Stores mean theta relative to orientation of magnetization
	std::vector<double> Theta_rel_to_M_2_;							///< Stores mean theta^2 relative to orientation of magnetization
	std::vector<double> Theta_rel_to_M_4_;							///< Stores mean theta^2 relative to orientation of magnetization
	std::vector<double> temperature_;								///< Stores temperature
	std::vector<double> temperature_squared_;						///< Stores temperature squared
	std::vector<double> temperature_omega_;							///< Stores spin momentum temperature
	std::vector<double> temperature_omega_squared_;					///< Stores spin momentum temperature squared
	std::vector<double> temperature_p_;								///< Stores linear momentum temperature
	std::vector<double> temperature_p_squared_;						///< Stores linear momentum temperature squared
	std::vector<topology::Vector2d> P_;								///< Stores momentum components
	std::vector<double> P_2_;										///< Stores momentum squared
	std::vector<double> P_4_;										///< Stores momentum to the fourth power

	std::vector<double> H_x_;										///< H_x as defined for the derivation of the helicity Upsilon
	std::vector<double> H_y_;										///< H_y as defined for the derivation of the helicity Upsilon
	std::vector<double> I_x_;										///< I_x as defined for the derivation of the helicity Upsilon
	std::vector<double> I_y_;										///< I_y as defined for the derivation of the helicity Upsilon
	std::vector<double> I_x_2_;										///< I_x^2 as defined for the derivation of the helicity Upsilon
	std::vector<double> I_y_2_;										///< I_y^2 as defined for the derivation of the helicity Upsilon
	std::vector<double> Upsilon_;									///< Helicity Upsilon

	std::vector<double> abs_vortices_;								///< Stores total vortex density
	std::vector<double> signed_vortices_;							///< Stores total signed vortex density (i.e. positive minus negative density).

	std::vector<double> coordination_number_;						///< Coordination number in mobile model (avg. number of neighbors in interaction radius)

	std::vector<topology::Vector2d> qvals_;							///< Values of q for evaluation of field fluctuations (changes over sampling process)

	std::vector<std::complex<double> > mxq_;						///< Stores the \f$m_{x,q}\f$ (x-spin field fluctuation)
	std::vector<std::complex<double> > myq_;						///< Stores the \f$m_{y,q}\f$ (y-spin field fluctuation)
	std::vector<std::complex<double> > mparq_;						///< Stores the \f$m_{\parallel,q}\f$ (spin field fluctuation parallel to spontaneous M)
	std::vector<std::complex<double> > mperpq_;						///< Stores the \f$m_{\perp,q}\f$ (spin field fluctuation perpendicular to spontaneous M)
	std::vector<std::complex<double> > eq_;							///< Stores the \f$e_{q}\f$ (energy field fluctuation)
	std::vector<std::complex<double> > wq_;							///< Stores the \f$w_{q}\f$ (omega field fluctuation)
	std::vector<std::complex<double> > teq_;						///< Stores the \f$\theta_{q}\f$ (theta field fluctuation)
	std::vector<std::complex<double> > rq_;							///< Stores the \f$\rho_{q}\f$ (density field fluctuation)
	std::vector<std::complex<double> > jparq_;						///< Stores the \f$j_{\parallel,q}\f$ (parallel momentum field fluctuation)
	std::vector<std::complex<double> > jperpq_;						///< Stores the \f$j_{\perp,q}\f$ (perpendicular momentum field fluctuation)
	std::vector<std::complex<double> > lq_;							///< Stores the \f$l_{q}\f$ (spatial angular momentum field fluctuation)

	std::vector<std::complex<double> > mxq_cur_;					///< Stores the \f$m_{x,q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > myq_cur_;					///< Stores the \f$m_{y,q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > mparq_cur_;					///< Stores the \f$m_{\parallel,q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > mperpq_cur_;					///< Stores the \f$m_{\perp,q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > eq_cur_;						///< Stores the \f$e_{q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > wq_cur_;						///< Stores the \f$w_{q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > teq_cur_;					///< Stores the \f$\theta_{q}\f$ at the current time (avoids repeated computation, which is costly)
	std::vector<std::complex<double> > rq_cur_;						///< Stores the \f$\rho_{q}\f$ at the current time
	std::vector<std::complex<double> > jparq_cur_;					///< Stores the \f$j_{\parallel,q}\f$ at the current time
	std::vector<std::complex<double> > jperpq_cur_;					///< Stores the \f$j_{\perp,q}\f$ at the current time
	std::vector<std::complex<double> > lq_cur_;						///< Stores the \f$l_{q}\f$ at the current time

	std::vector<std::complex<double> > convol_wmx_cur_;				///< Stores the convolution \f$w_{q} \star m_{x,q}\f$
	std::vector<std::complex<double> > convol_wmy_cur_;				///< Stores the convolution \f$w_{q} \star m_{y,q}\f$
	std::vector<std::complex<double> > convol_jparmx_cur_;			///< Stores the convolution \f$j_{L,q} \star m_{x,q}\f$
	std::vector<std::complex<double> > convol_jparmy_cur_;			///< Stores the convolution \f$j_{L,q} \star m_{y,q}\f$


	std::vector<std::complex<double> > mxq_initial_;				///< Stores the \f$m_{x,q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > myq_initial_;				///< Stores the \f$m_{y,q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > mparq_initial_;				///< Stores the \f$m_{\parallel,q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > mperpq_initial_;				///< Stores the \f$m_{\perp,q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > eq_initial_;					///< Stores the \f$e_{q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > wq_initial_;					///< Stores the \f$w_{q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > teq_initial_;				///< Stores the \f$\theta_{q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > rq_initial_;					///< Stores the \f$\rho_{q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > jparq_initial_;				///< Stores the \f$j_{\parallel,q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > jperpq_initial_;				///< Stores the \f$j_{\perp,q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)
	std::vector<std::complex<double> > lq_initial_;					///< Stores the \f$l_{q}\f$ at the initial time (avoids repeated computation, only usable if q is not refreshed)

	std::vector<std::complex<double> > convol_wmx_initial_;			///< Stores the convolution \f$w_{q} \star m_{x,q}\f$ at the initial time
	std::vector<std::complex<double> > convol_wmy_initial_;			///< Stores the convolution \f$w_{q} \star m_{y,q}\f$ at the initial time
	std::vector<std::complex<double> > convol_jparmx_initial_;		///< Stores the convolution \f$j_{\parallel,q} \star m_{x,q}\f$ at the initial time
	std::vector<std::complex<double> > convol_jparmy_initial_;		///< Stores the convolution \f$j_{\parallel,q} \star m_{y,q}\f$ at the initial time


	std::vector<double> chimxq_;									///< Stores the \f$\chi_{mx,q}\f$
	std::vector<double> chimyq_;									///< Stores the \f$\chi_{my,q}\f$
	std::vector<double> chimparq_;									///< Stores the \f$\chi_{m\parallel,q}\f$
	std::vector<double> chimperpq_;									///< Stores the \f$\chi_{m\perp,q}\f$
	std::vector<double> chieq_;										///< Stores the \f$\chi_{e,q}\f$
	std::vector<double> chiwq_;										///< Stores the \f$\chi_{w,q}\f$
	std::vector<double> chiteq_;									///< Stores the \f$\chi_{\theta,q}\f$
	std::vector<double> chirq_;										///< Stores the \f$\chi_{\rho,q}\f$
	std::vector<double> chijparq_;									///< Stores the \f$\chi_{jL,q}\f$
	std::vector<double> chijperpq_;									///< Stores the \f$\chi_{jT,q}\f$
	std::vector<double> chilq_;										///< Stores the \f$\chi_{l,q}\f$

	std::vector<std::complex<double> > SCFq_xy_;					///< Stores the \f$\langle m_{x,q}^* m_{y,q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_xw_;					///< Stores the \f$\langle m_{x,q}^* w_{q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_xe_;					///< Stores the \f$\langle m_{x,q}^* e_{q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_yw_;					///< Stores the \f$\langle m_{y,q}^* w_{q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_ye_;					///< Stores the \f$\langle m_{y,q}^* e_{q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_we_;					///< Stores the \f$\langle w_{q}^* e_{q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_mparmperp_;				///< Stores the \f$\langle m_{\parallel,q}^* m_{\perp,q} \rangle\f$ static correlation
	std::vector<std::complex<double> > SCFq_re_;					///< Stores the \f$\langle r_{q}^* e_{q} \rangle\f$ static correlation


	std::vector<double> SCF_Spin_;									///< Stores the static spin correlation function
	std::vector<double> SCF_Spin_par_;								///< Stores the static spin correlation function
	std::vector<double> SCF_Spin_perp_;								///< Stores the static spin correlation function
	std::vector<double> SCF_anglediff_;								///< Stores the static angle difference correlation function
	std::vector<double> SCF_W_;										///< Stores the static spin momentum correlation function
	std::vector<double> SCF_P_;										///< Stores the static linear momentum correlation function
	std::vector<double> SCF_g_;										///< Stores the pair distribution function
	std::vector<double> SCF_E_;										///< Stores the static total energy correlation function
	std::vector<double> SCF_Ekin_;									///< Stores the static kinetic energy correlation function
	std::vector<double> SCF_Eint_;									///< Stores the static interaction energy correlation function

	std::vector<double> ACF_Spin_;									///< Stores the spin autocorrelation function
	std::vector<double> ACF_anglediff_;								///< Stores the angle difference autocorrelation function


	std::vector<double> ACF_q0_M_;									///< Stores the magnetization autocorrelation function (q=0 limit)
	std::vector<double> ACF_q0_absM_;								///< Stores the autocorrelation function of the absolute magnetization (q=0 limit), basically the limit q=0 for mparq

	std::vector<double> ACF_Sx_;									///< Stores the spin autocorrelation function in x direction
	std::vector<double> ACF_Sy_;									///< Stores the spin autocorrelation function in y direction
	std::vector<double> ACF_W_;										///< Stores the omega autocorrelation function in y direction
	std::vector<double> ACF_E_;										///< Stores the energy autocorrelation function in y direction
	std::vector<double> ACF_Eint_;									///< Stores the interaction energy autocorrelation function in y direction
	std::vector<double> ACF_Ekin_;									///< Stores the kinetic energy autocorrelation function in y direction
	std::vector<double> ACF_P_;										///< Stores the momentum autocorrelation function in y direction
	std::vector<double> ACF_Ppar_;									///< Stores the momentum autocorrelation function parallel to the magnetization
	std::vector<double> ACF_Pperp_;									///< Stores the spin autocorrelation function perpendicular to the magnetization
	std::vector<double> ACF_MSD_;									///< Stores the positional mean squared displacement.
	std::vector<double> ACF_ang_MSD_;								///< Stores the angular mean squared displacement.
	std::vector<double> MSD_aux_;									///< Auxiliary MSD storage. Stores all distances and angles (depending on group type) that the particles have travelled.


	std::vector<std::complex<double> > gxx_;						///< Stores the \f$C_{xx}(q,t) = g_{x,x,q}(t) = \langle m_{x,q}^*m_{x,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gxy_;						///< Stores the \f$C_{xy}(q,t) = g_{x,y,q}(t) = \langle m_{x,q}^*m_{y,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gxw_;						///< Stores the \f$C_{xw}(q,t) = g_{x,w,q}(t) = \langle m_{x,q}^*w_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gxe_;						///< Stores the \f$C_{we}(q,t) = g_{x,e,q}(t) = \langle m_{x,q}^*e_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gyy_;						///< Stores the \f$C_{yy}(q,t) = g_{y,y,q}(t) = \langle m_{y,q}^*m_{y,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gyw_;						///< Stores the \f$C_{yw}(q,t) = g_{y,w,q}(t) = \langle m_{y,q}^*w_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gye_;						///< Stores the \f$C_{ye}(q,t) = g_{y,e,q}(t) = \langle m_{y,q}^*e_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gww_;						///< Stores the \f$C_{ww}(q,t) = g_{w,w,q}(t) = \langle w_{q}^*w_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gwe_;						///< Stores the \f$C_{we}(q,t) = g_{w,e,q}(t) = \langle w_{q}^*e_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gee_;						///< Stores the \f$C_{ee}(q,t) = g_{e,e,q}(t) = \langle e_{q}^*e_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gmparmpar_;					///< Stores the \f$C_{m\parallel,m\parallel}(q,t) = g_{m\parallel,m\parallel,q}(t) = \langle m_{\parallel,q}^*m_{\parallel,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gmperpmperp_; 				///< Stores the \f$C_{m\perp, m\perp}(q,t) = g_{m\perp,m\perp,q}(t) = \langle m_{\perp,q}^*m_{\perp,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gmparmperp_;					///< Stores the \f$C_{m\parallel, m\perp}(q,t) = g_{m\parallel,m\perp,q}(t) = \langle m_{\perp,q}^*m_{\parallel,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gre_;						///< Stores the \f$C_{\rho e}(q,t) = g_{r,e,q}(t) = \langle \rho_{q}^*e_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > grr_;						///< Stores the \f$C_{\rho \rho}(q,t) = g_{r,r,q}(t) = \langle \rho_{q}^*\rho_{q}(t)\rangle\f$ (time correlation function, intermediate scattering function)
	std::vector<std::complex<double> > gjparjpar_;					///< Stores the \f$C_{jL,jL}(q,t) = g_{j\parallel,j\parallel,q}(t) = \langle j_{L,q}^*j_{L,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gjperpjperp_; 				///< Stores the \f$C_{jT,jT}(q,t) = g_{j\perp,j\perp,q}(t) = \langle j_{T,q}^*j_{T,q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gll_;						///< Stores the \f$C_{ll}(q,t) = g_{l,l,q}(t) = \langle l_{q}^*l_{q}(t)\rangle\f$ (time correlation function)
	std::vector<std::complex<double> > gtt_;						///< Stores the \f$C_{\theta\theta}(q,t) = g_{\theta,\theta,q}(t) = \langle \theta_{q}^*\theta_{q}(t)\rangle\f$ (time correlation function)

	std::vector<double> TransCoeff_J1_;								///< Transport Coefficient: \f$\langle J(r_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_J1_cos1_;						///< Transport Coefficient: \f$\langle J(r_{ij}) cos(te_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_J1_cos2_;						///< Transport Coefficient: \f$\langle J(r_{ij}) cos^2(te_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_J1_cos2_r2_;						///< Transport Coefficient: \f$\langle J(r_{ij}) cos^2(te_{ij}) r_{ij}^2\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_J1_cos1_r2_;						///< Transport Coefficient: \f$\langle J(r_{ij}) cos(te_{ij}) r_{ij}^2\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_J1_sin1_te_;						///< Transport Coefficient: \f$\langle J(r_{ij}) sin(te_{ij}) te_{ij}\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Up_;								///< Transport Coefficient: \f$\langle {\cal U}'(r_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Up_rinv_;						///< Transport Coefficient: \f$\langle {\cal U}'(r_{ij}) r_{ij}^{-1}\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Upp_;							///< Transport Coefficient: \f$\langle {\cal U}''(r_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Up_cos1_;						///< Transport Coefficient: \f$\langle {\cal U}'(r_{ij}) cos(te_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Up_rinv_cos1_;					///< Transport Coefficient: \f$\langle {\cal U}'(r_{ij}) r_{ij}^{-1} cos(te_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Upp_cos1_;						///< Transport Coefficient: \f$\langle {\cal U}''(r_{ij}) cos(te_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_cos1_;							///< Transport Coefficient: \f$\langle cos(te_{ij})\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Up_rinv_te2_;					///< Transport Coefficient: \f$\langle {\cal U}'(r_{ij}) r_{ij}^{-1} te_{ij}^2\rangle\f$ (see old notes, not included in thesis)
	std::vector<double> TransCoeff_Upp_te2_;						///< Transport Coefficient: \f$\langle {\cal U}''(r_{ij}) te_{ij}^2\rangle\f$ (see old notes, not included in thesis)

	std::vector<double> TransCoeff_times_;							///< Stores times used in calculation of transport coefficients.

	std::vector<topology::Vector2d> tau_;							///< Stores momentum flux tau, which is the q to 0 limit of $\tau'_q$ in eq. (1.31b). Required for ZM transport Coefficients. (see old notes, not included in thesis)
	std::vector<topology::Vector2d> je_;							///< Stores energy flux je, which is the q to 0 limit of $\j'^e_q$ in eq. (1.31a). Required for ZM transport Coefficients. (see old notes, not included in thesis)
	topology::Vector2d tau_initial_;								///< Stores initial momentum flux tau (see old notes, not included in thesis)
	topology::Vector2d je_initial_;									///< Stores initial energy flux tau (see old notes, not included in thesis)
	std::vector<double> kappa_TransCoeff_;							///< Stores kappa transport coefficient, see Bissinger notes (sec. 1.2.7), (see old notes, not included in thesis)
	std::vector<double> Gamma_TransCoeff_;							///< Stores Gamma transport coefficient, see Bissinger notes (sec. 1.2.7), (see old notes, not included in thesis)
	std::vector<double> kappa_TransCoeff_new_;						///< Stores different kappa transport coefficient, see Bissinger notes (sec. 1.2.7), (see old notes, not included in thesis)
	std::vector<double> Gamma_TransCoeff_new_;						///< Stores different Gamma transport coefficient, see Bissinger notes (sec. 1.2.7), (see old notes, not included in thesis)

	std::vector<std::complex<double> > K_convol_wmx_;				///< Stores the memory kernel of the convolution \f$w_q \star m_{x,q}\f$ with itself
	std::vector<std::complex<double> > K_convol_wmy_;				///< Stores the memory kernel of the convolution \f$w_q \star m_{y,q}\f$ with itself
	std::vector<std::complex<double> > K_convol_jmx_;				///< Stores the memory kernel of the convolution \f$j_{L,q} \star m_{x,q}\f$ with itself
	std::vector<std::complex<double> > K_convol_jmy_;				///< Stores the memory kernel of the convolution \f$j_{L,q} \star m_{y,q}\f$ with itself
	std::vector<std::complex<double> > K_convol_wmx_jmy_;			///< Stores the memory kernel of the convolution \f$w_q \star m_{x,q}\f$ with \f$j_{L,q} \star m_{y,q}\f$
	std::vector<std::complex<double> > K_convol_wmy_jmx_;			///< Stores the memory kernel of the convolution \f$w_q \star mxy\f$ with \f$j_{L,q} \star m_{x,q}\f$
	std::vector<std::complex<double> > K_wmx_;						///< Stores the memory kernel of the quadratic form \f$w_q m_{x,q}\f$ with itself
	std::vector<std::complex<double> > K_wmy_;						///< Stores the memory kernel of the quadratic form \f$w_q m_{y,q}\f$ with itself
	std::vector<std::complex<double> > K_jmx_;						///< Stores the memory kernel of the quadratic form \f$jqpar m_{x,q}\f$ with itself
	std::vector<std::complex<double> > K_jmy_;						///< Stores the memory kernel of the quadratic form \f$jqpar m_{y,q}\f$ with itself
	std::vector<std::complex<double> > K_wmx_jmy_;					///< Stores the memory kernel of the quadratic form \f$w_q \star m_{x,q}\f$ with the quadratic form \f$j_{L,q} m_{y,q}\f$
	std::vector<std::complex<double> > K_wmy_jmx_;					///< Stores the memory kernel of the quadratic form \f$w_q \star mxy\f$ with the quadratic form \f$j_{L,q} m_{x,q}\f$


	bool store_static_ = 0;											///< Decides whether or not static quantities are calculated.

	bool store_vortices_ = 0;										///< Decides whether or not vortices_ is calculated.

	bool store_mxq_ = 0;											///< Decides whether or not \f$m_{x,q}\f$ is calculated.
	bool store_myq_ = 0;											///< Decides whether or not \f$m_{y,q}\f$ is calculated.
	bool store_eq_ = 0;												///< Decides whether or not \f$e_q\f$ is calculated.
	bool store_wq_ = 0;												///< Decides whether or not \f$w_q\f$ is calculated.
	bool store_rq_ = 0;												///< Decides whether or not \f$\rho_q\f$ is calculated.
	bool store_jparq_ = 0;											///< Decides whether or not \f$j_{L,q}\f$ is calculated.
	bool store_jperpq_ = 0;											///< Decides whether or not \f$j_{T,q}\f$ is calculated.
	bool store_lq_ = 0;												///< Decides whether or not \f$l_q\f$ is calculated.

	bool store_teq_ = 0;											///< Decides whether or not \f$\theta_q\f$ is calculated.

	bool store_SCF_ = 0;											///< Decides whether or not the static correlation functions are calculated.
	bool store_TCF_ = 0;											///< Decides whether or not the time correlation functions (autocorrelations and so forth) are calculated.
	bool refresh_q_ = 0;											///< Decides whether or not to choose new random q values for each sampling.

	bool store_TransCoeff_ = 0;										///< Decides whether or not the transport coefficients are calculated.
	bool store_MemoryKernels_ = 0;									///< Decides whether or not the memory kernels (MCT) are calculated.

	bool print_snapshots_ = 0;										///< Decides whether trajectories should be sampled.

public:
	//   ===================================================================================================
	//   Constructors
	//   ===================================================================================================
	/// Empty constructor
	sampler(){};
	/// Constructor. Assigns par_ to par.
	sampler(parameters par) : par_(par) {} ;
	/// Constructor. Assigns par_ to par and qvals_ to qvals.
	sampler(parameters par, std::vector<topology::Vector2d> qvals) : par_(par), qvals_(qvals) {} ;


	//   ===================================================================================================
	//   Intialization
	//   ===================================================================================================
	/// Sets parameter member variable
	void set_parameters(parameters par);
	/// Sets qvals (redundant)
	void set_qvals(std::vector<topology::Vector2d> qvals);

	//   ===================================================================================================
	//   Switches management
	//   ===================================================================================================

	/// Reads switches from vector
	void switches_from_vector(const std::vector<bool>& boolvec);
	/// Reads switches from file
	void switches_from_file(std::ifstream& infile, std::ofstream& stdoutfile);

	/// Turns all switches on
	void all_switches_on();
	/// Turns all switches off
	void all_switches_off();

	/// Turns print_snaphshots_ switch on
	void print_snapshots_on();
	/// Turns print_snaphshots_ switch off
	void print_snapshots_off();
	/// Turns all switches off if on_fly_sampling is off
	void check_on_fly_sampling();

	//   ===================================================================================================
	//   get functions
	//   ===================================================================================================
	std::vector<bool> get_switches() const; ///< Returns vector of sample switches. Is for checking.
	inline bool print_snapshots() const { return print_snapshots_; };
	std::vector<topology::Vector2d> get_qvals() const { return qvals_; };

	//   ===================================================================================================
	//   binning functions and qvalue handling
	//   ===================================================================================================
	/// Chooses new values for q (not recommended, relic of old implementation of qbin)
	void refresh_qvals(const topology::Vector2d& boxsize);
	/// Sorts the qvals at which correlations were computed into a bin (double input)
	std::vector<double> bin_qvals_to_q(std::vector<double> vals) const ;
	/// Sorts the qvals at which correlations were computed into a bin (complex input)
	std::vector<std::complex<double> > bin_qvals_to_q(std::vector<std::complex<double> > vals) const ;


	//   ===================================================================================================
	//   Data access (only for very useful data, shouldn't be overdone)
	//   ===================================================================================================
	/// Returns last temperature that was calculated. Careful, no check if calcultaion was performed.
	inline double get_last_temperature() { return temperature_.back(); };

	//   ===================================================================================================
	//   sampling functions
	//   ===================================================================================================
	/// Samples MSD.
	void sample_MSD(const std::vector<double>& MSD);

	/// Samples all static properties and dynamic properties according to switches.
	void sample(const group& G, const group& G_new, double t);
	/// Samples all static properties that have an activated sample switch.
	void sample_static(const group& G, double t);
	 /// Samples all time correlation functions that have an activated sample switch.
	void sample_TCF(const group& G_initial, const group& G_new, const double t);

	/// Samples (i.e. stores) a snapshot to a file.
	void sample_snapshots(const group& G, const double t, std::string snapfile_name, std::ofstream& outfile);

	/// appends averages to the last entry of each nonempty sampling vector
	void average();

	//   ===================================================================================================
	//   printing functions
	//   ===================================================================================================
	/// Prints in a format that can be directly read in from matlab.
	void print_matlab(std::ofstream& outfile);
	/// Prints averages in a format that can be directly read in from matlab.
	void print_averages_matlab(std::ofstream& outfile);


};
#endif /* SAMPLER_H_ */
