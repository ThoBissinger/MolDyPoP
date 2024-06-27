/*! \file routines.h
 *
 *  \brief Header-File to the namespace routines. The namespace contains simulation routines
 *  that manage setting up the simulation, running it and communicating the results.
 *
 *  These routines are called by the file main.cpp and are the crucial building blocks of the simulation.
 (
 *  \date Created: 2019-12-11
 *  \date Last Updated: 2023-08-06
 *
 *  \namespace  routines
 *  \brief Different calculation routines with xygroups and mxygroups. Carries out simulation tasks.
 *
 *  \author Thomas Bissinger
 *  \date Created: 2019-12-11
 *  \date Last Updated: 2023-08-06

 *
 */
#include <vector>
#include "sampler.h"
#include "inputoutput.h"
#include "integrator.h"
#include "computations.h"
#include "topology.h"
#include <iostream>
#include <fstream>

#include <chrono>

#ifndef ROUTINES_H_
#define ROUTINES_H_


namespace routines {

	//   ===================================================================================================
	//   General routines
	//   ===================================================================================================
	/*! \brief basic integration routine
	 *
	 *  Proceeds as follows:
	 *
	 *  1. Preliminary stuff (opening files, initial print)
	 *
	 *  2. Initializes all relevant objects (group, integrator, sampler)
	 *
	 *  3. Performs an equilibration run (typically with check to temperature, depends on switch)
	 *
	 *  4. Performs an integration run (used for sampling)
	 *
	 *  5. Cleanup, final prints
	 *
	 *  \note The routine equilibrate and the routine integrate_snapshots are used within this routine.
	 *
	 *  \note Equilibration is not really managed in an elegant way. It is recommended to check
	 *  manually whether data has been equilbrated and to use a fixed equilibration time.
	 */
	int integration(parameters par);

	/*! \brief basic sampling routine (no integration performed)
	 *
	 *  Reads data from snapshot files provided in a list file
	 *  and performs sampling on it. Proceeds as follows:
	 *
	 *  1. Preliminary stuff (opening files, initial print)
	 *
	 *  2. Initializes all relevant objects (group and sampler)
	 *
	 *  3. For each time step in the list file, reads out the group and performs
	 *  sampling on the group at that time instant.
	 *
	 *  4. Cleanup, final prints
	 *
	 *  \note This routine can only be used when snapshots are stored during another
	 *  integration/sampling run. Useful for explorative investigations, but large numbers
	 *  of snapshots should not be stored for a large sample and it is recommended
	 *  to perform on-fly sampling.
	 *
	 */
	int sampling(parameters par);


	//   ===================================================================================================
	//   General routines
	//   ===================================================================================================
	/*! \brief equilibration routine
	 *
	 *  Takes in a group and integrates it until Tmax (or another break condition is
	 *  met), then returns whether or not the group is equilibrated then.
	 *
	 *  Reads data from snapshot files provided in a list file
	 *  and performs sampling on it. Proceeds as follows:
	 *
	 *  1. Preliminary stuff (opening files, initial print)
	 *
	 *  2. Performs time integration (depending on which integrator chosen)
	 *
	 *  3. After a time set in par, the integration checks for equlibration
	 *  and decides whether or not to continue equlibrating
	 *
	 *  4. Cleanup, final prints
	 *
	 *
	 *
	 *  @param[inout] G The group that should be equilibrated
	 *  @param[in] par A set of simulation parameters
	 *  @param[inout] samp The sampler in which equilibration data should be stored (careful,
	 *  will be reset during run - TODO!)
	 *  @param[in] Tmax Maximum equilibration time
	 *  @param[in] t Current time
	 *  @param[in] breakcond Break condition. The following values can be taken:
	 *  <table>
	 *    <caption id="multi_row">Values of fluctname</caption>
	 *      <tr><th> breakcond value 	<th>meaning
	 *      <tr><td> "time"				<td>Wait until Tmax. After that, another equilibration check
	 *                                      is performed and the routine may be called again.
	 *      <tr><td> "time_hard"		<td>Wait until Tmax. No further equilibration performed.
	 *      <tr><td> "temperature"      <td>Breaks if the desired temperature is reached and maintained
	 *                                      for a certain amount of time.
	 *      <tr><td> "any" 				<td>Any of the above.
	 *   </table>
	 *   @param[in] stdoutfile File to which output is to be printed.
	 *
	 *
	 */
	int equilibrate(group& G, const parameters& par, sampler& samp,
				const double Tmax, double& t,
				const std::string breakcond, std::ofstream& stdoutfile);


	//   ===================================================================================================
	//   Routines for sampling
	//   ===================================================================================================

	/*! \brief integration routine (the one that does the work)
	 *
	 *  Takes in a group and integrates it until par.Tmax(). May store snapshots or perform
	 *  sampling on the fly, depending on parameters. Many details depend on the parameters
	 *  chosen and can be checked in the declaration of parameters.h
	 *
	 *  Proceeds as follows:
	 *
	 *  1. Preliminary stuff (opening files, initial print)
	 *
	 *  2. Initializes group (setting positions to lattice, initializing partition,
	 *  drawing random velocities etc.)
	 *
	 *  3. Performs time integration (depending on which integrator chosen)
	 *
	 *  4. During integration, samples and stores data
	 *
	 *  5. Cleanup, final prints - no print of sampled data, that is done
	 *  in the routine integration that typically calls for this function
	 *
	 *
	 *
	 *  @param[inout] G The group that should be equilibrated
	 *  @param[in] par A set of simulation parameters
	 *  @param[inout] samp The sampler in which the run data will be stored
	 *  @param[in] stdoutfile File to which output is to be printed.
	 *
	 *
	 */
	void integrate_snapshots(group& G, const parameters& par, sampler& samp,
				std::ofstream& stdoutfile);

	//   ===================================================================================================
	//   Internal routines (input/output etc)
	//   ===================================================================================================
	/// Initial print of each routine. States the routine name.
	void initprint(std::string routine_name, std::ofstream& outfile);
	/// Terminal print of each routine. States the routine name.
	void terminateprint(std::string routine_name, std::ofstream& outfile);
};
#endif /* ROUTINES_H_ */
