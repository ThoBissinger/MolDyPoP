/*! \file main.cpp
 *
 *  \brief Main-file. Every computation starts here.
 *
 *  The main-File progresses as follows:
 *
 *  1. Read the file at input/infile.in (called infile for short)
 *
 *  2. Initialize parameter class object par
 *
 *  3. Read the infile into par. In case of input errors, abort the run.
 *     Return value is the error code thrown by the read function.
 *
 *  4. Depending on the value of the input parameter mode, choose
 *
 *  	- mode = none. Do nothing.
 *
 *  	- mode = test. User-specified tests can be performed.
 *
 *  	- mode = integrate, mode=integrate_cont, mode=equilibrate. Run specific
 *  	  integration routines from the routines namespace.
 *
 *  	- mode = samp. Run specific sampling routines from the routines namespace.
 *
 *  5. End the run, return value is the return value of the routine performed.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: Mid-February 2017
 *  \date Last updated: 23-08-01
 */

#include"computations.h"
#include"partition.h"
#include"inputoutput.h"
#include"topology.h"
#include"routines.h"
#include"integrator.h"
#include"parameters.h"
#include"sampler.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "math.h"
#include <chrono>  // for high_resolution_clock
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <ctime>
#include <complex>

#include <chrono>


// ===========================================================================================
// BEGINNING MAIN
// ===========================================================================================
int main() {
	std::ifstream infile("input/infile.in", std::ios::in);

	parameters par;
	int errorcode = par.read_from_file(infile);
	if (errorcode == 1){
		std::cout << "Input errors. Abort computation." << std::endl;
		return errorcode;
	}
	int returnval = -1;
	if (par.mode() == "none"){
		std::cout << "Simulation mode: " << par.mode() << ". No output generated" << std::endl;
		returnval = 0;
	} else if (par.mode() == "test") {
		std::cout << "Simulation mode: " << par.mode() << ". Testing computation" << std::endl;
		returnval = 0;

		/*  In this space, you can insert test code for direct testing.
		 *
		 */
	} else {
		if ((par.mode() == "integrate") | (par.mode() == "integrate_cont") | (par.mode() == "equilibrate")){
			returnval = routines::integration(par);
		}
		if (par.mode() == "samp"){
			returnval = routines::sampling(par);
		}
	}


	std::cout << "Computation finished. \nExit value: " << returnval << std::endl;
	return returnval;

}
