/*! \file inputoutput.h
 *  \brief Provides routines for printing std::vectors.
 *
 *  Used to be a bigger suite of programs that would also handle data input.
 *  The task of reading input has been shifted to the parameter class that
 *  handles the entirety of input. Now it's a small set of routines
 *  suited for adjustable printing conditions.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: mid-2017
 *  \date Last Update: 23-08-02
 *
 */
#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include"topology.h"

#include "math.h"
#include<stdio.h>
#include<stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>


/*! \brief Prints a std::vector (doubles) to an outflow-stream.
 *
 *  The result is comma-separated.
 */
void print_stdvector(std::vector<double> vec, std::ostream& outfile);

/*! \brief Prints a std::vector (doubles) to an outflow-stream. Includes prefix and suffix.
 *
 *  The result is comma-separated. The startstring and endstring can be used to specify a format
 *  of the vector.
 */
void print_stdvector(std::vector<double> vec, std::ostream& outfile, std::string startstring, std::string endstring);

/*! \brief Prints a std::vector (doubles) to an outflow-stream. Includes prefix and suffix.
 *  The entries to be printed can be chosen manually.
 *
 *  Does not put a line-break at the end (if desired, include in endstring).
 *
 *  Does not check for errors -- e.g. a negative index_start will not be caught.
 *
 *  The result is comma-separated. The startstring and endstring can be used to specify a format
 *  of the vector. Output for the vector \f$v = (v_1,\ldots,v_N)\f$ is of the form:
 *
 *  "startstring" + "\f$v_{\texttt{index_start}}, v_{\texttt{index_start} + \texttt{stepskip},
 *  \ldots,  v_{\texttt{index_end}}\f$" + "endstring"
 *
 *  if \f$\texttt{index_end} \neq \texttt{index_estart}+ n \cdot\texttt{stepskip}\f$  for some \f$n\f$, the
 *  next-lowest value is chosen with an \f$n\f$ such that
 *  \f$\texttt{index_estart}+ n \cdot\texttt{stepskip} < \texttt{index_end} < \texttt{index_estart}+ (n+1) \cdot\texttt{stepskip}\f$.
 *
 *  \param[in] vec The vector to be printed.
 *
 *  \param[in] outfile The stream to the output file
 *
 *  \param[in] startstring String placed before the vector entries are printed.
 *
 *  \param[in] endstring String placed behind the printed vector entries.
 *
 *  \param[in] index_start First index to be printed.
 *
 *  \param[in] index_end Last index to be printed.
 *
 *  \param[in] stepskip Separation between entries to be printed.
 */
void print_stdvector(std::vector<double> vec, std::ostream& outfile, std::string startstring, std::string endstring,
		int index_start, int index_end, int stepskip);

/*! \brief Prints a std::vector (complex numbers) to an outflow-stream. Includes prefix and suffix.
 *  The entries to be printed can be chosen manually.
 *
 *  Does not put a line-break at the end (if desired, include in endstring).
 *
 *  Does not check for errors -- e.g. a negative index_start will not be caught.
 *
 *  The result is comma-separated. The startstring and endstring can be used to specify a format
 *  of the vector. Output for the vector \f$v = (v_1,\ldots,v_N)\f$ is of the form:
 *
 *  "startstring" + "\f$v_{\texttt{index_start}}, v_{\texttt{index_start} + \texttt{stepskip},
 *  \ldots,  v_{\texttt{index_end}}\f$" + "endstring"
 *
 *  if \f$\texttt{index_end} \neq \texttt{index_estart}+ n \cdot\texttt{stepskip}\f$  for some \f$n\f$, the
 *  next-lowest value is chosen with an \f$n\f$ such that
 *  \f$\texttt{index_estart}+ n \cdot\texttt{stepskip} < \texttt{index_end} < \texttt{index_estart}+ (n+1) \cdot\texttt{stepskip}\f$.
 *
 *  \param[in] vec The vector to be printed.
 *
 *  \param[in] outfile The stream to the output file
 *
 *  \param[in] startstring String placed before the vector entries are printed.
 *
 *  \param[in] endstring String placed behind the printed vector entries.
 *
 *  \param[in] index_start First index to be printed.
 *
 *  \param[in] index_end Last index to be printed.
 *
 *  \param[in] stepskip Separation between entries to be printed.
 */
void print_stdvector(std::vector<std::complex<double> > vec, std::ostream& outfile, std::string startstring, std::string endstring,
		int index_start, int index_end, int stepskip);


/*! \brief Prints a std::vector (integers) to an outflow-stream. Includes prefix and suffix.
 *
 *  The result is comma-separated. The startstring and endstring can be used to specify a format
 *  of the vector.
 */
void print_stdvector(std::vector<int> vec, std::ostream& outfile, std::string startstring, std::string endstring);

#endif
