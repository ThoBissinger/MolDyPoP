/*! \file computations.h
 *  \brief Contains various computation methods that do not belong to a particular class.
 *
 *  Main functionalities: Provides the modulus function, provides code for
 *  averaging over data stored in a std::vector with some selection rules,
 *  handles treating std::vectors as representing matrices and performing
 *  column and row averages. Also has a few methods for computing random
 *  variables.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: mid-2017
 *  \date Last Update: 23-08-02
 *
 */


#ifndef COMPUTATIONS_H
#define COMPUTATIONS_H
#include<stdlib.h>
#include <vector>
#include <complex>
#include<math.h>
#include <numeric>
#include <algorithm>
#include "topology.h"


/*! \brief Manual introduction of M_PI, signifying \f$\pi\f$ up to computer precision.
 *
 *  Some (mostly windows) compilers have problems with the definition of M_PI in math.h.
 *  For them, it is manually introduced here.
 */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif




// =======================================================================================================
// general functions
// =======================================================================================================
/// Modulus function between min and max.
double mod(const double val, const double min, const double max);

// =======================================================================================================
// functions for std::vector
// =======================================================================================================
/// Calculates mean of vector. Template, instantiated with double, int and std::complex<double> in computations.cpp
template<typename T, typename A>
T vector_mean(std::vector<T,A> v);

/// Calculates mean of row of matrix (vector is M matrix rows of length N together). Template, instantiated with double, int and std::complex<double> in computations.cpp
template<typename T, typename A>
T row_mean(std::vector<T,A> v, int M, int N, int i);

/// Calculates mean of column of matrix (vector is M matrix rows of length N together). Template, instantiated with double, int and std::complex<double> in computations.cpp
template<typename T, typename A>
T column_mean(std::vector<T,A> v, int M, int N, int j);

/// Calculates mean of vector, but only for entries that are at least xsep_min apart.
template<typename T, typename A>
T selective_vector_mean(std::vector<T,A> yvals, std::vector<double> xvals, double xsep_min);

/// Calculates row mean of vector, but only for entries that are at least xsep_min apart.
template<typename T, typename A>
T selective_row_mean(std::vector<T,A> yvals, std::vector<double> xvals, double xsep_min, int M, int N, int i);

template<typename T, typename A>
T selective_column_mean(std::vector<T,A> yvals, std::vector<double> xvals, double xsep_min, int M, int N, int j); ///< Calculates column mean of vector, but only for entries that are at least xsep_min apart.

/// Calculates variance of vector, also stores mean.
/*! Does not make much sense to use templates here, since
 *  the variance of an int vector is not much use and
 *  the variance of std::complex would require complex conjugation.
 */
double vector_variance(std::vector<double> v, double& mean);

/// Scales vector by factor a.
template<typename T, typename A>
std::vector<T,A> vector_scale(const std::vector<T,A>& v, double a);

/// Adds two vectors v and w.
template<typename T, typename A>
std::vector<T,A> vector_sum(const std::vector<T,A>& v, const std::vector<T,A>& w);

/// Inner product of two vectors v and w.
template<typename T, typename A>
T vector_inpr(const std::vector<T,A>& v, const std::vector<T,A>& w);

/// Inner product of two complex vectors v and w, where the complex conjugate of the components of v is used.
std::complex<double> vector_inpr(const std::vector<std::complex<double> >& v, const std::vector<std::complex<double> >& w);

/// Pointwise multiplication of two vectors v and w.
template<typename T, typename A>
std::vector<T,A> vector_pw_mult(const std::vector<T,A>& v, const std::vector<T,A>& w);

/// Pointwise multiplication of two complex vectors v and w, where the complex conjugate of the components of v is used.
std::vector<std::complex<double> > vector_pw_mult(const std::vector<std::complex<double> >& v, const std::vector<std::complex<double> >& w);

/// Pointwise multiplication of four complex vectors v and w, where the complex conjugate of the components of v_1 and v_2 is used.
std::vector<std::complex<double> > vector_pw_mult(const std::vector<std::complex<double> >& v_1, const std::vector<std::complex<double> >& v_2,
		const std::vector<std::complex<double> >& v_3, const std::vector<std::complex<double> >& v_4);

/// Pointwise norm value of a complex vector v, i.e. sum v[i]^* v[i].
std::vector<double > vector_pw_norm(const std::vector<std::complex<double> >& v);

/*! \brief Gives vector index if a vector of size \f$M \cdot N\f$ stands for an \f$M \times N\f$ matrix.
 *
 *  In the current implementation, a matrix \f$\mathbf{A} \in \mathcal{R}^{M\times N}\f$ is related
 *  to a vector \f$\mathbf{a} \in \mathcal{R}^{M\cdot N}\f$ by identifying
 *  \f[
 *      a_k = A_{ij}
 *  \f]
 *  with the prescription \f$k(i,j) = i + N\cdot k\f$.
 *
 *  Therefore, the value of M is not important. For completeness, it is taken as input
 *  (if anyone wishes to change the mapping between \f$k\f$ and \f$i,j\f$.
 *
 *  Careful, numbering starts at 0, i.e. last index is actually \f$i = M-1\f$, \f$j = N-1\f$.
 */
inline int matr_index(int M, int N, int i, int j) { return i * N + j; };

// =======================================================================================================
// functions for binning
// =======================================================================================================
/*! \fn log_bin()
 *
 *  \brief Returns a logarithmic bin.
 *
 *  \param binmin	Minimum value of bin (i.e. starting value of bin) [double]
 *  \param binmax	Maximum value of bin (will not be exceeded, may not be reached) [double]
 *  \param diffmin	Minimal seperation between two bins [double]
 *  \param Nbins		number of bins [int]
 *  \return			Vector of bins (can be interpreted as left boundaries of bins. [std::vector<double>]
 */
std::vector<double> log_bin(double binmin, double binmax, double diffmin, int Nbins); ///< Makes logarithmic bin.




//std::vector<double> lin_log_bin(double qmin, double qmax, double linsep, int Nlin, int Nlog); ///< Makes linear-logarithmic bin.
// /*!<
// * \fn std::vector<double> lin_log_bin()
// * \param qmin		Minimum q (starting value) [double]
// * \param qmax		Maximum q (will not be exceeded, may not be reached) [double]
// * \param linsep	Separation for the linear part of the bin [double]
// * \param Nlin		number of linear points [int]
// * \param Nlog		number of logarithmic points. Due to rounding errors, the logarithmic length can be one shorter than Nlog. [int]
// * \return			Vector of bins (can be interpreted as left boundaries of bins or mid points). [std::vector<double>]
// */

/// Returns bin index for value within bin. If value < bin[0], return 0, otherwise returns
/// i for which bin[i - 1] < value < bin[i]
int find_bin(const std::vector<double>& bin, const double& value); ///< Finds index of bin for data to be put in.
/*!<
 * Any entry that is smaller than the first entry of bin enters there,
 * then the resulting index is the index i for which bin[i - 1] < value < bin[i].
 * Note that the bin has to have one extra element for this to work,
 * see also comments for lin_log_bin().
 * \fn std::vector<double> lin_log_bin()
 */

// =======================================================================================================
// =======================================================================================================
// Random functions
// =======================================================================================================
// =======================================================================================================
/// Random int between min and max.
inline int random_int(int min, int max) { return rand() % (max - min) + min; };


/// Random double between min and max. min < max is not required.
inline double random_double(double min, double max) { return min + (double)rand() / RAND_MAX * (max - min); };
/// Random double between 0 and max.
inline double random_double(double max) { return (double)rand() / RAND_MAX * (max); };

/// Return random boltzmann-distributed double (Box-Muller method, commented out is a version using the polar method
inline double random_boltzmann_double(double kbT){	return sqrt(kbT)*sqrt(-2*log(random_double(1)))*cos(2*M_PI*random_double(1)); };

/// functions for random angle between -pi and pi.
inline double random_angle() { return random_double(-M_PI, M_PI); };
/// functions for random angle, designed such that only -pi:pi is allowed (automatical limits when values are not allowed)
double random_angle(double min, double max);
/// Random angle between -max and max (or -max and 0 if max is negative.)
double random_angle(double max);

#endif
