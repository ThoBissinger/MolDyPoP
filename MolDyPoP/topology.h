#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include<stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include<vector>
#include<array>
#include"math.h"

/*! \file topology.h
 *
 *  \brief Header-File to declaration of namespace topology. Defines
 *  classes Vector2d and angle2d (the latter redundant).
 *
 *  Provides features for handling spatial vectors both for position
 *  and spin vectors.
 *
 *
 *  \namespace topology.
 *  \brief Contains the vector classes and vector functions and operations. So far, only
 *  the class Vector2d is used in further routines and fully implemented.
 *
 *  Generalizations are possible. The namespace was originallz intended to be extensible
 *  to 3d vectors and spherical coordinates. Right now, however, no such routines are
 *  implemented. The class Angle2d exists with basic routines and some conversion
 *  and other relations with Vector2d are available, but it is not fully developed.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: early 2017
 *  \date Last Updated: 2023-08-01
 *
 */
namespace topology
{

/*! \class Vector2d.
 *  \brief Mathematical 2d vectors. Can be added, multiplied by a scalar, norm computation is possible.
 *  There are print-to-file and print-to-command-line functions available.
 *
 *  The vector can be accessed by the standard "VEC[]" command. Some additional functions are included
 *  to handle periodic boundary conditions in a box for position vectors, angles with respect to the
 *  x-axis, inner products, rotations and the like.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: early 2017
 *  \date Last Updated: 2023-08-01
 */
class Vector2d{
protected:
	std::array<double,2> v_;	///< The two components of the 2d vector.

public:
	//   ===================================================================================================
	//   Constructors
	//   ===================================================================================================
	Vector2d(); ///< Constructor without argument.
	Vector2d(double* v); ///< Constructor from double array.
	Vector2d(const double& x); ///< Constructor from a double (same arguments).
	Vector2d(const int& x); ///< Constructor from an int (same arguments).
	Vector2d(const double x, const double y); ///< Constructor from two doubles, x and y argument.

	Vector2d(const Vector2d& w); ///< Copy constructor from two doubles.

	//   ===================================================================================================
	//   Print functions
	//   ===================================================================================================
	/// Print to file in xyz format.
	void print(std::ostream &outfile) const;
	/// Print to command line in xy format.
	void print() const;
	/// Print to command line in xyz format.
	void print_xyz() const;
	/// Return x component.
	inline double get_x() const { return v_[0]; };
	/// Return y component.
	inline double get_y() const { return v_[1]; };

	void set_x(double x); ///< Assign value x to component x
	void set_y(double y); ///< Assign value y to component y

	/// Returns size of the vector (i.e. 2).
	inline const int size() const { return 2; }
	/// Returns 1 if vector is zero-vector
	inline bool is_zero() const { return ( v_[0] == 0 && v_[1] == 0 ); };

	//   ===================================================================================================
	//   Pointers to element
	//   ===================================================================================================
	double& operator[](int index); ///< Element access - only recommended in sepcific situation, use get_x(), get_y() whenever possible.
	const double& operator[](int index) const; ///< const element access - only recommended in sepcific situation, use get_x(), get_y() whenever possible.

	//   ===================================================================================================
	//   Arithmetic operations
	//   ===================================================================================================
	Vector2d& operator+=(const Vector2d& w);	///< Vector-Vector addition
	Vector2d& operator-=(const Vector2d& w);	///< Vector-Vector subtraction
	Vector2d& operator*=(const double a);		///< Vector-scalar multiplication (double)
	Vector2d& operator*=(const int a);			///< Vector-scalar multiplication (int)
	// double operator*=(const Vector2d& w);	///< Scalar product (not recommended for use)
	Vector2d& operator/=(const double a);		///< Vector-scalar division (double)
	Vector2d& operator/=(const int a);			///< Vector-scalar division (int)
	Vector2d operator-() const;					///< Unary additive inversion

	//   ===================================================================================================
	//   Rotation, norms, angles, normalization
	//   ===================================================================================================
	/// Rotates vector by angle theta.
	void rotate(const double theta);

	/// Returns squared L2 norm of vector
	inline double norm2() const { return v_[0] * v_[0] + v_[1] * v_[1]; }

	/// Returns the orientation angle of a vector.
	inline double angle() { return std::atan2(this->v_[1],this->v_[0]); };

	/// Normalizes vector.
	void normalized();

	//   ===================================================================================================
	//   Boundary conditions handling
	//   ===================================================================================================

	/// Periodic boundary conditions in box.
	/*! PERIODIC_BOX: Changes the vector according to periodic boundary conditions in a box. Can either take
	 *  two Vector2d arguments, indicating the minima and maxima of the box. If only one argument is provided,
	 *  the minima are set to zero.
	 */
	void periodic_box(const Vector2d& minima, const Vector2d& maxima);
	/// Periodic boundary conditions in box from 0 to maxima.
	void periodic_box(const Vector2d& maxima);

	/// Boundary handler within [0,maxima]. For more details, see the namesake function with arbitrary boundaries.
	Vector2d get_boundary_handler_periodic_box(const Vector2d& maxima, const double cutoff) const;

	/// Boundary handler within [minima, maxima].
	/*! A boundary handler is a Vector2d that determines whether the
	 *  vector is close to the boundary of a periodic box. If it is within the cutoff radius of the minimum
	 *  (or 0) in one component, this component will read 0. If it is within the cutoff radius to the maximum,
	 *  the component will read 1. If neither is the case, i.e. the vector does is well within the volume, the
	 *  component reads 0.
	 */
	Vector2d get_boundary_handler_periodic_box(const Vector2d& minima,const Vector2d& maxima, const double cutoff) const;

};

//   ===================================================================================================
//   Typical linear and algebraic operations.
//   ===================================================================================================
Vector2d operator+(const Vector2d& v, const Vector2d& w); ///< + operator
Vector2d operator-(const Vector2d& v, const Vector2d& w); ///< - operator
Vector2d operator*(const Vector2d& v, const double a); ///< scalar * operator
Vector2d operator*(const double a, const Vector2d& v); ///< scalar * operator
Vector2d operator*(const Vector2d& v, const int a); ///< scalar * operator
Vector2d operator*(const int a, const Vector2d& v); ///< scalar * operator
Vector2d operator/(const Vector2d& v, const double a); ///< scalar / operator
Vector2d operator/(const Vector2d& v, const int a); ///< scalar / operator


//   ===================================================================================================
//   Norms, rotations, distances
//   ===================================================================================================
/// Returns the L2 norm of a vector
inline double norm2(Vector2d v) {
  return v.norm2();
}

/// Normalizes a vector.
inline Vector2d normalized(Vector2d v) {
  return v/sqrt(v.norm2());
}

/// Returns squared distance |w - v|^2 considering periodic boundaries (square box, length L). Squared function faster to calculate.
double periodic_distance_squared(const Vector2d& v, const Vector2d& w, const double& L);
/// Returns distance |w - v| considering periodic boundaries (square box, length L). Taking sqrt takes more time than returning the squared quantity by dist_periodic_squared.
inline double periodic_distance(const Vector2d& v, const Vector2d& w, const double& L) { return std::sqrt(periodic_distance_squared(v,w,L)); };
/// Returns distance vector w - v considering periodic boundaries (square box, length L).
Vector2d periodic_distance_vector(const Vector2d& v, const Vector2d& w, const double& L);

/// Returns squared distance |w - v|^2 considering periodic boundaries (rectangular box, widths stored in L). Squared function faster to calculate.
double periodic_distance_squared(const Vector2d& v, const Vector2d& w, const Vector2d& L);
/// Returns squared distance |w - v| considering periodic boundaries (rectangular box, widths stored in L). Taking sqrt takes more time than returning the squared quantity by dist_periodic_squared.
inline double periodic_distance(const Vector2d& v, const Vector2d& w, const Vector2d& L) { return std::sqrt(periodic_distance_squared(v,w,L)); };
/// Returns distance vector w - v considering periodic boundaries (rectangular box, widths stored in L).
Vector2d periodic_distance_vector(const Vector2d& v, const Vector2d& w, const Vector2d& L);

/// Rotates a vector
inline Vector2d rotate(Vector2d v, double theta) {
  return Vector2d(cos(theta)*v.get_x()-sin(theta)*v.get_y(),
		  sin(theta)*v.get_x()+cos(theta)*v.get_y());
}

/// Rotates 90 degrees
inline Vector2d rotate_orthogonal(Vector2d v) {
	return Vector2d(-v.get_y(),v.get_x());
}

/// Inner product
inline double innerproduct(Vector2d v, Vector2d w) {
  return v.get_x() * w.get_x() + v.get_y() * w.get_y();
}

/// Parallel projection
inline double parallel_projection(Vector2d v, Vector2d w) {
	return innerproduct(v,normalized(w));
}

/// Parallel projection
inline double orthogonal_projection(Vector2d v, Vector2d w) {
	return parallel_projection(v,rotate_orthogonal(w));
}
/// Returns a random vector within a volume [minima, maxima].
Vector2d random_vector(const Vector2d& minima, const Vector2d& maxima);
/// Returns a random vector within a volume [0,maxima].
Vector2d random_vector(const Vector2d& maxima);
/// Returns a random vector with Gaussian distribution in each component.
Vector2d random_gaussian_vector(const double& sigma_squared);


/// Returns nearest neighbor to 2D-Vector v on grid of grid point separation gridsep
Vector2d nearest_gridvec(Vector2d v, double gridsep);
/// Returns nearest neighbor to 2D-Vector v on grid of anisotropic grid point separation gridsep
Vector2d nearest_gridvec(Vector2d v, Vector2d gridseps);
/// Creates a std::vector containing qsamps_per_bin 2D-vectors that lie on a grid with modulus between qmin and qmax. No vector appears double.
std::vector<Vector2d> qvalues_within_radius(double qmin, double qmax, Vector2d gridseps, int qsamps_per_bin);

/// Returns a random vector in the first quadrant.
Vector2d random_vector_first_quadrant(const double length);

/// Returns a random vector within a sector given by thetamin, thetamax.
Vector2d random_vector_sector(const double length, const double thetamin, const double thetamax);


/// Returns a random vector on a sphere surface of radius r.
Vector2d random_velocity(const double r);

/// Returns a vector of length r and orientation angle.
Vector2d vector_from_angle(const double angle, const double r);

/// Returns a vector of unit length and orientation angle.
inline Vector2d vector_from_angle(const double angle) { return vector_from_angle(angle,1.0); };

/// Returns the orientation angle of a vector.
inline double angle_from_vector(const topology::Vector2d& v) { return std::atan2(v.get_y(),v.get_x()); };

/// Returns a spin vector, i.e. unit vector, with given angle to x-axis
inline Vector2d spin(double theta) { return Vector2d(cos(theta), sin(theta)); };
/// Returns an orthogonal spin vector, i.e. a unit vector rotated 90^° counterclockwise from the vector of spin(theta)
inline Vector2d orthospin(double theta) { return Vector2d(-sin(theta), cos(theta)); };

Vector2d vector_on_squarelattice(int index, int Nx, int Ny, double spacing); ///< Index-dependent position vector for spin at index, square lattice
Vector2d vector_on_trigonallattice(int index, int Nx, int Ny, double spacing); ///< Index-dependent position vector for spin at index, trigonal lattice


/// Prints a list of vectors to in matlab-readable form
void print_matlab(const std::vector<Vector2d>& v, std::string name, std::ostream &outfile);


// =============================================================================================================================

/*! \class angle2d.
 *  \brief A (double-valued) angle in 2d space with the possibility of identifying it with its corresponding
 *  Vector2d unit vector. In the current state of the simulation, this is redundant.
 *
 *  Mostly incomplete and unnecessary. The idea was to have a specialized double-like class
 *  that can easily be converted to Vector2d and back for simplified calcultion. But the
 *  spin and orthospin functions defined for double -> Vector2d actually do the trick perfectly fine.
 *  Leaving this here for someone feeling a bit freaky.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: a somewhat uneventful weekend in late 2019
 *  \date Last Updated: 2023-08-01
 *
 */
class angle2d {
protected:
	/// Angle theta, to be interpreted as an angle with respect to the x-axis in a 2d xy-plane.
	double theta_;

public:
	/// Constructor without argument.
	angle2d() {} ;
	/// Constructor from a double.
	angle2d(const double& x);
	/// Copy constructor.
	angle2d(const angle2d& w);

	///< Conversion to double.
	inline operator double() const { return theta_ ;}
	/// Conversion to Vector2d - creates a unit vector with angle theta_ to x-axis.
	inline operator Vector2d() const { return Vector2d(cos(theta_), sin(theta_)) ;}

	/// Returns a spin vector, i.e. unit vector, with given angle to x-axis
	inline Vector2d spin() { return Vector2d(cos(theta_), sin(theta_)); };
	/// Returns an orthogonal spin vector, i.e. a unit vector rotated 90^° counterclockwise from the vector of spin(theta)
	inline Vector2d orthospin() { return Vector2d(-sin(theta_), cos(theta_)); };

	/// Resets the angle to be within (-pi, pi]
	void boundary();

	/// Addition by double
	angle2d& operator+=(const double& a);
	/// Subtraction of double
	angle2d& operator-=(const double& a);
	/// Multiplication by double
	angle2d& operator*=(const double a);
	/// Division by double
	angle2d& operator/=(const double a);
	/// Unary additive inversion
	angle2d operator-() const;
};

/// Addition operator
angle2d operator+(const angle2d& v, const angle2d& a);
/// Subtraction operator
angle2d operator-(const angle2d& v, const angle2d& w);
/// Multiplication operator (double, right)
angle2d operator*(const angle2d& v, const double a);
/// Multiplication operator (double, left)
angle2d operator*(const double a, const angle2d& v);
/// Division operator (double)
angle2d operator/(const angle2d& v, const double a);
}
#endif
