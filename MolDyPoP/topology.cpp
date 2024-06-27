/*! \file topology.cpp
 *
 *  \brief Cpp-File to declaration of namespace topology. Implements routines.
 *
 *  Provides features for handling spatial vectors both for position
 *  and spin vectors. Not many detailed comments, most comments can
 *  be found in the corresponding file topology.h.
 */
#include "topology.h"
#include "computations.h"
#include "math.h"
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
// =====================================================================
// =====================================================================
// CLASS VECTOR2D
// =====================================================================
// =====================================================================

// =====================================================================
// Constructors, copy constructor
// =====================================================================
// Constructor with no input, one-element vector initialized.
topology::Vector2d::Vector2d(){}

// Constructor with one double for both elements. Basically only used
// to get a zero-vector.
topology::Vector2d::Vector2d(const double& x) {
  v_[0]=x;
  v_[1]=x;
}

// Constructor with one int for both elements. Basically only used
// to get a zero-vector.
topology::Vector2d::Vector2d(const int& x) {
  v_[0]=x;
  v_[1]=x;
}

// Constructor with two doubles, x and y.
topology::Vector2d::Vector2d(const double x, const double y) {
  v_[0]=x;
  v_[1]=y;
}


// Copy constructor
topology::Vector2d::Vector2d(const topology::Vector2d& w) {
  v_[0]=w.get_x();
  v_[1]=w.get_y();
}
// =====================================================================

// Print function
void topology::Vector2d::print(std::ostream &outfile) const {
  outfile.setf(std::ios_base::fixed, std::ios_base::floatfield);
  outfile.precision(7);
  outfile << std::setw(13) << v_[0] << " " << std::setw(13) << v_[1] <<  std::setw(13) << "0.0000" << std::endl;
}

void topology::Vector2d::print() const {
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  std::cout.precision(7);
  std::cout << std::setw(13) << v_[0] << " " << std::setw(13) << v_[1] << std::endl;
}

void topology::Vector2d::print_xyz() const {
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  std::cout.precision(7);
  std::cout << std::setw(13) << v_[0] << " " << std::setw(13) << v_[1] <<  std::setw(13) << "0.0000" << std::endl;
}

// =====================================================================
// Standard vector assessment

// Returns pointer to element.
double& topology::Vector2d::operator[](int index) {
	return v_[index];
}
// =====================================================================

// Returns a pointer to const element.
const double& topology::Vector2d::operator[](int index) const {
	return v_[index];
}

// =====================================================================
// Value assignment
void topology::Vector2d::set_x(double x) {
	v_[0] = x;
}

void topology::Vector2d::set_y(double y) {
	v_[1] = y;
}

// =====================================================================
// MEMBER-OPERATORS
// =====================================================================
// Add another vector
topology::Vector2d& topology::Vector2d::operator+=(const topology::Vector2d& w){
	v_[0]+=w.get_x();
	v_[1]+=w.get_y();
	return *this;
}

// Subtract another vector
topology::Vector2d& topology::Vector2d::operator-=(const topology::Vector2d& w){
	v_[0]-=w.get_x();
	v_[1]-=w.get_y();
	return *this;
}

// Multiply by a scalar (double)
topology::Vector2d& topology::Vector2d::operator*=(const double a){
	v_[0]*=a;
	v_[1]*=a;
	return *this;
}

// Multiply by a scalar (int)
topology::Vector2d& topology::Vector2d::operator*=(const int a){
	v_[0]*=a;
	v_[1]*=a;
	return *this;
}


// Divide by a scalar (double)
topology::Vector2d& topology::Vector2d::operator/=(const double a){
	v_[0]/=a;
	v_[1]/=a;
	return *this;
}

// Divide by a scalar (int)
topology::Vector2d& topology::Vector2d::operator/=(const int a){
	v_[0]/=a;
	v_[1]/=a;
	return *this;
}

// Unary Minus
topology::Vector2d topology::Vector2d::operator-() const{
	return topology::Vector2d(-v_[0],-v_[1]);
}
// =====================================================================

// =====================================================================
// NON-MEMBER OPERATORS
// =====================================================================
topology::Vector2d topology::operator+(const topology::Vector2d& v, const topology::Vector2d& w)
{
  return topology::Vector2d(v) += w;
}

topology::Vector2d topology::operator-(const topology::Vector2d& v, const topology::Vector2d& w)
{
  return topology::Vector2d(v) -= w;
}

topology::Vector2d topology::operator*(const topology::Vector2d& v, const double a)
{
  return topology::Vector2d(v) *= a;
}

topology::Vector2d topology::operator*(const double a, const topology::Vector2d& v)
{
  return v * a;
}

topology::Vector2d topology::operator*(const topology::Vector2d& v, const int a)
{
  return topology::Vector2d(v) *= a;
}

topology::Vector2d topology::operator*(const int a, const topology::Vector2d& v)
{
  return v * a;
}

topology::Vector2d topology::operator/(const topology::Vector2d& v, const double a)
{
  return topology::Vector2d(v) /= a;
}

topology::Vector2d topology::operator/(const topology::Vector2d& v, const int a)
{
  return topology::Vector2d(v) /= a;
}
// =====================================================================


// =====================================================================
// PERIODIC DISTANCES IN A SQUARE BOX OF LENGTH L
// =====================================================================
double topology::periodic_distance_squared(const Vector2d& v, const Vector2d& w, const double& L){
	double abs_x = std::abs(v.get_x() - w.get_x());
	double abs_y = std::abs(v.get_y() - w.get_y());
	if (abs_x > .5 * L){
		abs_x -= L;
	}
	if (abs_y > .5 * L){
		abs_y -= L;
	}
	return abs_x * abs_x + abs_y * abs_y;
}

topology::Vector2d topology::periodic_distance_vector(const Vector2d& v, const Vector2d& w, const double& L){
	topology::Vector2d distvec = w - v;
	if (distvec.get_x() > .5 * L){
		distvec -= topology::Vector2d(L, 0 );
	} else if (distvec.get_x() < -.5 * L){
		distvec += topology::Vector2d(L, 0 );
	}
	if (distvec.get_y() > .5 * L){
		distvec -= topology::Vector2d(0, L);
	} else if (distvec.get_y() < -.5 * L){
		distvec += topology::Vector2d(0, L);
	}
	return distvec;
};


// PERIODIC DISTANCES IN A RECTANGULAR BOX
double topology::periodic_distance_squared(const Vector2d& v, const Vector2d& w, const Vector2d& L){
	double abs_x = std::abs(v.get_x() - w.get_x());
	double abs_y = std::abs(v.get_y() - w.get_y());
	if (abs_x > .5 * L.get_x()){
		abs_x -= L.get_x();
	}
	if (abs_y > .5 * L.get_y()){
		abs_y -= L.get_y();
	}
	return abs_x * abs_x + abs_y * abs_y;
}

topology::Vector2d topology::periodic_distance_vector(const Vector2d& v, const Vector2d& w, const Vector2d& L){
	topology::Vector2d distvec = w - v;
	if (distvec.get_x() > .5 * L.get_x()){
		distvec -= topology::Vector2d(L.get_x(), 0 );
	} else if (distvec.get_x() < -.5 * L.get_x()){
		distvec += topology::Vector2d(L.get_x(), 0 );
	}
	if (distvec.get_y() > .5 * L.get_y()){
		distvec -= topology::Vector2d(0, L.get_y());
	} else if (distvec.get_y() < -.5 * L.get_y()){
		distvec += topology::Vector2d(0, L.get_y());
	}
	return distvec;
};
// =====================================================================







// =====================================================================
// ROTATIONS
// =====================================================================
// Rotate. Member function, rotates by given angle
void topology::Vector2d::rotate(const double theta){
  double oldv0 = v_[0];
  v_[0]=cos(theta)*v_[0]-sin(theta)*v_[1];
  v_[1]=sin(theta)*oldv0+cos(theta)*v_[1];
}
// =====================================================================

// =====================================================================
// BOUNDARY CONDITIONS
// =====================================================================
// periodic_box: Enforces periodic boundary conditions (in a box)
  void topology::Vector2d::periodic_box(const Vector2d& minima, const Vector2d& maxima){
    if (v_[0] < minima.get_x()){
      v_[0] = maxima.get_x() - minima.get_x() +  v_[0];
    } else if (v_[0] >= maxima.get_x()){
      v_[0] = - maxima.get_x() + minima.get_x() +  v_[0];
    }
    if (v_[1] < minima.get_y()){
      v_[1] = maxima.get_y() - minima.get_y() +  v_[1];
    } else if (v_[1] >= maxima.get_y()){
      v_[1] = - maxima.get_y() + minima.get_y() +  v_[1];
    }
  }

  topology::Vector2d topology::Vector2d::get_boundary_handler_periodic_box
     (const topology::Vector2d& minima, const topology::Vector2d& maxima,
      const double cutoff) const
  {
    topology::Vector2d handler=0.0;
    for (int i=0; i < 2; i++){
      if (v_[i] - minima[i] < cutoff){
        handler[i]=-1;
      } else if (maxima[i] - v_[i] < cutoff){
        handler[i]=1;
      }
    }
    return handler;
  }
// =====================================================================

// =====================================================================
// BOUNDARY CONDITIONS, ZERO BORDER
// =====================================================================
// periodic_box: Enforces periodic boundary conditions (in a box)
  void topology::Vector2d::periodic_box(const Vector2d& maxima){
    if (v_[0] < 0){
      v_[0] = maxima.get_x() + v_[0];
    } else if (v_[0] >= maxima.get_x()){
      v_[0] = - maxima.get_x() +  v_[0];
    }
    if (v_[1] < 0){
      v_[1] = maxima.get_y() +  v_[1];
    } else if (v_[1] >= maxima.get_y()){
      v_[1] = - maxima.get_y() +  v_[1];
    }
  }

  topology::Vector2d topology::Vector2d::get_boundary_handler_periodic_box
     (const topology::Vector2d& maxima, const double cutoff) const {
    topology::Vector2d handler=0.0;
    for (int i=0; i < 2; i++){
      if (v_[i] < cutoff){
        handler[i]=-1;
      } else if (maxima[i] - v_[i] < cutoff){
        handler[i]=1;
      }
    }
    return handler;
  }
// =====================================================================


// =====================================================================
// NORM FUNCTIONS
// =====================================================================
// normalized. Changes to normalized value.
void topology::Vector2d::normalized() {
  *this/=sqrt(norm2());
}

topology::Vector2d topology::vector_from_angle(const double angle, const double r){
  return topology::Vector2d(r*cos(angle),r*sin(angle));
}
topology::Vector2d topology::vector_on_squarelattice(int index, int Nx, int Ny, double spacing) {
	int x = index % Nx;
	int y = (index - x) / Ny;
	return spacing * topology::Vector2d(x,y);
}
topology::Vector2d topology::vector_on_trigonallattice(int index, int Nx, int Ny, double spacing) {
	int x = index % Nx;
	int y = (index - x) / Ny;
	if ( y % 2 == 0){
	return spacing * topology::Vector2d(x, .5 * sqrt(3) * y);
	} else {
	return spacing * topology::Vector2d(x + .5 ,.5 * sqrt(3) * y);
	}
}

topology::Vector2d topology::random_vector(const topology::Vector2d& minima, const topology::Vector2d& maxima){
  return topology::Vector2d(random_double(minima.get_x(),maxima.get_x()),random_double(minima.get_y(),maxima.get_y()));
}

topology::Vector2d topology::random_vector(const topology::Vector2d& maxima){
  return topology::Vector2d(random_double(maxima.get_x()),random_double(maxima.get_y()));
}
topology::Vector2d topology::random_gaussian_vector(const double& sigma_squared){
	return topology::Vector2d(random_boltzmann_double(sigma_squared),random_boltzmann_double(sigma_squared));
}

topology::Vector2d topology::random_vector_first_quadrant(const double length){
  double theta=random_double(0,.5*M_PI);
  return topology::Vector2d(length*cos(theta),length*sin(theta));
}

topology::Vector2d topology::random_vector_sector(const double length, const double thetamin, const double thetamax){
  double theta=random_double(thetamin,thetamax);
  return topology::Vector2d(length*cos(theta),length*sin(theta));
}


topology::Vector2d topology::random_velocity(const double vel){
  double theta=random_double(0,2*M_PI);
  return topology::Vector2d(vel*cos(theta),vel*sin(theta));
}


topology::Vector2d topology::nearest_gridvec(topology::Vector2d v, double gridsep){
	topology::Vector2d w = 0.0;
	for(int i = 0; i < 2; i++){
		w[i] = std::round(v[i] / gridsep) * gridsep;
	}
	return w;
}

topology::Vector2d topology::nearest_gridvec(Vector2d v, Vector2d gridseps){
	topology::Vector2d w = 0.0;
	for(int i = 0; i < 2; i++){
		w[i] = std::round(v[i] / gridseps[i]) * gridseps[i];
	}
	return w;
}

std::vector<topology::Vector2d> topology::qvalues_within_radius(double qmin, double qmax, topology::Vector2d gridseps, int qsamps_per_bin){
	std::vector<topology::Vector2d> w;
	topology::Vector2d candidatevec;
	bool fillcheck;
	int count = 0; int trials = 0;
	while( (count < qsamps_per_bin) & (trials < 10 * qsamps_per_bin) ){ // Breaks if more than qsamps_per_bin trials have been made (unlikely to find another value).
		candidatevec = nearest_gridvec(random_vector_sector(random_double(qmin, qmax),-.5*M_PI,.5*M_PI), gridseps);
		if ( (candidatevec.norm2() < qmax * qmax) & (candidatevec.norm2() > qmin * qmin)){
			fillcheck = 0;
			for(size_t i = 0; ( i < w.size() ) & (!fillcheck); i++){
				if ( (w[i].get_x() == candidatevec.get_x()) & (w[i].get_y() == candidatevec.get_y())){
					fillcheck = 1;
				}
			}
			if ( !fillcheck){
				w.push_back(candidatevec);
				count++;
			}
		}
		trials++;
	}
	return w;
}


void topology::print_matlab(const std::vector<Vector2d>& v, std::string name, std::ostream &outfile) {
	for (int i = 0; i < 2; i++){
		outfile << name << "(" << i+1 << ",:) = [" << (v[0])[i];
		for (size_t k = 1; k < v.size(); k++){
			outfile << ", " << (v[k])[i] ;
		}
		outfile << "];" << std::endl;
	}
}

// =====================================================================
// =====================================================================
// CLASS ANGLE2D
// =====================================================================
// =====================================================================

// =====================================================================
// Constructors, copy constructor
// =====================================================================
// Constructor from double.
topology::angle2d::angle2d(const double& theta) {
  theta_=theta;
}

// Copy constructor
topology::angle2d::angle2d(const topology::angle2d& ang) {
  theta_=ang.theta_;
}

// =====================================================================
// Boundary Conditions, i.e. angle in (-pi,pi]
// =====================================================================
  void topology::angle2d::boundary() {
    while (theta_ > M_PI){
      theta_ -= 2 * M_PI;
    }
    while (theta_ <= -M_PI){
      theta_ += 2 * M_PI;
    }
  }

  // =====================================================================
  // MEMBER-OPERATORS
  // =====================================================================
  // Add another vector
  topology::angle2d& topology::angle2d::operator+=(const double& a){  // return type topology::angle2d&, & damit dass return das objekt nicht kopiert sondern ändert
    theta_ += a;
    return *this;
  }

  // Subtract another vector
  topology::angle2d& topology::angle2d::operator-=(const double& a){
	theta_ -= a;
    return *this;
  }

  // Multiply by a scalar (double)
  topology::angle2d& topology::angle2d::operator*=(const double a){
	theta_ *= a;
    return *this;
  }

  // Divide by a scalar (double)
  topology::angle2d& topology::angle2d::operator/=(const double a){
	theta_ /= a;
    return *this;
  }

  // Unary Minus
  topology::angle2d topology::angle2d::operator-() const{
    return angle2d(-theta_);
  }
  // =====================================================================

  // =====================================================================
  // NON-MEMBER OPERATORS
  // =====================================================================
  topology::angle2d topology::operator+(const topology::angle2d& ang1, const angle2d& ang2)
  {
    return topology::angle2d(ang1) += ang2;
  }

  topology::angle2d topology::operator-(const topology::angle2d& ang1, const topology::angle2d& ang2)
  {
    return topology::angle2d(ang1) -= ang2;
  }

  topology::angle2d topology::operator*(const topology::angle2d& ang, const double a)
  {
    return topology::angle2d(ang) *= a;
  }

  topology::angle2d topology::operator*(const double a, const topology::angle2d& ang)
  {
    return ang * a;
  }

  topology::angle2d topology::operator/(const topology::angle2d& v, const double a)
  {
    return topology::angle2d(v) /= a;
  }
