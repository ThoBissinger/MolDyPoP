/*! \file computations.cpp
 *  \brief cpp-File to computations.h, implementation of the functions.
 *
 *  For more details, see \link computations.h computations.h\endlink.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: mid-2017
 *  \date Last Update: 23-08-02
 *
 */

#include"computations.h"
#include<iostream>
#include<stdlib.h>

// =======================================================================================================
// general functions
// =======================================================================================================
double mod(double val, const double min, const double max){
	while (val > max){
		val -= max - min;
	}
	while (val <= min){
		val += max - min;
	}
	return val;
}
// =======================================================================================================
// vector functions
// =======================================================================================================

// =============================
// Non-selective means
// =============================
template<typename T, typename A>
T vector_mean(std::vector<T,A> v){
	T meanval = 0;
	for (size_t i = 0; i < v.size(); i++){
		meanval += v[i];
	}
	return meanval / (double)v.size();
}
// Instantiations. Otherwise, the code would have to be written in the h file.
template double vector_mean(std::vector<double, std::allocator<double> > v);
template int vector_mean(std::vector<int, std::allocator<int> > v);
template std::complex<double> vector_mean(std::vector<std::complex<double>, std::allocator<std::complex<double> > > v);
template topology::Vector2d vector_mean(std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > v);


template<typename T, typename A>
T row_mean(std::vector<T,A> v, int M, int N, int i){
	T meanval = 0;
	for (int k = 0; k < N; k++){
		meanval += v[matr_index(M, N, i, k)];
	}
	return meanval * (1.0/ (double)N); // conversion necessary to avoid problems with tempalte for std::complex<double>
}
// Instantiation. As above
template double row_mean(std::vector<double, std::allocator<double> > v, int M, int N, int i);
template int row_mean(std::vector<int, std::allocator<int> > v, int M, int N, int i);
template std::complex<double> row_mean(std::vector<std::complex<double>, std::allocator<std::complex<double> > > v, int M, int N, int i);
template topology::Vector2d row_mean(std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > v, int M, int N, int i);

template<typename T, typename A>
T column_mean(std::vector<T,A> v, int M, int N, int j){
	T meanval = 0;
	for (int k = 0; k < M; k++){
		meanval += v[matr_index(M, N, k, j)];
	}
	return meanval * (1.0/ (double)M); // conversion necessary to avoid problems with tempalte for std::complex<double>
}
// Instantiation. As above
template double column_mean(std::vector<double, std::allocator<double> > v, int M, int N, int j);
template int column_mean(std::vector<int, std::allocator<int> > v, int M, int N, int j);
template std::complex<double> column_mean(std::vector<std::complex<double>, std::allocator<std::complex<double> > > v, int M, int N, int j);
template topology::Vector2d column_mean(std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > v, int M, int N, int j);





// =============================
// selective means
// =============================
template<typename T, typename A>
T selective_vector_mean(std::vector<T, A> yvals, std::vector<double> xvals, double xsep_min){
	T meanval = 0;
	int meancount = 0;
	double xlast = xvals[0] - xsep_min;
	for (size_t i = 0; i < xvals.size(); i++){
		if (xvals[i] - xlast >= xsep_min){
			meanval += yvals[i];
			xlast = xvals[i];
			meancount ++;
		}
	}
	return meanval / (double)meancount;
}
// Instantiations. Otherwise, the code would have to be written in the h file.
template double selective_vector_mean(std::vector<double, std::allocator<double> > yvals, std::vector<double> xvals, double xsep_min);
template int selective_vector_mean(std::vector<int, std::allocator<int> > yvals, std::vector<double> xvals, double xsep_min);
template std::complex<double> selective_vector_mean(std::vector<std::complex<double>, std::allocator<std::complex<double> > > yvals, std::vector<double> xvals, double xsep_min);
template topology::Vector2d selective_vector_mean(std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > yvals, std::vector<double> xvals, double xsep_min);

template<typename T, typename A>
T selective_row_mean(std::vector<T,A> yvals, std::vector<double> xvals, double xsep_min, int M, int N, int i){
	T meanval = 0;
	int meancount = 0;
	double xlast = xvals[0] - xsep_min;
	for (size_t k = 0; k < xvals.size(); k++){
		if (xvals[k] - xlast >= xsep_min){
			meanval += yvals[matr_index(M, N, i, k)];
			xlast = xvals[k];
			meancount ++;
		}
	}
	return meanval / (double)meancount;
}
// Instantiation. As above
template double selective_row_mean(std::vector<double, std::allocator<double> > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int i);
template int selective_row_mean(std::vector<int, std::allocator<int> > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int i);
template std::complex<double> selective_row_mean(std::vector<std::complex<double>, std::allocator<std::complex<double> > > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int i);
template topology::Vector2d selective_row_mean(std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int i);

template<typename T, typename A>
T selective_column_mean(std::vector<T,A> yvals, std::vector<double> xvals, double xsep_min, int M, int N, int j){
	T meanval = 0;
	int meancount = 0;
	double xlast = xvals[0] - xsep_min;
	for (size_t k = 0; k < xvals.size(); k++){
		if (xvals[k] - xlast >= xsep_min){
			meanval += yvals[matr_index(M, N, k, j)];
			xlast = xvals[k];
			meancount ++;
		}
	}
	return meanval / (double)meancount;
}
// Instantiation. As above
template double selective_column_mean(std::vector<double, std::allocator<double> > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int j);
template int selective_column_mean(std::vector<int, std::allocator<int> > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int j);
template std::complex<double> selective_column_mean(std::vector<std::complex<double>, std::allocator<std::complex<double> > > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int j);
template topology::Vector2d selective_column_mean(std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > yvals, std::vector<double> xvals, double xsep_min, int M, int N, int j);



double vector_variance(std::vector<double> v, double& mean){
	mean = 0;
	double variance = 0;
	for (size_t i = 0; i < v.size(); i++){
		mean += v[i];
		variance += v[i]*v[i];
	}
	mean /= v.size();
	return variance / v.size() - mean * mean;
}

template<typename T, typename A>
std::vector<T,A> vector_scale(const std::vector<T,A>& v, double a){
	std::vector<T,A> w = v;
	for (size_t i = 0; i < v.size(); i++){
		w[i] *= a;
	}
	return w;
}
// Instantiation. As above
template std::vector<double, std::allocator<double> > vector_scale(const std::vector<double, std::allocator<double> >& v, double a);
template std::vector<int, std::allocator<int> > vector_scale(const std::vector<int, std::allocator<int> >& v, double a);
template std::vector<std::complex<double>, std::allocator<std::complex<double> > > vector_scale(const std::vector<std::complex<double>, std::allocator<std::complex<double> > >& v, double a);
template std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > vector_scale(const std::vector<topology::Vector2d, std::allocator<topology::Vector2d > >& v, double a);

template<typename T, typename A>
std::vector<T,A> vector_sum(const std::vector<T,A>& v, const std::vector<T,A>& w){
	std::vector<T,A> u = v;
	for (size_t i = 0; i < v.size(); i++){
		u[i] += w[i];
	}
	return u;
}
// Instantiation. As above
template std::vector<double, std::allocator<double> > vector_sum(const std::vector<double, std::allocator<double> >& v, const std::vector<double, std::allocator<double> >& w);
template std::vector<int, std::allocator<int> > vector_sum(const std::vector<int, std::allocator<int> >& v, const std::vector<int, std::allocator<int> >& w);
template std::vector<std::complex<double>, std::allocator<std::complex<double> > > vector_sum(const std::vector<std::complex<double>, std::allocator<std::complex<double> > >& v, const std::vector<std::complex<double>, std::allocator<std::complex<double> > >& w);
template std::vector<topology::Vector2d, std::allocator<topology::Vector2d > > vector_sum(const std::vector<topology::Vector2d, std::allocator<topology::Vector2d > >& v, const std::vector<topology::Vector2d, std::allocator<topology::Vector2d > >& w);


template<typename T, typename A>
T vector_inpr(const std::vector<T,A>& v, const std::vector<T,A>& w){
	T mult_sum = 0;
	for (size_t i = 0; i < v.size(); i++){
		mult_sum += v[i] * w[i];
	}
	return mult_sum;
}
template double vector_inpr(const std::vector<double, std::allocator<double> >& v, const std::vector<double, std::allocator<double> >& w);
template int vector_inpr(const std::vector<int, std::allocator<int> >& v, const std::vector<int, std::allocator<int> >& w);

std::complex<double> vector_inpr(const std::vector<std::complex<double> >& v, const std::vector<std::complex<double> >& w){
	std::complex<double> mult_sum = 0;
	for (size_t i = 0; i < v.size(); i++){
		mult_sum += std::conj(v[i]) * w[i];
	}
	return mult_sum;
}

template<typename T, typename A>
std::vector<T,A> vector_pw_mult(const std::vector<T,A>& v, const std::vector<T,A>& w){
	std::vector<T,A> mult_vec = w;
	for (size_t i = 0; i < v.size(); i++){
		mult_vec[i] *= v[i];
	}
	return mult_vec;
}
template std::vector<double, std::allocator<double> > vector_pw_mult(const std::vector<double, std::allocator<double> >& v, const std::vector<double, std::allocator<double> >& w);
template std::vector<int, std::allocator<int> > vector_pw_mult(const std::vector<int, std::allocator<int> >& v, const std::vector<int, std::allocator<int> >& w);

std::vector<std::complex<double> > vector_pw_mult(const std::vector<std::complex<double> >& v, const std::vector<std::complex<double> >& w){
	std::vector<std::complex<double> > mult_vec = w;
	for (size_t i = 0; i < v.size(); i++){
		mult_vec[i] *= std::conj(v[i]);
	}
	return mult_vec;
}

std::vector<std::complex<double> > vector_pw_mult(const std::vector<std::complex<double> >& v_1, const std::vector<std::complex<double> >& v_2,
		const std::vector<std::complex<double> >& v_3, const std::vector<std::complex<double> >& v_4){
	std::vector<std::complex<double> > mult_vec = v_3;
	for (size_t i = 0; i < v_3.size(); i++){
		mult_vec[i] *= std::conj(v_1[i]) * std::conj(v_2[i]) * v_4[i];
	}
	return mult_vec;
}

std::vector<double> vector_pw_norm(const std::vector<std::complex<double> >& v){
	std::vector<double> mult_vec;
	for (size_t i = 0; i < v.size(); i++){
		mult_vec.push_back(std::norm(v[i]));
	}
	return mult_vec;
}

// =======================================================================================================
// functions for binning
// =======================================================================================================
std::vector<double> log_bin(double binmin, double binmax, double diffmin, int Nbins){
	std::vector<double> binvals;
	double nu = (binmax - binmin) / diffmin;
	double exp_factor = 2;
	while ( abs(std::pow(exp_factor, Nbins) - nu * exp_factor + nu - 1) > 1e-3 ){
		exp_factor -= (std::pow(exp_factor, Nbins) - nu * exp_factor + nu - 1)/(Nbins * std::pow(exp_factor,Nbins-1) - nu);
	}
	binvals.push_back(binmin + diffmin);
	binvals.push_back(binmin + (1 + exp_factor) * diffmin);
	if ( ( exp_factor > 1 )  &  ( binvals.back() < binmax ) ){
		while (binvals.back() * exp_factor < binmax ){
			binvals.push_back(binvals.back() + (binvals.back() - binvals[binvals.size()-2]) * exp_factor);
//			binvals.push_back((binvals.back() ) * exp_factor);
		}
	} else {
		std::cout << "PROBLEM in computations, log_bin: exp_factor = " << exp_factor << ", binmax = " << binmax
				<< ", binmin + diffmin = " << binmin + diffmin
				<< ".\n"
				<< "Check input values!" << std::endl;
	}
	if (binvals.back() < binmax){
		// The way we handle bins demands that we have to add the max value at the end.
		// The value itself is a dummy. The bin ist just filled with anything larger than
		// the pre-to-last entry.
		binvals.push_back(binmax);
	} else if ( binvals.back() > binmax){
		binvals.back() = binmax;
	}
	return binvals;
}


int find_bin(const std::vector<double>& bin, const double& value){
	int index = 0;
	// Any entry that is smaller than the first entry of bin enters there,
	// then the resulting index is the index i for which bin[i - 1] < value < bin[i].
	// Note that the bin has to have one extra element for this to work,
	// see also comment at the end of lin_log_bin().
	while ( (index < (int)(bin.size()) - 1) & ( value > bin[index])){
		index++;
	}
	return index;
}



// =======================================================================================================
// random functions
// =======================================================================================================
double random_angle(double min, double max) {
	if(max>M_PI){
		max=M_PI;
	}
	if(min<-M_PI){
		min=-M_PI;
	}
  return min + (double)rand() / RAND_MAX * (max - min);
}

double random_angle(double max) {
	if(max>M_PI){
		max=M_PI;
		  return (double)rand() / RAND_MAX * 2*(max)-M_PI;
	}
	else{
		if(max<-M_PI){
			max=-M_PI;
			return (double)rand() / RAND_MAX * 2*(-max)-M_PI;
		}
	}

	if(max>0){
			  return (double)rand() / RAND_MAX * 2*(max)-max;
		}
		else{
		if(max<0){
			  return (double)rand() / RAND_MAX * 2*(-max)+max;
		}
		}
	return 0;
}

