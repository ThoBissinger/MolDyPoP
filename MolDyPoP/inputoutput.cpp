/*! \file inputoutput.cpp
 *  \brief Implements the routines declared in inputoutput.h
 *
 *  Detailed description can be found in \link inputoutput.h inputoutput.h\endlink.
 *
 *  \author Thomas Bissinger
 *
 *  \date Created: mid-2017
 *  \date Last Update: 23-08-02
 *
 */
#include"inputoutput.h"
// TODO: Comma Separation optional (with additional string as input).
void print_stdvector(std::vector<double> vec, std::ostream& outfile){
	for (size_t i = 0; i < vec.size(); i++){
		outfile << vec[i] << ", ";
	}
	outfile << std::endl;
}

void print_stdvector(std::vector<double> vec, std::ostream& outfile, std::string startstring, std::string endstring){
	if (vec.size() > 0) {
		outfile << startstring;
		outfile << vec[0];
		for (size_t i = 1; i < vec.size(); i++){
				outfile << ", " << vec[i] ;
			}
		outfile << endstring;
	}
}

void print_stdvector(std::vector<double> vec, std::ostream& outfile, std::string startstring, std::string endstring,
		int index_start, int index_end, int stepskip){
	if (vec.size() > 0) {
		outfile << startstring;
		outfile << vec[index_start];
		for (int i = index_start + stepskip; (i < index_end) | (i < (int)vec.size()); i+= stepskip){
				outfile << ", " << vec[i] ;
			}
		outfile << endstring;
	}
}


void print_stdvector(std::vector<int> vec, std::ostream& outfile, std::string startstring, std::string endstring){
	if (vec.size() > 0) {
		outfile << startstring;
		outfile << vec[0];
		for (size_t i = 1; i < vec.size(); i++){
				outfile << ", " << vec[i] ;
			}
		outfile << endstring;
	}
}

void print_stdvector(std::vector<std::complex<double> > vec, std::ostream& outfile, std::string startstring, std::string endstring,
		int index_start, int index_end, int stepskip){
	if (vec.size() > 0) {
		outfile << startstring;
		outfile << vec[index_start].real() << " + i * " << vec[index_start].imag();
		for (int i = index_start + stepskip; (i < index_end) | (i < (int)vec.size()); i+= stepskip){
			outfile << ", " << vec[i].real() << " + i * " << vec[i].imag() ;
		}
		outfile << endstring;
	}
}
