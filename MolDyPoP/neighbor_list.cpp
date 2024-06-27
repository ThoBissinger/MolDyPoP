/*! \file neighbor_list.cpp
 *
 *  \brief Cpp-file to class declaration of neighbor_list. Implements routines declared
 *  in neighbor_list.h
 *
 *  For more details, see neighbor_list.h
 *
 *
 *  \author Thomas Bissinger

 *  \date Created: mid 2019
 *  \date Last Updated: 2023-08-06
 */

#include "neighbor_list.h"
///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
neighbor_list::neighbor_list(const group& G) {
	nb_indices_.clear();
	nb_first_.clear();
	nb_dist_.clear();
	std::vector<int> nb;
	std::vector<double> distances;
	int index_first = 0;
	for (int index = 0; index < G.get_N(); index++){
		nb = G.get_neighbors(index, distances);
		nb_first_.push_back(index_first);
		index_first += nb.size();
		nb_indices_.insert(nb_indices_.end(),nb.begin(),nb.end());
		nb_dist_.insert(nb_dist_.end(),distances.begin(),distances.end());
	}
	nb_first_.push_back(nb_first_.back() + nb.size()); // Final index is set to N to avoid out-of-bounds errors
	N_ = G.get_N(); // Set at the very end because only then should isempty return false.
}

void neighbor_list::clear() {
	nb_indices_.clear();
	nb_first_.clear();
	nb_dist_.clear();
	N_ = 0;
}

///////////////////////////////////////////////////////////////////////////
// get_functions
///////////////////////////////////////////////////////////////////////////
std::vector<int> neighbor_list::get_neighbors(int i) const {
	// Error messages. Should be faster, but the other one is there for checking
	//return std::vector<int>(nb_indices_.begin() + nb_first_[i], nb_indices_.begin() + nb_first_[i+1]);
	std::vector<int> returnvec;
	for (int j = nb_first_[i]; j < nb_first_[i+1]; j++){
		returnvec.push_back(nb_indices_[j]);
	}
	return returnvec;
}
std::vector<double> neighbor_list::get_dist(int i) const {
	// Error messages. Should be faster, but the other one is there for checking
	// return std::vector<double>(nb_dist_.begin() + nb_first_[i], nb_dist_.begin() + nb_first_[i+1]);
	std::vector<double> returnvec;
	for (int j = nb_first_[i]; j < nb_first_[i+1]; j++){
		returnvec.push_back(nb_dist_[j]);
	}
	return returnvec;
}



