/*! \file partition.cpp
 *
 *  \brief Cpp-file to class declaration of partition. Implements routines
 *  defined in partition.h.
 *
 *  For details, see partition.h. Most code here is more or less eslf-explanatory.
 *  Some use of standard library functions is made for vector handling etc.
 *
 *  \author Thomas Bissinger

 *  \date Created: late 2017
 *  \date Last Updated: 2023-08-06
 */

#include "partition.h"

#include <iomanip>
#include <fstream>
#include <chrono>


partition::partition(const int& N, const double& cutoff, const topology::Vector2d& L):
N_(N), L_(L) {
	cellnum_ = 1;
	for (int i=0; i < L.size(); i++){
		M_[i] = floor(L_[i] / cutoff);
		cellwidth_[i] = L_[i] / M_[i];
		cellnum_ *= M_[i];
	}
   cellelem_ = std::vector<int>(cellnum_, 0);
   firsts_ = std::vector<int>(cellnum_, 0);
   cellvec_.reserve(N_);
}

partition::partition(const int& N, const double& cutoff,
		const topology::Vector2d& L,
		const std::vector<topology::Vector2d>& positions)
: partition(N, cutoff, L) {
  clear();
  fill(positions);
}



partition::~partition() {
  clear();
}


void partition::clear() {
	cellvec_.clear();
	cellelem_.clear();
	firsts_.clear();
	for (int i=0; i < cellnum_ ; i++){
		cellelem_.push_back(0);
		firsts_.push_back(0);
	}
}




void partition::fill(const std::vector<topology::Vector2d>& positions) {
	clear();
	std::vector<int> cell_index_vector; // stores the cell indices for each particle
	std::vector<int> filling_vector(cellnum_,0); // stores how far the cells are filled (used later)
	int index;
	for (int i = 0; i < N_; i++){
		index = find_cell(positions[i]);
		cell_index_vector.push_back(index);
		cellelem_[index]++;
	}
	firsts_[0] = 0;
	for (int i = 1; i < cellnum_; i++){
		firsts_[i] = firsts_[i-1] + cellelem_[i-1];
	}
	cellvec_ = std::vector<int>(N_);
	for (int i = 0; i < N_; i++){
		index = cell_index_vector[i];
		cellvec_[firsts_[index] + filling_vector[index]] = i;
		filling_vector[index]++;
	}
}


int partition::find_cell(const topology::Vector2d& r) const{
   double factor = 1;
   int k=1;
   int index=0;
   double curr;
   for(int j=0; j < r.size(); j++){
      curr = r[j];
      while (k * cellwidth_[j] < curr){
         k++;
      }
      index+=factor*(k-1);
      factor*=M_[j];
      k=1;
   }
   if (index <= N_ - 1 && index >= 0){
	   return index;
   } else {
	   std::cout << "++ ERROR in partition find_cell ++" << std::endl;
	   std::cout << "++       cell index:    " << index << std::endl;
	   std::cout << "++       particle position:    " << std::endl;
	   r.print();
	   exit(-2);
   }
}

// Computes neighboring cell to cell at index, lr in left-right (x) and ud in up-down (y) direction, i,j are in {-1,0,1}. Uses periodic boundary conditions.
int partition::neighbor_cell(int index, int lr, int ud) const {
   int boundary = 0;
   if (lr == -1 && index % M_[0] == 0){
      boundary = M_[0];
   } else if (lr == 1 && index % M_[0] == M_[0] - 1){
      boundary = -M_[0];
   }
   if (index + M_[0] * ud < 0 ){
      boundary += cellnum_;
   } else if (index + M_[0] * ud >= cellnum_){
      boundary -= cellnum_;
   }

  return index + lr + M_[0]*ud + boundary;
}

// Computes neighboring cell to cell at index, lr in left-right (x) and ud in up-down (y) direction, i,j are in {-1,0,1}. Uses periodic boundary conditions. Includes a shift.
int partition::neighbor_cell(int index, int lr, int ud, topology::Vector2d& shift) const {
   shift=0;
   int boundary = 0;
   if (lr == -1 && index % M_[0] == 0){
      boundary = M_[0];
      shift[0] = -L_[0];
   } else if (lr == 1 && index % M_[0] == M_[0] - 1){
      boundary = -M_[0];
      shift[0] = L_[0];
   }
   if (index + M_[0] * ud < 0 ){
      boundary += cellnum_;
      shift[1] = -L_[1];
   } else if (index + M_[0] * ud >= cellnum_){
      boundary -= cellnum_;
      shift[1] = L_[1];
   }

  return index + lr + M_[0]*ud + boundary;
}

std::vector<int> partition::nb_cells_all(int index, std::vector<topology::Vector2d>& shifts) const {
	shifts.clear();
	topology::Vector2d shift_cur;
	std::vector<int> neivec;
	for (int lr = -1; lr < 2; lr++){
		for (int ud = -1; ud < 2; ud++){
			neivec.push_back(neighbor_cell(index, lr, ud, shift_cur));
			shifts.push_back(shift_cur);
		}
	}
	return neivec;
}

std::vector<int> partition::nb_cells_all(const topology::Vector2d& r, double cutoffsquared,
		std::vector<topology::Vector2d>& shifts) const {
	shifts.clear();
	topology::Vector2d shift_cur;
	int index = find_cell(r);
	std::vector<int> neivec;
	int nei_candidate;
	// the cell itself
	if (cellelem_[index] > 1) {
		neivec.push_back(index);
		shifts.push_back(0.0);
	}
	for (int lr = -1; lr < 2; lr++){
		for (int ud = -1; ud < 2; ud++){
			if ( (lr != 0) || (ud != 0) ){ // ignores the cell itself, which was already covered
				nei_candidate = neighbor_cell(index, lr, ud, shift_cur);
				if ( cellelem_[nei_candidate] > 0){ // does not include empty cells
					if ( ((lr == 0) || (ud == 0)) // cells at edges always included
							|| (norm2(r - corner(index,lr,ud)) < cutoffsquared ) ){
						// if the particle under consideration is too far from the closest point in the cell
						// (the corner) it is not included.
						neivec.push_back(neighbor_cell(index, lr, ud, shift_cur));
						shifts.push_back(shift_cur);
					}
				}
			}

		}
	}
	return neivec;
}

std::vector<int> partition::nb_cells_ur(int index, std::vector<topology::Vector2d>& shifts) const {
	shifts.clear();
	topology::Vector2d shift_cur = 0;
	std::vector<int> neivec;

	neivec.push_back(index); // cell
	shifts.push_back(0.0);
	neivec.push_back(neighbor_cell(index, 1, 0, shift_cur)); // right cell
	shifts.push_back(shift_cur);
	neivec.push_back(neighbor_cell(index, -1, 1, shift_cur)); // top left cell

	shifts.push_back(shift_cur);
	neivec.push_back(neighbor_cell(index, 0, 1, shift_cur)); // top center cell
	shifts.push_back(shift_cur);
	neivec.push_back(neighbor_cell(index, 1, 1, shift_cur)); // top right cell
	shifts.push_back(shift_cur);

	return neivec;
}


std::vector<int> partition::nb_cells_ur(const topology::Vector2d& r, double cutoffsquared,
		std::vector<topology::Vector2d>& shifts) const {
	shifts.clear();
	topology::Vector2d shift_cur = 0;
	int index = find_cell(r);
	std::vector<int> neivec;
	int nei_candidate;

	// The cell itself (only contains a relevant neighbor if it has more than one entry)
	if (cellelem_[index] > 1 ){
		neivec.push_back(index);
		shifts.push_back(shift_cur);
	}
	// The neighboring cell to the right
	nei_candidate = neighbor_cell(index, 1, 0, shift_cur);
	if (cellelem_[nei_candidate] > 0 ){
		neivec.push_back(nei_candidate);
		shifts.push_back(shift_cur);
	}

	// The neighboring cell above
	nei_candidate = neighbor_cell(index, 0, 1, shift_cur);
	if (cellelem_[nei_candidate] > 0 ){
		neivec.push_back(nei_candidate);
		shifts.push_back(shift_cur);
	}

	// The neighboring cell to the top left and top right
	for (int lr = -1; lr < 2; lr+= 2){
		nei_candidate = neighbor_cell(index, lr, 1, shift_cur);
		if ( (cellelem_[nei_candidate] > 0) && (norm2(r - corner(index,lr,1)) < cutoffsquared )){
			neivec.push_back(nei_candidate);
			shifts.push_back(shift_cur);
		}
	}

	return neivec;
}


void partition::print() const {
  int index=0;
  for(int i=0; i < cellnum_; i++){
     std::cout << "------- CELL " << i << " -------" << std::endl;
     if (cellelem_[i] == 0){
        std::cout << "    EMPTY    " << std::endl;
     } else {
        std::cout << "cellelem_[" << i << "]=" << cellelem_[i] << std::endl;
        for (int k=0; k < cellelem_[i]; k++){
           std::cout << index << "   " << cellvec_[index] << std::endl;
           index++;
        }
     }
  std::cout << "---------------------" << std::endl;
  }
}


std::vector<int> partition::part_in_cell(int index) const {
   int first = find_first(index);
   // Returns the vector starting at cellvec_[first], ending et cellvec[first + cellelem_[index]].
   return std::vector<int>(cellvec_.begin() + first, cellvec_.begin() + first + cellelem_[index]);
}



std::vector<int> partition::nb_in_cell_index(int index, const topology::Vector2d& r,
		double cutoffsquared, const std::vector<topology::Vector2d >& positions,
		std::vector<double>& distances, const topology::Vector2d shift) const {
	distances.clear();
	std::vector<int> neighborvec;
	double dist_squared;
	int first = find_first(index);
	for (int i=0; i < cellelem_[index]; i++){
		dist_squared = norm2(positions[cellvec_[first+i]] - r + shift);
		if (dist_squared < cutoffsquared){
			neighborvec.push_back(cellvec_[first+i]);
			distances.push_back(std::sqrt(dist_squared));
		}
	}
	return neighborvec;
}


std::vector<int> partition::nb_in_cell_index_above(int index, const topology::Vector2d& r,
		double cutoffsquared, const std::vector<topology::Vector2d >& positions,
		std::vector<double>& distances, const topology::Vector2d shift) const {
	distances.clear();
   std::vector<int> neighborvec;
   topology::Vector2d dist;
   double dist_squared;
   int first = find_first(index);
   for (int i=0; i < cellelem_[index]; i++){
		dist = positions[cellvec_[first+i]] - r + shift;
		if (dist.get_y() > 0 || ( dist.get_y() == 0 && dist.get_x() > 0 )) {
			dist_squared = norm2(dist);
			if (dist_squared < cutoffsquared){
				neighborvec.push_back(cellvec_[first+i]);
				distances.push_back(std::sqrt(dist_squared));
			}
      }
   }
   return neighborvec;
}


//==================================================================================================================================
// Cluster analysis functions.
// cluster_analysis return a vector of cluster sizes and the corresponding vector of the number of particles (NoP) per cluster.
// The variable rho_min determines the minimum number a cell has to contain to be considered part of a cluster.

void partition::cluster_analysis(std::vector<int>& cluster_sizes, std::vector<int>& cluster_NoP, const int& rho_min) const {
  cluster_sizes.clear();
  cluster_NoP.clear();
  int NoP_counter;
  std::vector<int> current_elements, all_elements;
  for (int i=0; i < cellnum_; i++){
    if (cellelem_[i] >= rho_min && std::find(all_elements.begin(), all_elements.end(), i) == all_elements.end()){
      current_elements.push_back(i);
      current_elements=cluster_recursion(i, current_elements, rho_min);
      cluster_sizes.push_back(current_elements.size());
      NoP_counter=0;
      for (size_t j=0; j<current_elements.size(); j++) {
        NoP_counter += cellelem_[current_elements[j]];
        all_elements.push_back(current_elements[j]);
      }
      cluster_NoP.push_back(NoP_counter);
      current_elements.clear();
    }
  }
}

// recursive function cluster_recursion, return vector containing all the indices being part of the cluster.
std::vector<int> partition::cluster_recursion(int current_index, std::vector<int> current_elements, const int& rho_min) const {
  topology::Vector2d shift;
  int index;
  for (int i=-1; i < 2; i++){
    for (int j=-1; j < 2; j++){
      if (i != 0 || j != 0){
        index=neighbor_cell(current_index, i, j);
        if (cellelem_[index] >= rho_min && std::find(current_elements.begin(), current_elements.end(), index) == current_elements.end()){
//          std::cout << "index " << index << std::endl;
          current_elements.push_back(index);
          current_elements=cluster_recursion(index, current_elements, rho_min);
        }
      }
    }
  } 
  return current_elements;
}
