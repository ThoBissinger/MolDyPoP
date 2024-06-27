/*! \file neighbor_list.h
 *
 *  \brief Header-file to class declaration of neighbor_list. Introduces a data structure
 *  storing the neighbors to each particle. Useful when neighborhoods stay static.
 *
 *  Works especially well for disordered XY models and the like. Can also be used for
 *  lattice-based models, though the memory access for calls to the neighbor_list
 *  may be more expensive than the straightforward modulo calculations associated
 *  with neighbor calculations in lattice-based simulations
 *
 *
 *  \author Thomas Bissinger

 *  \date Created: mid 2019
 *  \date Last Updated: 2023-08-06
 *
 *
 *  \class neighbor_list
 *
 *  \brief Defines the neighbor_list class. Can be used to extract all information about
 *  neighborhood in a group, that is the neighbor pairs and all the distances.
 *  Neighbors are those other particles within the cutoff radius or the nearest
 *  neighbors in case of a lattice system.
 *  
 *  \author Thomas Bissinger
 *  
 *
 *  
 *  
 */



#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H
#include "topology.h"
#include "partition.h"
#include "group.h"

#include <fstream>
#include <list>
#include <vector>
#include <algorithm>


class group ;

class neighbor_list {
	protected:

		int N_ = 0; ///< Number of particles.

		std::vector<int> nb_indices_;		///< Indices of neighbors. Only efficient for systems with fixed neighbors
		std::vector<int> nb_first_; 		///< Index of first neighbor of each particle. Helpful when accessing nb_indices.
		std::vector<double> nb_dist_; 		///< Distances to neighbors. Only efficient for systems with fixed neighbors
	public:
		// ==============================================================================================
		// Constructors, destructors.
		// ==============================================================================================
		/// Empty constructor.
		neighbor_list() {} ;

		/*! \brief Standard constructor.
		 *
		 *  @param[in] group the group based on which the neighbor_list has to be intialized. WARNING: Off-lattice systems need a full partition.
		*/
		neighbor_list(const group& G);


		// ==============================================================================================
		// Clear function
		// ==============================================================================================
		/// Clears the neighbor_list
		void clear();

		// ==============================================================================================
		// Get functions.
		// ==============================================================================================
		/// Returns whether or not the neighbor_list is empty
		inline bool is_empty() const { return N_ == 0;} ;
		/// Returns all neighbors of particle i (elements of nb_indices_)
		std::vector<int> get_neighbors(int i) const ;
		/// Returns all distances between neighbors and particle i (elements of nb_dist_)
		std::vector<double> get_dist(int i) const ;



} ;

#endif
