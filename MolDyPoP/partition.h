/*! \file partition.h
 *
 *  \brief Header-file to class declaration of partition. Introduces the partition,
 *  a cell list data structure that greatly facilitates neighbor calculation.
 *
 *  The partition is the preferred method for neighborhood calculation. The class neighbor_list
 *  is an alternative, but not recommended when neighborhoods change dynamically in each time
 *  step.
 *
 *  \author Thomas Bissinger

 *  \date Created: late 2017
 *  \date Last Updated: 2023-08-06
 *
 *  \class partition.
 *  \brief Defines the partition class. The simulation box is partitioned into cells. The indices of a 
 *  vector of particles (used for initialization) are sorted according to the cell they belong to.
 *  Has functions for printing and computing average velocities in a neighborhood.
 *  
 *  \author Thomas Bissinger
 *  
 *  \note A previous version of this class used pointers to the actual particles instead of their
 *  indices. This was changed due to memory access problems. The code is
 *  also simpler. On the other hand, one has to store the particle indices and combine the stored
 *  indices with the actual vector to access elements. Therefore, some functions need the original
 *  particle vector as a call-by-reference input.
 *
 *  \date Created: late 2017, presumably
 *  \date Last Updated: 2023-08-06
 *
 */


#ifndef PARTITION_H
#define PARTITION_H
#include "topology.h"
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>


class partition {
   protected:

      int N_; ///< Size of the group.

      /*! \brief Size of simulation box.
       *
       *  Box has to be at least twice the size of the cutoff radius, else particles can have
       *  double-interaction. This case is not an error caught by the simulation so far.
       *  Note that the simulation box is from 0 to L_.
       */
      topology::Vector2d L_;

      topology::Vector2d cellwidth_; ///< Width of a cell in the cell list. Typically the cutoff radius of the interaction.

      
      std::vector<int> M_ = std::vector<int>(2); ///< Number of cells per dimension. Dim = d = 2.
      int cellnum_; ///< Total number of cells.
      /*! \brief Indices of the first element pertaining to a cell. Dim = N_;
	   *
	   *  Set to the first element of the next cell if the cell is empty.
	   */
      std::vector<int> firsts_;

                                 
      /// Integers indicating the number of elements per cell. Dim = M_[0] * ... * M_[d];
      std::vector<int> cellelem_; 

      /// Vector of particle indices. Dim = N_. Easier to access than the list, and since it has to be accessed a lot, we save it here.
      std::vector<int> cellvec_;

   public:
		//   ===================================================================================================
		//   Constructors, destructors.
		//   ===================================================================================================
		/// Empty constructor.
		partition() : N_(0), L_(0.0), cellnum_(0) { clear(); } ;


		/*! \brief Standard constructor.
		*
		*  @param[in] N Number of particles.
		*  @param[in] cutoff Cutoff radius and width of cell.
		*  @param[in] L Size of simulation box, [0, L].
		*/
		partition(const int& N, const double& cutoff, const topology::Vector2d& L);

		/*! \brief Standard constructor.
		*
		*  @param[in] N Number of particles.
		*  @param[in] cutoff Cutoff radius and width of cell.
		*  @param[in] L Size of simulation box, [0, L].
		*  @param[in] positions Vector of positions to be sorted.
		*/
		partition(const int& N, const double& cutoff,
			  const topology::Vector2d& L,
			  const std::vector<topology::Vector2d>& positions);


		~partition(); ///< Destructor.
		void clear(); ///< Clears all entries.

		/*! \brief Fills the partition.
		*
		*  Assigns each particle to a cell in the partition.
		*
		*  @param[in] positions Vector of positions to be sorted.
		*/
		void fill(const std::vector<topology::Vector2d>& positions);


		//   ===================================================================================================
		//   get-functions
		//   ===================================================================================================
		/// Returns int cellelem_[m]
		inline int get_cellelem(int m) const { return cellelem_[m];};
		/// Returns vector cellelem_
		inline std::vector<int> get_cellelem() const { return cellelem_;};
		/// Returns int cellnum_
		inline int get_cellnum() const { return cellnum_;};
		/// Returns vector cellvec
		inline std::vector<int> get_cellvec() const { return cellvec_;};


		//   ===================================================================================================
		//   print-function
		//   ===================================================================================================
		void print() const ; ///< Print partition (only for troubleshooting).
    
		//   ===================================================================================================
		//   finding cells corresponding to a coordinate
		//   ===================================================================================================
		/*! \brief Returns the index of the cell containing r.
		 *
		 *  \warning Has no error output if r is not in the simulation box.
		 */
		int find_cell(const topology::Vector2d& r) const ;
		/// Returns index of the first element in cellvec_ pertaining to the cell with index index.
		inline int find_first(int index) const { return firsts_[index];};

		//   ===================================================================================================
		//   finding adjacent cells
		//   ===================================================================================================
		/*! \brief Returns index of neighboring cell to cell at index. Uses periodic boundary conditions.
		 *
		 *  lr and ur are used to describe where the neighboring cell is. E.g. a cell to the top left
		 *  would correspond to lr = -1 (left) and ud = 1 (up).
		 *  @param[in] index index of cell.
		 *  @param[in] lr neighboring cell position in left-right (x) direction. In {-1,0,1}. -1 is left.
		 *  @param[in] du neighboring cell position in down-up (y) direction. In {-1,0,1}. -1 is down.
		 */
		int neighbor_cell(int index, int lr, int du) const ;

		/*! \brief Returns index of neighboring cell to cell at index. Uses periodic boundary conditions. Contains a shift.
		 *
		 *  lr and ur are used to describe where the neighboring cell is. E.g. a cell to the top left
		 *  would correspond to lr = -1 (left) and ud = 1 (up).
		 *  @param[in] index index of cell.
		 *  @param[in] lr neighboring cell position in left-right (x) direction. In {-1,0,1} = (left, center, right).
		 *  @param[in] du neighboring cell position in down-up (y) direction. In {-1,0,1} = (below, center, above).
		 *  @param[inout] shift if the neighboring cell is on the other side of the periodic box, a nonzero shift is stored here.
		 *                The shift is designed such that the cell at index is next to its neighboring cell if shift is
		 *                subtracted from the positions of all particles in the cell at index. (Later function calls with r-shift)
		 */
		int neighbor_cell(int index, int lr, int du, topology::Vector2d& shift) const ;
		/*!  \brief Returns indices of all neighboring cells (including shifts)
		*
		*  These functions could be used if the old index is known. But the comparisons it
		*  saves are probably similarly expensive as the modulus calculations necessary to
		*  return to the original index and check whether the index is still accurate.
		*/
		std::vector<int> nb_cells_all(int index, std::vector<topology::Vector2d>& shifts) const;
		/// Returns indices of all neighboring cells (including shifts). Empty cells and cells out of reach of r are not included.
		std::vector<int> nb_cells_all(const topology::Vector2d& r, double cutoffsquared, std::vector<topology::Vector2d>& shifts) const;
		/// Returns indices the neighboring cells above and right of the cell and the cell itself (including shifts). Can be used to avoid double counting.
		std::vector<int> nb_cells_ur(int index, std::vector<topology::Vector2d>& shifts) const;
		/// Returns indices the neighboring cells above and right of the cell and the cell itself (including shifts). Empty cells and cells out of reach of r are not included. Can be used to avoid double counting.
		std::vector<int> nb_cells_ur(const topology::Vector2d& r, double cutoffsquared, std::vector<topology::Vector2d>& shifts) const;


		/// Get position of a corner of a cell. Which corner is determined by lr, ud.
		inline topology::Vector2d corner(int index, int lr, int ud) const { return topology::Vector2d(cellwidth_.get_x() * ((index % M_[0])+( 1 + lr )/2), cellwidth_.get_y() * (floor((double)(index)/M_[0])+( 1 + ud )/2)); };

		/// Returns a std::vector containing the indices of the particles in the box with index 'index'.
		std::vector<int> part_in_cell(int index) const ;

		/*! \brief Returns the indices of the particles in a cell for all particles within a cutoff radius of r.
		 *
		 *  @param[in] index Cell index
		 *  @param[in] r Position of particle.
		 *  @param[in] cutoffsquared Square of the cutoff radius (square because no roots have to be taken).
		 *  @param[in] positions Vector of positions to be accessed.
		 *  @param[out] distances Stores the distances to all neighbors. Saves computation time.
		 *  @param[in] shift Shifts position of vector r. Optional, default is 0.
		 *
		 *  ATTENTION: The integer neighbors is not re-initialized, so the number of neighbors is added to the
		 *             value that is already stored in neighbors. This simplifies calculations with more than one
		 *             group, but one has to bear that in mind.
		 */
		std::vector<int> nb_in_cell_index(int index, const topology::Vector2d& r, double cutoffsquared,
			  const std::vector<topology::Vector2d >& positions, std::vector<double>& distances,
			  const topology::Vector2d shift = 0) const ;

		std::vector<int> nb_in_cell_index_above(int index, const topology::Vector2d& r, double cutoffsquared,
			  const std::vector<topology::Vector2d >& positions, std::vector<double>& distances,
			  const topology::Vector2d shift = 0) const ; ///< Calculates only neighbors above particle (reduces computation time).

		//   ===================================================================================================
		//   Cluster analysis functions (haven't been used in a while,)
		//   ===================================================================================================
		/*! \brief Returns a vector of cluster sizes and the corresponding vector of the number of
		*  particles (NoP) per cluster. The variable rho_min determines the minimum number
		*  a cell has to contain to be considered part of a cluster.
		*
		*  @param[in] cluster_sizes   Stores size (number of pertaining cells) of each cluster. std::vector<int>
		*  @param[in] cluster_NoP     Stores the number of particles of each cluster. std::vector<int>
		*  @param[in] rho_min         Minimal density of cluster, i.e. minimal number of particles in a cell to be
		*                             considered part of the cluster. int
		*
		*  \note The function was originally implemented to check for large number fluctuations in Visek
		*  swarms. Hasn't been used since preliminary tests in the early 2018s. Might be useful to someone,
		*  but not recommended for use by the author.
		*/
		void cluster_analysis(std::vector<int>& cluster_sizes, std::vector<int>& cluster_NoP, const int& rho_min) const ;
		///< Analyze particle clusters contained in the partition.

		/*! \brief Recursive function, returns vector containing all the indices being part of a cluster.
		*  In principle, current_index is not necessary, one could also use the last entry of indices. But it doesn't disturb
		*  performance and is a little more readable.
		*
		*  @param[in] current_index  Index that is currently being viewed. int
		*  @param[in] indices        Indices that are already part of the cluster. std::vector<int>
		*  @param[in] rho_min        Minimal density of cluster, i.e. minimal number of particles in a cell to be
		*                            considered part of the cluster. int
		*/
		std::vector<int> cluster_recursion(int current_index, std::vector<int> indices, const int& rho_min) const ;
		///< Recursive function, finds all cells being part of a cluster.

		//==================================================================================================================================


} ;


#endif
