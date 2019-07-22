/**
 * @brief find the most optimal domain decomposition 
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */
#include <vector>
#include <map>
#include <numeric>
#include <limits>
#include <functional>
#include "MultiArrayIter.h"

#ifndef CUBE_DECOMP_H
#define CUBE_DECOMP_H

class CubeDecomp {

public:

  /**
   * No argument constructor
   */
  CubeDecomp() {}

  /**
   * Constructor
   * @param nprocs number of processes (communicator size)
   * @param dims global dimensions 
   */
  CubeDecomp(int nprocs, const std::vector<size_t>& dims) {

    this->build(nprocs, dims);
  }

  /** 
   * Build the decomposition
   * @param nprocs number of processes (communicator size)
   * @param dims global dimensions 
   * @return true if a valid decomposition was found, false otherwise
   */
  bool build(int nprocs, const std::vector<size_t>& dims);

/** 
 * Compute the decomposition that minimizes the ratio of surface to sub-domain volume
 * @return decomposition
 * @note returned array can be empty in which case no valid decomposition was found
 */
  std::vector<size_t> computeOptimalDecomp();

/**
 * Get the beginning index set for given rank
 * @param rk MPI rank
 * @return index set
 */
	std::vector<size_t> getBegIndices(int rk);

/**
 * Get the ending index set for given rank
 * @param rk MPI rank
 * @return index set
 */
  std::vector<size_t> getEndIndices(int rk);

/**
 * Get the rank for the sub-domain in given direction
 * @param rk base rank
 * @param direction away from the above base rank
 * @return rank
 * @note periodic boundary conditions are applied to ensure that
 *       the returned rank is valid
 */
  int getNeighborRank(int rk, const std::vector<int>& dir);
/**
 * Get the decomposition
 * @return domain decomposition
 */
  std::vector<size_t> getDecomp() const;

private:

  int nprocs;
  size_t ndims;
  std::vector<size_t> decomp;
  std::vector<size_t> dims;
  std::vector<size_t> localDims;
  MultiArrayIter mit;

/**
 * Get all the prime factors of given integer
 * @param n integer
 * @return list
 */
  std::vector<size_t> getPrimeFactors(size_t n) const;

/** 
 * Compute the cost (ratio of surface to volume)
 * @param decomp input decomposition
 * @return cost
 */
  double computeCost(const std::vector<size_t>& decomp) const;

};

#endif // CUBE_DECOMP_H
