/**
 * @brief Apply a differencing filter to input data distributed across multiple processes
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */
#include <ostream>
#include <cstdlib>
#include <mpi.h>
#include <vector>
#include <map>
#include <cmath>
#include <string>

#include <mpi.h>
#include "CubeDecomp.h"
#include "MultiArrayIter.h"
#include "writeVTK.h"

#ifndef FILTER_H
#define FILTER_H

class Filter {

public:

  /**
   * Constructor
   * @param globalDims global dimensions
   * @param xmins low domain corner
   * @param xmaxs high domain corner
   * @param stencil discretization stencil (offset -> value map)
   */
  Filter(const std::vector<size_t>& globalDims, 
    const std::vector<double>& xmins, 
    const std::vector<double>& xmaxs,
    const std::map< std::vector<int>, double >& stencil);

  /**
   * Destructor
   */
  ~Filter();
    
/**
 * Get the side from the offset
 * @param offset stencil offset vector from base cell
 * @return side 
 */
    std::vector<int> getSideFromOffset(const std::vector<int>& offset);

/**
 * Get the position of a cell
 * @param globalInds global index set
 * @return position
 */
  std::vector<double> getPosition(const std::vector<size_t>& globalInds);

/**
 * Get remote data and store the result in dstData
 * @param otherRank rank from which the data will be fetched
 * @param side window side as defined with respect to the data source window
 */
 double* fetch(int otherRank, const std::vector<int>& side);

/**
 * Set input data
 * @param f function pointer taking a position and returning a value
 */
  void setInData( double (*f)(const std::vector<double>&) );

/**
 * Set input data by indices
 * @param f function pointer taking a global index set and returning a value
 */
  void setInDataByIndices( double (*f)(const std::vector<size_t>&) );


/** 
 * Apply filter to input data
 */
  void applyFilter();
/** 
 * Get the rank on this process
 * @return rank
 */
  int getRank() const;

/**
 * Get the number of processes (communicator size)
 * @return number
 */
  int getNumProcs() const;

/** ;
 * Get the number of ghosts for an offset vector
 * @param offset displacement vector
 * @return number
 * @note assumes offset is along one of the principal axes (only one non-zero value)
 */
  size_t getNumberOfGhostsFromOffset(const std::vector<int>& offset) const;

/**
 * Get remote window interator 
 * @param rk MPI rank of the remote window
 * @param side window side
 * @return window interator
 * @note returns the domain iterator if side == {0, 0, 0,..}
 */
  MultiArrayIter getWindowIter(int rk, const std::vector<int>& side);

/**
 * Create and allocate data structures associated with an MPI window
 * @param side array containing 0s, -1s, and +1s used to uniquely denote a window
 * @param numGhosts number of ghosts
 */
  void newWindow(const std::vector<int>& side);

/** 
 * Delete data structures associated with an MPI window
 * @param side array containing 0s, -1s, and +1s used to uniquely denote a window
 */
  void deleteWindow(const std::vector<int>& side);

/**
 * Display the content of a data window on the source rank
 * @param side window identifier
 */
  void printSrcDataWindow(const std::vector<int>& side);

/**
 * Display the content of a data window on the destination rank
 * @param side window identifier
 */
  void printDstDataWindow(const std::vector<int>& side);

/**
 * Display all the input data on this process
 */
  void printInData();

/**
 * Display all the output data on this process
 */
  void printOutData();

/**
 * Copy the output data into the input data container
 */
  void copyOutToIn();

/**
 * Compute check sums of the data
 * @param inOrOut either "input" or "output"
 * @return check sum
 */
  double computeCheckSum(const std::string& inOrOut);

/** 
 * Save outout data to VTK file
 * @param filename file name
 */
  void saveVTK(const std::string& filename);

/**
 * Is the decomposition valid?
 * @param return true if valid, false otherwise
 */
 bool isDecompValid() const;

private:

  size_t getNumberOfNonZeroElements(const std::vector<int>& v);

  // side to MPI window association
  std::map< std::vector<int>, MPI_Win > windows;

  // side to src/dst data for each sub-domain
  std::map< std::vector<int>, std::pair<double*, double*> > winData;

  // side to iterator for the window
  std::map< std::vector<int>, MultiArrayIter > winIter;

  // side to the number of ghosts
  std::map< std::vector<int>, size_t > winNumGhosts;

  // stencil 
  std::map< std::vector<int>, double > stencil;

  // the data before applying the stencil
  std::vector<double> inData;

  // the data after applyng the stencil
  std::vector<double> outData;

  // the grid dimensions prior to the domain decomposition
  std::vector<size_t> globalDims;

  // the low/high end corner values of the domain
  std::vector<double> xmins;
  std::vector<double> xmaxs;

  // the low/one past the last local index sets
  std::vector<size_t> lo;
  std::vector<size_t> hi;

  // the domain decomposition object
  CubeDecomp dcmp;

  // local sub-domain cell iterator
  MultiArrayIter mit;

  // the MPI rank
  int rk;

  // number of processes
  int sz;

  // number of space dimensions
  size_t ndims;

  // true if a valid decomposition  was found
  bool validDecomp;

};

#endif // FILTER_H
