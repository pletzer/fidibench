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
    const std::map< std::vector<int>, double >& stencil) {

    this->validDecomp = false;
    this->ndims = globalDims.size();

    this->globalDims = globalDims;
    this->xmins = xmins;
    this->xmaxs = xmaxs;

    MPI_Comm_rank(MPI_COMM_WORLD, &this->rk);
    MPI_Comm_size(MPI_COMM_WORLD, &this->sz);

    // domain decomp
    this->validDecomp = this->dcmp.build(this->sz, globalDims);
    if (!this->validDecomp) {
      if (this->rk == 0) {
        std::cerr << "ERROR: No valid domain decomposition could be found. Adjust the number\n";
        std::cerr << "of processes and/or the domain dimensions.\n";
      }
      return;
    }

    this->lo = this->dcmp.getBegIndices(this->rk);
    this->hi = this->dcmp.getEndIndices(this->rk);
    std::vector<size_t> decomp = this->dcmp.getDecomp();
    if (this->rk == 0) {
      std::cout << "Number of procs: " << this->sz << "\nglobal dimensions: ";
      for (size_t i = 0; i < this->ndims; ++i)
        std::cout << globalDims[i] << ' ';
      std::cout << "\nDomain decomp ";
      for (size_t i = 0; i < this->ndims; ++i)
        std::cout << decomp[i] << ' ';
      std::cout << '\n';
    }

    // build the iterator over the entire sub-domain
    this->mit.build(this->lo, this->hi, false);
    size_t n = this->mit.getNumberOfTerms();
    this->inData.resize(n);
    this->outData.resize(n);

    // save the stencil
    this->stencil = stencil;

    // create all the windows
    for (std::map< std::vector<int>, double >::const_iterator it = stencil.begin(); 
      it != stencil.end(); ++it) {

      const std::vector<int>& offset = it->first;
      if (this->getNumberOfNonZeroElements(offset) == 0) {
        continue;
      }

      // side is the negative of offset
      std::vector<int> side(this->ndims);
      for (size_t i = 0; i < this->ndims; ++i) {
        side[i] = -offset[i];
      }

      this->newWindow(side);
    }
  }

  /**
   * Destructor
   */
  ~Filter() {

    // destroy windows
    for (std::map< std::vector<int>, MPI_Win >::const_iterator it = this->windows.begin(); 
      it != this->windows.end(); ++it) {
      this->deleteWindow(it->first);
    }
  }

/**
 * Get the position of a cell
 * @param globalInds global index set
 * @return position
 */
  std::vector<double> getPosition(const std::vector<size_t>& globalInds) {
    std::vector<double> pos(this->ndims);
    for (size_t i = 0; i < this->ndims; ++i) {
      double delta = (this->xmaxs[i] - this->xmins[i])/ double(this->globalDims[i]);
      // cell centered
      pos[i] = this->xmins[i] + (globalInds[i] + 0.5) * delta;
    }
    return pos;
  }

/**
 * Get remote data and store the result in dstData
 * @param otherRank rank from which the data will be fetched
 * @param side window side as defined with respect to the data source window
 */
 double* fetch(int otherRank, const std::vector<int>& side) {

  int n = (int) this->winIter.find(side)->second.getNumberOfTerms();
  MPI_Win win = this->windows.find(side)->second;

  double* dstData = this->winData.find(side)->second.second;
    
  MPI_Win_fence((MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE), win);
      
  MPI_Get(dstData, n, MPI_DOUBLE, otherRank, 0, n, MPI_DOUBLE, win);

  MPI_Win_fence(MPI_MODE_NOSUCCEED, win); 

  return dstData;
 }

/**
 * Set input data
 * @param f function pointer taking a position and returning a value
 */
  void setInData( double (*f)(const std::vector<double>&) ) {

// set the data in the main sub-domain
    this->mit.begin();
    for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {
      std::vector<size_t> inds = this->mit.getIndices();
      std::vector<double> pos = this->getPosition(inds);
      this->inData[i] = f(pos);
      this->mit.next();
    }

// set the data in the windows
    for (std::map< std::vector<int>, std::pair<double*, double*> >::const_iterator it = this->winData.begin(); 
      it != winData.end(); ++it) {

      const std::vector<int>& side = it->first;
      double* srcData = this->winData.find(side)->second.first;
      MultiArrayIter& wit = this->winIter.find(side)->second;

      wit.begin();
      for (size_t i = 0; i < wit.getNumberOfTerms(); ++i) {
        std::vector<size_t> inds = wit.getIndices();
        size_t bi = this->mit.computeBigIndex(inds);
        srcData[i] = this->inData[bi];
        wit.next();
      }
    }
  }

  /**
 * Set input data by indices
 * @param f function pointer taking a global index set and returning a value
 */
  void setInDataByIndices( double (*f)(const std::vector<size_t>&) ) {

// set the data in the main sub-domain
    this->mit.begin();
    for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {
      std::vector<size_t> inds = this->mit.getIndices();
      this->inData[i] = f(inds);
      this->mit.next();
    }

// set the data in the windows
    for (std::map< std::vector<int>, std::pair<double*, double*> >::const_iterator it = this->winData.begin(); 
      it != winData.end(); ++it) {

      const std::vector<int>& side = it->first;
      double* srcData = this->winData.find(side)->second.first;
      MultiArrayIter& wit = this->winIter.find(side)->second;

      wit.begin();
      for (size_t i = 0; i < wit.getNumberOfTerms(); ++i) {
        std::vector<size_t> inds = wit.getIndices();
        size_t bi = this->mit.computeBigIndex(inds);
        srcData[i] = this->inData[bi];
        wit.next();
      }
    }
  }


/** 
 * Apply filter to input data
 */
  void applyFilter() {

    // initialize output data to zero
    this->mit.begin();
    for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {
      this->outData[i] = 0;
      this->mit.next();
    }

    // send data from the remote src to the local dst container
    for (std::map< std::vector<int>, MPI_Win >::const_iterator it = this->windows.begin(); 
      it != this->windows.end(); ++it) {

      const std::vector<int>& side = it->first;
      MPI_Win win = it->second;

      MultiArrayIter& wit = this->winIter.find(side)->second;
      int n = (int) wit.getNumberOfTerms();

      std::vector<int> offset(this->ndims);
      for (size_t i = 0; i < this->ndims; ++i) {
        offset[i] = -side[i];
      }

      // apply periodic BCs, sub-domains wrap around to find the neighbor rank
      int neighborRk = this->dcmp.getNeighborRank(this->rk, offset);

      // fetch the remote data
      this->fetch(neighborRk, side);

    }

    // the ghosts have been sent to the local process. We're ready to apply the filter
    for (std::map< std::vector<int>, double >::const_iterator it = this->stencil.begin(); 
      it != this->stencil.end(); ++it) {

      const std::vector<int>& offset = it->first;
      double val = it->second;

      // window side is the negative of offset
      std::vector<int> side(this->ndims);
      for (size_t i = 0; i < this->ndims; ++i) {
        side[i] = -offset[i];
      }

      double* dstData = this->winData.find(side)->second.second;
      MultiArrayIter& wit = this->winIter.find(side)->second;

      // iterate over all the elements of the window
      this->mit.begin();
      for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {

        std::vector<size_t> inds = this->mit.getIndices();
        std::vector<int> indOffset(this->ndims);

        // apply the stencil offset & check that we are still in the 
        // sub-domain
        bool insideLocalDomain = true;
        for (size_t j = 0; j < this->ndims; ++j) {
          indOffset[j] = (int) inds[j] + offset[j];
          insideLocalDomain &= (indOffset[j] >= (int) this->lo[j]);
          insideLocalDomain &= (indOffset[j] < (int) this->hi[j]);
        }

        if (insideLocalDomain) {
          // local update
          int bi2 = this->mit.computeBigIndex(indOffset);
          this->outData[i] += val * this->inData[ (size_t) bi2];
        }
        else {
          // update using halo data
          size_t bi2 = wit.computeBigIndex(inds);
          this->outData[i] += val * dstData[bi2];
        }

        this->mit.next();
      }
    }

  }

/** 
 * Get the rank on this process
 * @return rank
 */
  int getRank() const {
    return this->rk;
  }

/**
 * Get the number of processes (communicator size)
 * @return number
 */
  int getNumProcs() const {
    return this->sz;
  }

/**
 * Create and allocate data structures associated with an MPI window
 * @param side array containing 0s, -1s, and +1s used to uniquely denote a window
 * @param numGhosts number of ghosts
 */
  void newWindow(const std::vector<int>& side, size_t numGhosts = 1) {

    // compute the size of the window
    std::vector<size_t> sLo = this->lo;
    std::vector<size_t> sHi = this->hi;

    size_t n = 1;
    for (size_t i = 0; i < this->ndims; ++i) {
      if (side[i] != 0) {
        n *= numGhosts;
        if (side[i] < 0) {
          // low side
          sHi[i] = sLo[i] + numGhosts;
        }
        else {
          // hi side
          sLo[i] = sHi[i] - numGhosts;
        }
      }
      else {
        n *= (this->hi[i] - this->lo[i]);
      }
    }

    // allocate data and create window
    double* srcData = 0;
    double* dstData = 0;
    MPI_Alloc_mem(sizeof(double)*n, MPI_INFO_NULL, &srcData);
    MPI_Alloc_mem(sizeof(double)*n, MPI_INFO_NULL, &dstData);

    MPI_Win win;
    MPI_Win_create(srcData, n*sizeof(double), sizeof(double), 
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    MultiArrayIter sit(sLo, sHi, false);
    std::pair<double*, double*> pSrcDst(srcData, dstData);
    this->winData.insert( std::pair< std::vector<int>, std::pair<double*, double*> >(side, pSrcDst) );
    this->windows.insert( std::pair< std::vector<int>, MPI_Win >(side, win) );
    this->winIter.insert( std::pair< std::vector<int>, MultiArrayIter >(side, sit) );
  }

/** 
 * Delete data structures associated with an MPI window
 * @param side array containing 0s, -1s, and +1s used to uniquely denote a window
 */
  void deleteWindow(const std::vector<int>& side) {

    std::map< std::vector<int>, std::pair<double*, double*> >::iterator dataIt = this->winData.find(side);
    if (dataIt != this->winData.end()) {
      MPI_Free_mem(dataIt->second.first);
      MPI_Free_mem(dataIt->second.second);
    }
    
    std::map< std::vector<int>, MPI_Win >::iterator winIt = this->windows.find(side);
    if (winIt != this->windows.end()) {
      MPI_Win_free(&winIt->second);
    }

  }

/**
 * Display the content of a data window on the source rank
 * @param side window identifier
 */
  void printSrcDataWindow(const std::vector<int>& side) {

    double* data = this->winData.find(side)->second.first;
    MultiArrayIter& wit = this->winIter.find(side)->second;
    wit.begin();
    std::cerr << "[" << this->rk << "] srcData for side: ";
    for (size_t j = 0; j < this->ndims; ++j)
      std::cerr << side[j] << ' ';
    std::cerr << '\n';
    for (size_t i = 0; i < wit.getNumberOfTerms(); ++i) {

      std::vector<size_t> inds = wit.getIndices();

      std::cerr  << "\t[" << this->rk << "] inds = ";
      for (size_t j = 0; j < this->ndims; ++j) 
        std::cerr << inds[j] << ' ';
      std::cerr << " srcData = " << data[i] << '\n';

      wit.next();
    }
  }

/**
 * Display the content of a data window on the destination rank
 * @param side window identifier
 */
  void printDstDataWindow(const std::vector<int>& side) {

    double* data = this->winData.find(side)->second.second;
    MultiArrayIter& wit = this->winIter.find(side)->second;
    wit.begin();
    std::cerr << "[" << this->rk << "] dstData for side: \n";
    for (size_t j = 0; j < this->ndims; ++j)
      std::cerr << side[j] << ' ';
    std::cerr << '\n';
    for (size_t i = 0; i < wit.getNumberOfTerms(); ++i) {

      std::vector<size_t> inds = wit.getIndices();

      std::cerr  << "\t[" << this->rk << "] inds = ";
      for (size_t j = 0; j < this->ndims; ++j) 
        std::cerr << inds[j] << ' ';
      std::cerr << " dstData = " << data[i] << '\n';

      wit.next();
    }
  }

/**
 * Display all the input data on this process
 */
  void printInData() {

    this->mit.begin();
    std::cerr << "[" << this->rk << "] inData: \n";
    for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {

      std::vector<size_t> inds = this->mit.getIndices();

      std::cerr  << "\t[" << this->rk << "] inds = ";
      for (size_t j = 0; j < this->ndims; ++j) 
        std::cerr << inds[j] << ' ';
      std::cerr << " inData = " << this->inData[i] << '\n';

      this->mit.next();
    }
  }

/**
 * Display all the output data on this process
 */
  void printOutData() {

    this->mit.begin();
    std::cerr << "[" << this->rk << "] outData: \n";
    for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {

      std::vector<size_t> inds = this->mit.getIndices();

      std::cerr  << "\t[" << this->rk << "] inds = ";
      for (size_t j = 0; j < this->ndims; ++j) 
        std::cerr << inds[j] << ' ';
      std::cerr << " outData = " << this->outData[i] << '\n';

      this->mit.next();
    }
  }

/**
 * Copy the output data into the input data container
 */
  void copyOutToIn() {
    for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {
      this->inData[i] = this->outData[i];
    }
  }

/**
 * Compute check sums of the data
 * @param inOrOut either "input" or "output"
 * @return check sum
 */
  double computeCheckSum(const std::string& inOrOut) {

    // local sum then reduce
    double localSum = 0, globalSum = 0;
    this->mit.begin();
    if (inOrOut == "input") {
      for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {
        localSum += this->inData[i];
        this->mit.next();
      }
    }
    else {
      for (size_t i = 0; i < this->mit.getNumberOfTerms(); ++i) {
        localSum += this->outData[i];
        this->mit.next();
      }
    }
    MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return globalSum;
  }

/** 
 * Save outout data to VTK file
 * @param filename file name
 */
  void saveVTK(const std::string& filename) {

    size_t ntot = 1;
    for (size_t i = 0; i < this->ndims; ++i) {
      ntot *= globalDims[i];
    }

    int nLocal = (int) this->outData.size();
    int n = (int) ntot;

    std::vector<double> recvBuffer;
    if (this->rk == 0) {
      recvBuffer.resize(ntot);
    }

    // gather the data
    MPI_Gather(&this->outData.front(), nLocal, MPI_DOUBLE,
               &recvBuffer.front(), nLocal, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    std::vector<double> globalField;
    std::vector<size_t> globalInds(this->ndims);

    if (this->rk == 0) {

      // re-assemble the field and write to VTK file
      globalField.resize(ntot); 

      std::vector<size_t> zeros(this->ndims, 0);
      MultiArrayIter globalIt(zeros, globalDims, false);

      for (int rk = 0; rk < this->sz; ++rk) {

        // get the global start indices for this chunk 
        std::vector<size_t> begInds = this->dcmp.getBegIndices(rk);

        // assume all sub-domains to have the same size!
        this->mit.begin();
        for (size_t iLocal = 0; iLocal < this->mit.getNumberOfTerms(); ++iLocal) {
          std::vector<size_t> localInds = this->mit.getIndices();
          for (size_t i = 0; i < this->ndims; ++i) {
            globalInds[i] = begInds[i] + localInds[i];
          }
          size_t globalBigIndex = globalIt.computeBigIndex(globalInds);
          globalField[globalBigIndex] = recvBuffer[rk*nLocal + iLocal];
          this->mit.next();
        }
      }
      // write 
      writeVTK(filename, this->globalDims, this->xmins, this->xmaxs, globalField, "outData");
    }
  }

/**
 * Id the decomposition valid?
 * @param return true if valid, false otherwise
 */
 bool isDecompValid() const {
  return this->validDecomp;
 }

private:

  size_t getNumberOfNonZeroElements(const std::vector<int>& v) {
    size_t res = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      if (v[i] != 0) {
        res++;
      }
    }
    return res;
  }

  // side to MPI window association
  std::map< std::vector<int>, MPI_Win > windows;

  // src/dst data for each sub-domain
  std::map< std::vector<int>, std::pair<double*, double*> > winData;

  // iterators for each sub-domain
  std::map< std::vector<int>, MultiArrayIter > winIter;

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
