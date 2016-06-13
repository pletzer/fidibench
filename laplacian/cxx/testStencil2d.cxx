/**
 * @brief test driver for general stencil finite differencing
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include <ostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <string>

#include "Filter.h"
#include "CmdLineArgParser.h"

/**
 * Input function
 */
double func(const std::vector<double>& pos) {
  double res = 1;
  for (size_t i = 0; i < pos.size(); ++i) {
    res *= sin(2.0 * M_PI * pos[i]);
  }
  return res;
}

double func(const std::vector<size_t>& inds) {
    // zero everywhere except in first cell
    size_t prod = 1;
    for (size_t i = 0; i < inds.size(); ++i) prod *= inds[i];
    if (prod == 0) {
        return 1.0;
    }
    return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {


  MPI_Init(&argc, &argv);

  int rk;
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);

  CmdLineArgParser args;
  args.setPurpose("Purpose: benchmark finite difference operations.");
  args.set("-numCells", 8, "Number of cells along each axis");
  args.set("-vtk", false, "Write output to VTK file");

  bool success = args.parse(argc, argv);
  bool help = args.get<bool>("-h");

  if (success && !help) {

    size_t numCells = (size_t) args.get<int>("-numCells");
    size_t numDims = 2;
    bool writeVTK = args.get<bool>("-vtk");

    // assemble the stencil
    std::map< std::vector<int>, double > stencil;
    std::vector<int> offset(numDims, 0);
    // diagonal
    stencil.insert( std::pair< std::vector<int>, double >(offset, 0.0) );
    // in the x direction
    offset[0] = 1;
    stencil.insert( std::pair< std::vector<int>, double >(offset, 1.0) );
    offset[0] = 0;
    // in the y direction
    offset[1] = -1;
    stencil.insert( std::pair< std::vector<int>, double >(offset, -1.0) );
    offset[1] = 0;

    // set the global dimensions
    std::vector<size_t> globalDims(numDims, numCells);

    // set the domain boundaries
    std::vector<double> xmins(numDims, 0.0);
    std::vector<double> xmaxs(numDims, 1.0);

    // create filter object
    Filter fltr(globalDims, xmins, xmaxs, stencil);

    if (fltr.getRank() == 0 && !fltr.isDecompValid()) {
      std::cerr << "Decomposition is invalid\n";
    }
    if (fltr.isDecompValid()) {
        fltr.setInDataByIndices(func);
        fltr.applyFilter();
    }
    
    //fltr.printInData();
    fltr.printOutData();

    if (writeVTK) {
        // write data to VTK file
        if (fltr.getRank() == 0) {
          std::cout << "Data will be written to file stencil2d.vtk\n";
        }
        fltr.saveVTK("stencil2d.vtk");
    }

  }
  else {
      // error when parsing command line arguments
      if (!success && rk == 0) {
          std::cerr << "ERROR when parsing command line arguments\n";
      }
      args.help();
  }

  MPI_Finalize();
  return 0;
}

