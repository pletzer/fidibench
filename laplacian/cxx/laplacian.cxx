/**
 * @brief test driver for Laplacian-type differencing
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

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {


  MPI_Init(&argc, &argv);

  int rk;
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);

  CmdLineArgParser args;
  args.setPurpose("Purpose: benchmark finite difference operations.");
  args.set("-numCells", 8000, "Number of cells along each axis");
  args.set("-numDims", 2, "Number of dimensions");
  args.set("-vtk", false, "Write output to VTK file");

  bool success = args.parse(argc, argv);
  bool help = args.get<bool>("-h");

  if (success && !help) {

    size_t numCells = (size_t) args.get<int>("-numCells");
    size_t numDims = (size_t) args.get<int>("-numDims");
    bool writeVTK = args.get<bool>("-vtk");

    // assemble the stencil
    std::map< std::vector<int>, double > stencil;
    std::vector<int> offset(numDims, 0);
    // diagonal
    stencil.insert( std::pair< std::vector<int>, double >(offset, -2.0 * numDims) );
    for (size_t i = 0; i < numDims; ++i) {
      offset[i] = 1;
      stencil.insert( std::pair< std::vector<int>, double >(offset, 1.0) );
      offset[i] = -1;
      stencil.insert( std::pair< std::vector<int>, double >(offset, 1.0) );
      offset[i] = 0;
    }

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

      fltr.setInData(func);
      double tic = MPI_Wtime();

      // repeat to improve statistics
      size_t numIter = 10;
      for (size_t i = 0; i < numIter; ++i) {
        fltr.applyFilter();
        fltr.copyOutToIn();
      }

      double toc = MPI_Wtime();
      double walltime = toc - tic;
      double sumTime, minTime, maxTime;
      MPI_Reduce(&walltime, &sumTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&walltime, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&walltime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      if (writeVTK) {
        // write data to VTK file
        if (fltr.getRank() == 0) {
          std::cout << "Data will be written to file laplacian.vtk\n";
        }
        fltr.saveVTK("laplacian.vtk");
      }

      double inSum = fltr.computeCheckSum("input");
      double outSum = fltr.computeCheckSum("output");

      if (fltr.getRank() == 0) {
        std::cout << "Laplace times min/max/avg: " << minTime << '/' << maxTime << '/' << 
          sumTime/fltr.getNumProcs() << " [seconds]\n";
        std::cout << "Check sums: input = " << inSum << " output = " << outSum << '\n';
      }
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

