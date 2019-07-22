/**
 * @brief test driver for upwind-type differencing
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
#include <limits>
#include <functional>

#include "Filter.h"
#include "CmdLineArgParser.h"

double initialCondition(const std::vector<size_t>& inds) {
  double res = 0;
  if (std::accumulate(inds.begin(), inds.end(), 0) == 0) {
    res = 1;
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
  args.set("-numCells", 128, "Number of cells along each axis");
  args.set("-numSteps", 10, "Number of time steps");
  args.set("-vtk", false, "Write output to VTK file");

  bool success = args.parse(argc, argv);
  bool help = args.get<bool>("-h");

  if (success && !help) {

    // number of space dimensions
    const size_t numDims = 3;

    size_t numCells = (size_t) args.get<int>("-numCells");
    size_t numSteps = (size_t) args.get<int>("-numSteps");
    bool writeVTK = args.get<bool>("-vtk");

    // velocity field, constant
    const std::vector<double> velocities(numDims, 1.);
    // domain lengths
    const std::vector<double> lengths(numDims, 1.);
    // time step
    double courant = 0.1;
    std::vector<double> deltas(numDims);
    std::vector<int> signs(numDims);
    double dt = std::numeric_limits<double>::max();
    for (size_t j = 0; j < velocities.size(); ++j) {
      double dx = lengths[j] / (double) numCells;
      deltas[j] = dx;
      double val = courant * dx / velocities[j];
      dt = (val < dt? val: dt);
      signs[j] = (velocities[j] > 0? 1: -1);
    }

    // assemble the stencil
    std::map< std::vector<int>, double > stencil;
    std::vector<int> offset(numDims, 0);

    // diagonal
    double val = 1.0;
    for (size_t i = 0; i < velocities.size(); ++i) {
      val -= signs[i] * dt * velocities[i] / deltas[i];
    }
    stencil.insert( std::pair< std::vector<int>, double >(offset, val) );

    // upwind contributions
    for (size_t i = 0; i < numDims; ++i) {
      offset[i] = -signs[i];
      stencil.insert( std::pair< std::vector<int>, double >(offset, 
                      signs[i] * dt * velocities[i] / deltas[i]) );
      // reset
      offset[i] = 0;
    }

    // set the global dimensions
    std::vector<size_t> globalDims(numDims, numCells);

    // set the domain boundaries
    std::vector<double> xmins(numDims, 0.0);

    // create filter object
    Filter fltr(globalDims, xmins, lengths, stencil);

    if (fltr.getRank() == 0 && !fltr.isDecompValid()) {
      std::cerr << "Decomposition is invalid\n";
    }
    if (fltr.isDecompValid()) {

      double tic = MPI_Wtime();

      // set the initial condition
      fltr.setInDataByIndices(initialCondition);
      //fltr.printInData();

      // advect in time
      for (size_t i = 0; i < numSteps; ++i) {
        fltr.applyFilter();
#define CHECK_NAN
#ifdef CHECK_NAN
        double inSum = fltr.computeCheckSum("input");
        double outSum = fltr.computeCheckSum("output");
        if (fltr.getRank() == 0) {
        	std::cout << "iter " << i << " check sum  in/out = " << inSum << " / " << outSum << '\n';
        }
#endif
        fltr.copyOutToIn();
      }
      //fltr.printOutData();

      double toc = MPI_Wtime();
      double walltime = toc - tic;
      double sumTime, minTime, maxTime;
      MPI_Reduce(&walltime, &sumTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&walltime, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&walltime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      if (writeVTK) {
        // write data to VTK file
        if (fltr.getRank() == 0) {
          std::cout << "Data will be written to file upMpi.vtk\n";
        }
        fltr.saveVTK("upMpi.vtk");
      }

      double outSum = fltr.computeCheckSum("output");

      if (fltr.getRank() == 0) {
        std::cout << " times min/max/avg: " << minTime << '/' << maxTime << '/' << 
          sumTime/fltr.getNumProcs() << " [seconds]\n";
        std::cout << "Check sum: " << outSum << '\n';
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


