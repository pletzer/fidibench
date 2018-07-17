#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include <limits>
#include <iostream>
#include <sstream>
#include <cstring>
#include "CmdLineArgParser.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

template <size_t NDIMS>
class Upwind {
public:
  // ctor
  Upwind(const std::vector<double>& velocity, 
         const std::vector<double>& lengths, 
         const std::vector<size_t>& numCells) {

    this->numCells = numCells;
    this->deltas.resize(NDIMS);
    this->upDirection.resize(NDIMS);
    this->v = velocity;
    this->lengths = lengths;
    this->ntot = 1;
    for (size_t j = 0; j < NDIMS; ++j) {
      this->upDirection[j] = -1;
      if (velocity[j] < 0.) this->upDirection[j] = +1;
      this->deltas[j] = lengths[j] / numCells[j];
      this->ntot *= numCells[j];
    }
    this->dimProd.resize(NDIMS);
    this->dimProd[NDIMS - 1] = 1;
    for (int i = (int) NDIMS - 2; i >= 0; --i) {
        // last index varies fastest
        this->dimProd[i] =  this->dimProd[i + 1] * this->numCells[i + 1];
      }
    this->f.resize(this->ntot, 0.0);
    this->oldF.resize(this->ntot, 0.0);
    // initialize lower corner to one
    this->f[0] = 1;
  }

  void advect(double deltaTime) {

    // copy
    this->oldF = this->f;
    const double* __restrict__ oldFPtr = &this->oldF[0];
    double* __restrict__ fPtr = &this->f[0];

#pragma omp parallel for
    for (int i = 0; i < (int) this->ntot; ++i) {

      std::vector<int> inds = this->getIndexSet(i);

      for (size_t j = 0; j < NDIMS; ++j) {

        int oldIndex = inds[j];
        const double coeff = deltaTime * this->v[j] * this->upDirection[j] / this->deltas[j];

        // periodic BCs
        inds[j] += this->upDirection[j] + this->numCells[j];
        inds[j] %= this->numCells[j];

        size_t upI = this->getFlatIndex(inds);

        fPtr[i] -= coeff * (oldFPtr[upI] - oldFPtr[i]);

        inds[j] = oldIndex;
      }
    }
  }

#include "saveVTK.h"


  double checksum() const {
    return std::accumulate(this->f.begin(), this->f.end(), 0.0, std::plus<double>());
  }

  void print() const {
    for (size_t i = 0; i < this->f.size(); ++i) {
      std::cout << i << " " << this->f[i] << '\n';
    }
  }

private:
  std::vector<double> v;
  std::vector<double> lengths;
  std::vector<double> deltas;
  std::vector<double> f;
  std::vector<double> oldF;
  std::vector<int> upDirection;
  std::vector<size_t> dimProd;
  std::vector<size_t> numCells;
  size_t ntot;

  inline
  std::vector<int> getIndexSet(size_t flatIndex) const {
    std::vector<int> res(NDIMS);
    for (size_t i = 0; i < NDIMS; ++i) {
      res[i] = flatIndex / this->dimProd[i] % this->numCells[i];
    }
    return res;
  }

  inline 
  size_t getFlatIndex(const std::vector<int>& inds) const {
    return std::inner_product(this->dimProd.begin(), this->dimProd.end(), inds.begin(), 0);
  }
};

///////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

  const int ndims = 3;

  CmdLineArgParser args;
  args.setPurpose("Purpose: benchmark finite difference operations.");
  args.set("-numCells", 128, "Number of cells along each axis");
  args.set("-numSteps", 10, "Number of time steps");
  args.set("-vtk", false, "Write output to VTK file");

  bool success = args.parse(argc, argv);
  bool help = args.get<bool>("-h");

  if (success && !help) {

    int numTimeSteps = args.get<int>("-numSteps");
    bool doVtk = args.get<bool>("-vtk");

#pragma omp parallel
    {
      int numThreads = 1;
      int maxNumThreads = 1;
      int threadId = 0;
#ifdef HAVE_OPENMP
      numThreads = omp_get_num_threads();
      maxNumThreads = omp_get_max_threads();
      threadId = omp_get_thread_num();
      if (threadId == 0)
        std::cout << "Running with OpenMP enabled\n";
#endif
      if (threadId == 0)
        std::cout << "number of threads: " << 
           numThreads << " max number of threads: " << maxNumThreads << '\n';
    }

    // same resolution in each direction
    std::vector<size_t> numCells(ndims, args.get<int>("-numCells"));
    std::cout << "number of cells: ";
    for (size_t i = 0; i < numCells.size(); ++i) {
      std::cout << ' ' << numCells[i];
    }
    std::cout << '\n';
    std::cout << "number of time steps: " << numTimeSteps << '\n';

    std::vector<double> velocity(numCells.size(), 1.0);
    std::vector<double> lengths(numCells.size(), 1.0);

    // compute dt 
    double courant = 0.1;
    double dt = std::numeric_limits<double>::max();
    for (size_t j = 0; j < velocity.size(); ++j) {
      double dx = lengths[j]/numCells[j];
      double val = courant * dx / velocity[j];
      dt = (val < dt? val: dt);
    }

    Upwind<ndims> up(velocity, lengths, numCells);
    if (doVtk) {
      up.saveVTK("up0.vtk");
    }
    for (int i = 0; i < numTimeSteps; ++i) {
      up.advect(dt);
    }
    std::cout << "check sum: " << up.checksum() << '\n';
    if (doVtk) {
      up.saveVTK("up1.vtk");
    }
  }
  else {
    // error when parsing command line arguments
    if (!success) {
      std::cerr << "ERROR when parsing command line arguments\n";
    }
    args.help();
  }


  return 0;
}
