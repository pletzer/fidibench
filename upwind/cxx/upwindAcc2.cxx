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

template <size_t NDIMS>
class Upwind {
public:
  // ctor
  Upwind(const std::vector<double>& velocity, 
         const std::vector<double>& lengths, 
         const std::vector<int>& numCells) {

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
    // initialize lower corner to one
    this->f[0] = 1;

    this->fOld.resize(this->ntot);

    this->coeff.resize(NDIMS);
    for (size_t j = 0; j < NDIMS; ++j) {
      this->coeff[j] = this->v[j] * this->upDirection[j] / this->deltas[j];
    }
  }

  void advect(double deltaTime) {

    // OpenACC works with primitive arrays
    double* fPtr = &this->f.front();
    double* fOldPtr = &this->fOld.front();
    double* coeffPtr = &this->coeff.front();
    int* upDirectionPtr = &this->upDirection.front();
    int* dimProdPtr = &this->dimProd.front();
    int* numCellsPtr = &this->numCells.front();
    int ntot = this->ntot;
    int upI;

#pragma acc data \
  copy(fPtr[ntot]) \
  create(fOldPtr[ntot]) \
  copyin(coeffPtr[NDIMS], dimProdPtr[NDIMS],	\
  upDirectionPtr[NDIMS], numCellsPtr[NDIMS], deltaTime)
   {

#pragma acc parallel loop
    for (int i = 0; i < ntot; ++i) {
      fOldPtr[i] = fPtr[i];
    }

#pragma acc parallel loop
    for (int i = 0; i < ntot; ++i) {

      int inds[NDIMS];

#include "compute_index_set.h"

#include "compute_flat_index_offset_x.h"
fPtr[i] -= deltaTime * coeffPtr[0] * (fOldPtr[upI] - fOldPtr[i]);

#include "compute_flat_index_offset_y.h"
fPtr[i] -= deltaTime * coeffPtr[1] * (fOldPtr[upI] - fOldPtr[i]);

#include "compute_flat_index_offset_z.h"
fPtr[i] -= deltaTime * coeffPtr[2] * (fOldPtr[upI] - fOldPtr[i]);

     } // acc parallel loop
    } // acc data

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
  std::vector<double> f;
  std::vector<double> fOld;
  std::vector<double> v;
  std::vector<double> lengths;
  std::vector<double> deltas;
  std::vector<double> coeff;
  std::vector<int> upDirection;
  std::vector<int> dimProd;
  std::vector<int> numCells;
  int ntot;
};


int main(int argc, char** argv) {

  const int ndims = 3;
  CmdLineArgParser args;
  args.setPurpose("Purpose: benchmark finite difference operations.");
  args.set("-numCells", 128, "Number of cells along each axis");
  args.set("-numSteps", 10, "Number of time steps");
  args.set("-vtk", false, "Write output to VTK file");

  bool success = args.parse(argc, argv);
  bool help = args.get<bool>("-h");

  if (!success) {
    std::cerr << "ERROR: wrong command line arguments\n";
    args.help();
    return 1;
  } else if (help) {
    args.help();
    return 0;
  }

  int numCellsXYZ = args.get<int>("-numCells");
  int numTimeSteps = args.get<int>("-numSteps");
  bool doVtk = args.get<bool>("-vtk");

  // same resolution in each direction
  std::vector<int> numCells(ndims, numCellsXYZ);
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
  //up.saveVTK("up0.vtk");
  for (int i = 0; i < numTimeSteps; ++i) {
    up.advect(dt);
  }
  std::cout << "check sum: " << up.checksum() << '\n';
  if (doVtk) {
    up.saveVTK("up.vtk");
  }
}
