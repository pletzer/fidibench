// clang -I/usr/local/openmp-llvm/include/ -O3 upwindTemplated.cxx -L/usr/local/openmp-llvm/lib/ -lstdc++ -liomp5
// g++ -fopenmp -O3 upwindTemplated.cxx

#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include <limits>
#include <iostream>
#include <sstream>

#include <omp.h>

template <size_t NDIMS>
class Upwind {
public:
  // ctor
  Upwind(const std::vector<double>& velocity, const std::vector<double>& lengths, const std::vector<size_t>& numCells) {
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
  }

  void advect(double deltaTime) {

    std::vector<double> oldF(this->f);

#pragma omp parallel for
    for (int i = 0; i < (int) this->ntot; ++i) {

      std::vector<int> inds = this->getIndexSet(i);

      for (size_t j = 0; j < NDIMS; ++j) {

        int oldIndex = inds[j];

        // periodic BCs
        inds[j] += this->upDirection[j] + this->numCells[j];
        inds[j] %= this->numCells[j];

        size_t upI = this->getFlatIndex(inds);

        this->f[i] -= deltaTime * this->v[j] * this->upDirection[j] * (oldF[upI] - oldF[i])/this->deltas[j];

        inds[j] = oldIndex;
      }
    }
  }

  void saveVTK(const std::string& filename) {
    std::fstream file;
    file.open(filename.c_str(), std::ios_base::out);
    file << "# vtk DataFile Version 2.0\n";
    file << "upwind.cxx\n";
    file << "ASCII\n";
    file << "DATASET RECTILINEAR_GRID\n";
    file << "DIMENSIONS";
    // in VTK the first dimension varies fastest so need 
    // to invert the order of the dimensions
    if (NDIMS > 2) {
      file << ' ' << this->numCells[2] + 1;
    }
    else {
      file << " 1";
    }
    if (NDIMS > 1) {
      file << ' ' << this->numCells[1] + 1;
    }
    else {
      file << " 1";
    }
    file << ' ' << this->numCells[0] + 1 << '\n';
    file << "X_COORDINATES ";
    if (NDIMS > 2) {
      file << this->numCells[2] + 1 << " double\n";
      for (size_t i = 0; i < this->numCells[2] + 1; ++i) {
        file << ' ' << 0.0 + this->deltas[2] * i;
      }      
    }
    else {
      file << "1 double\n";
      file << "0.0\n";
    }
    file << "Y_COORDINATES ";
    if (NDIMS > 1) {
      file << this->numCells[1] + 1 << " double\n";
      for (size_t i = 0; i < this->numCells[1] + 1; ++i) {
        file << ' ' << 0.0 + this->deltas[1] * i;
      }      
    }
    else {
      file << "1 double\n";
      file << "0.0\n";
    }
    file << "Z_COORDINATES ";
    file << this->numCells[0] + 1 << " double\n";
    for (size_t i = 0; i < this->numCells[0] + 1; ++i) {
      file << ' ' << 0.0 + this->deltas[0] * i;
    }
    file << "CELL_DATA " << this->ntot << '\n';
    file << "SCALARS f double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < this->ntot; ++i) {
      file << this->f[i] << '\n';
    }
    file.close();
  }

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
  std::vector<int> upDirection;
  std::vector<size_t> dimProd;
  std::vector<size_t> numCells;
  size_t ntot;

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


int main(int argc, char** argv) {

  const int ndims = 3;

  if (argc <= 1) {
    std::cout << "must specify number of cells in each direction.\n";
    return 1;
  }

  int numTimeSteps = 100;
  if (argc > 2) {
    numTimeSteps = atoi(argv[2]);
  }


#pragma omp parallel
  {
    int numThreads = omp_get_num_threads();
    int threadId = omp_get_thread_num();
    if (threadId == 0)
      std::cout << "number of threads: " << numThreads << '\n';
  }

  // same resolution in each direction
  std::vector<size_t> numCells(ndims, atoi(argv[1]));
  std::cout << "number of cells: ";
  for (size_t i = 0; i < numCells.size(); ++i) {
    std::cout << ' ' << numCells[i];
  }
  std::cout << '\n';

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
      //if (i % 10 == 0) {
      // std::ostringstream ss;
      // ss << i;
      // up.saveVTK(std::string("up") + ss.str() + std::string(".vtk"));
      //}
  }
  //up.print();
  std::cout << "check sum: " << up.checksum() << '\n';
  //up.saveVTK("up.vtk");
}
