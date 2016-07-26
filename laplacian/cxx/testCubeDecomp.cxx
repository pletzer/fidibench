/**
 * @brief test driver for CubeDecomp class
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include <vector>
#include <iostream>
#include "CubeDecomp.h"

void test(int sz, const std::vector<size_t>& dims) {

  CubeDecomp dc(sz, dims);
  std::vector<size_t> decomp = dc.getDecomp();
  std::cout << " [" << sz << "] decomp: ";
  for (size_t i = 0; i < decomp.size(); ++i) 
    std::cout << decomp[i] << ' ';
  std::cout << '\n';

  for (int rk = 0; rk < sz; ++rk) {

    std::cout << "\trank: " << rk << '\n';

    std::vector<size_t> lo = dc.getBegIndices(rk);
    std::vector<size_t> hi = dc.getEndIndices(rk);
    std::cout << "\tlo = ";
    for (size_t i = 0; i < lo.size(); ++i) 
      std::cout << lo[i] << ' ';
    std::cout << '\n';
    std::cout << "\thi = ";
    for (size_t i = 0; i < hi.size(); ++i) 
      std::cout << hi[i] << ' ';
    std::cout << '\n';

    std::vector<int> dir(dims.size(), 0);
    for (size_t dim = 0; dim < dims.size(); ++dim) {
      for (int pm = -1; pm <= 1; pm += 2) {
        dir[dim] = pm;
        std::cout << "\t\tneighbor in direction ";
        for (size_t i = 0; i < dir.size(); ++i) 
          std::cout << dir[i] << ' ';

        int neighborRk = dc.getNeighborRank(rk, dir);
        std::cout << " is " << neighborRk << '\n';

        // reset
        dir[dim] = 0;
      }
    }
  }
}

int main() {

  int nprocs;
  std::vector<size_t> dims;

// 1d test
  dims.resize(1);
  nprocs = 2;
  dims[0] = 10;
  test(nprocs, dims);

// 2d test
  dims.resize(2);
  nprocs = 6;
  dims[0] = 3; dims[1] = 2;
  test(nprocs, dims);

// 3d test 
  dims.resize(3);
  nprocs = 6;
  dims[0] = 2; dims[1] = 3; dims[2] = 4;
  test(nprocs, dims);
}

