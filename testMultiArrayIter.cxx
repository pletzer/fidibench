/**
 * @brief test driver for a MultiArrayIter class
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include <vector>
#include <iostream>
#include "MultiArrayIter.h"

void test(bool rowMajor) {
  std::vector<size_t> lo(3);
  lo[0] = 1; lo[1] = 2; lo[2] = 3;
  std::vector<size_t> hi(3);
  hi[0] = 4; hi[1] = 7; hi[2] = 8;
  MultiArrayIter it(lo, hi, rowMajor);
  for (size_t i = 0; i < it.getNumberOfTerms(); ++i) {
    std::vector<size_t> inds = it.getIndices();
    std::cout << i << ": ";
    for (size_t j = 0; j < inds.size(); ++j) {
      std::cout << inds[j] << ", ";
    }
    std::cout << '\n';
    it.next();
  }

}

int main() {
  std::cout << "row major test:\n";
  test(true);
  std::cout << "column major test:\n";
  test(false);
}

