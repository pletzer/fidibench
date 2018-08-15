/**
 * @brief find the most optimal domain decomposition 
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */
#include "CubeDecomp.h"

bool 
CubeDecomp::build(int nprocs, const std::vector<size_t>& dims) {

    this->nprocs = nprocs;
    this->ndims = dims.size();
    this->dims = dims;

    this->decomp = this->computeOptimalDecomp();

    if (this->decomp.size() == 0) {
      // could not find a valid decomposition
      return false;
    }

    this->localDims.resize(decomp.size());
    for (size_t i = 0; i < this->localDims.size(); ++i) {
      this->localDims[i] = dims[i] / this->decomp[i];
    }

    return true;
}

std::vector<size_t> 
CubeDecomp::computeOptimalDecomp() {

    std::vector< std::vector<size_t> > primeNumbers;
    for (size_t i = 0; i < this->ndims; ++i) {
      std::vector<size_t> primeFacts = this->getPrimeFactors(dims[i]);
      primeNumbers.push_back(primeFacts);
    }

    std::vector<size_t> zeros(this->ndims);
    std::vector<size_t> ns(this->ndims);
    for (size_t i = 0; i < this->ndims; ++i) {
      ns[i] = primeNumbers[i].size();
    }

    // gather astd::vector<size_t> decompll the valid decompositions
    std::vector< std::vector<size_t> > validDecomps;
    MultiArrayIter it(zeros, ns, false);
    std::vector<size_t> decomp(this->ndims);
    for (size_t i = 0; i < it.getNumberOfTerms(); ++i) {
      std::vector<size_t> decompInds = it.getIndices();
      std::vector<size_t> decomp(this->ndims);
      for (size_t j = 0; j < this->ndims; ++j) {
        decomp[j] = primeNumbers[j][ decompInds[j] ];
      }
      size_t n = std::accumulate(decomp.begin(), decomp.end(), 1, std::multiplies<size_t>());
      if (n == this->nprocs) {
        validDecomps.push_back(decomp);
      }
      it.next();
    }

    // no decomp was found
    if (validDecomps.size() == 0) {
      // return empty vector
      return std::vector<size_t>();
    }

    // find the optimal decomp
    std::vector<size_t> bestDecomp = validDecomps[0];
    double minCost = std::numeric_limits<double>::max();
    for (size_t i = 1; i < validDecomps.size(); ++i) {
      double cost = this->computeCost(validDecomps[i]);
      if (cost < minCost) {
        minCost = cost;
        bestDecomp = validDecomps[i];
      }
    }

    // create decomp iterator
    this->mit.build(zeros, bestDecomp, true);

    return bestDecomp;
}

std::vector<size_t> 
CubeDecomp::getBegIndices(int rk) {
    size_t bi = (size_t) rk;
    this->mit.move(bi);
    std::vector<size_t> inds = this->mit.getIndices();
    std::vector<size_t> res(this->ndims);
    for (size_t i = 0; i < this->ndims; ++i) {
      res[i] = inds[i] * localDims[i];
    }
    return res;
}

std::vector<size_t> 
CubeDecomp::getEndIndices(int rk) {
    this->mit.move( (size_t) rk );
    std::vector<size_t> inds = this->mit.getIndices();
    std::vector<size_t> res(this->ndims);
    for (size_t i = 0; i < this->ndims; ++i) {
      res[i] = (inds[i] + 1) * localDims[i];
    }
    return res;
}

int 
CubeDecomp::getNeighborRank(int rk, const std::vector<int>& dir) {

    this->mit.move( (size_t) rk );
    std::vector<size_t> inds = this->mit.getIndices();
    std::vector<int> indOffset(inds.size());
    for (size_t i = 0; i < this->ndims; ++i) {

      indOffset[i] = inds[i] + dir[i];

      // apply periodic BCs outside of domain
      if (indOffset[i] < 0) {
        indOffset[i] += this->decomp[i];
      }
      else if (indOffset[i] >= this->decomp[i]) {
        indOffset[i] -= this->decomp[i];
      }

    }
    return (int) this->mit.computeBigIndex(indOffset);
}

std::vector<size_t> 
CubeDecomp::getDecomp() const {
    return this->decomp;
}

std::vector<size_t> 
CubeDecomp::getPrimeFactors(size_t n) const {
    std::vector<size_t> res(1, 1);
    size_t n2 = n / 2;
    for (size_t k = 2; k <= n2; ++k) {
      if ((n / k)*k == n) {
        res.push_back(k);
      }
    }
    res.push_back(n);
    return res;
}

double 
CubeDecomp::computeCost(const std::vector<size_t>& decomp) const {

    std::vector<size_t> localDims(this->ndims);
    for (size_t i = 0; i < this->ndims; ++i) {
      localDims[i] = this->dims[i] / decomp[i];
    }

    size_t volume = 1;
    size_t surface = 0;
    for (size_t i = 0; i < this->ndims; ++i) {
      volume *= localDims[i];
      size_t area = 1;
      for (size_t j = 0; j < this->ndims; ++j) {
        if (j != i) area *= localDims[j];
      }
      surface += area;
    }
    return (double) surface / (double) volume;
}


