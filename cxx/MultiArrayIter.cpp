/**
 * @brief Multi-dimensional array iterator. Write a single loop and iterate 
 *        over as many dimensions as you like.
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include "MultiArrayIter.h"

MultiArrayIter::MultiArrayIter(const std::vector<size_t>& lo, const std::vector<size_t>& hi, bool rowMajor) {
	this->build(lo, hi, rowMajor);
}

void 
MultiArrayIter::build(const std::vector<size_t>& lo, const std::vector<size_t>& hi, bool rowMajor) {

	this->lo = lo;
	this->hi = hi;
	this->ndims = lo.size();

	this->dims.resize(this->ndims);
	for (size_t i = 0; i < this->ndims; ++i) {
		this->dims[i] = hi[i] - lo[i];
	}

	this->rowMajor = rowMajor;

	// total number of elements
	this->ntot = 1;
	for (size_t i = 0; i < this->ndims; ++i) {
		this->ntot *= dims[i];
	}

	this->dimProd.resize(this->ndims);
	if (rowMajor) {
		this->dimProd[this->ndims - 1] = 1;
		for (int i = (int) this->ndims - 2; i >= 0; --i) {
			this->dimProd[i] = this->dimProd[i + 1] * this->dims[i + 1];
		}
	}
	else {
		this->dimProd[0] = 1;
		for (size_t i = 1; i < this->ndims; ++i) {
			this->dimProd[i] =  this->dimProd[i - 1] * this->dims[i - 1];
		}
	}

	this->begin();
}

size_t 
MultiArrayIter::getNumberOfTerms() const {
	return this->ntot;
}

void 
MultiArrayIter::begin() {
	this->bigIndex = 0;
}

void 
MultiArrayIter::next() {
	this->bigIndex++;
}

void 
MultiArrayIter::move(size_t bi) {
	this->bigIndex = bi;
}

const std::vector<size_t>& 
MultiArrayIter::getLBegIndices() const {
	return this->lo;
}

const std::vector<size_t>& 
MultiArrayIter::getEndIndices() const {
 	return this->hi;
}
