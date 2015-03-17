/**
 * @brief Multi-dimensional array iterator. Write a single loop and iterate 
 *        over as many dimensions as you like.
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to pletzer@psu.edu.
 */

#include <vector>

#ifndef MULTI_ARRAY_ITER_H
#define MULTI_ARRAY_ITER_H

class MultiArrayIter {

public:

/**
 * No argument constructor
 */
	MultiArrayIter() {}

/**
 * Constructor
 * @param lo low end of the index set
 * @param hi high end (one past last element) of the index set
 * @param rowMajor set this to true for row major indexing, false for column major indexing
 */
	MultiArrayIter(const std::vector<size_t>& lo, const std::vector<size_t>& hi, bool rowMajor) {
		this->build(lo, hi, rowMajor);
	}

/**
 * Build iterator
 * @param lo low end of the index set
 * @param hi high end (one past last element) of the index set
 * @param rowMajor set this to true for row major indexing, false for column major indexing
 */
	void build(const std::vector<size_t>& lo, const std::vector<size_t>& hi, bool rowMajor) {

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

/** 
 * Get the number of terms
 * @return number
 */
	size_t getNumberOfTerms() const {
		return this->ntot;
	}

/**
 * Get te current indices
 * @return indices
 */
	std::vector<size_t> getIndices() const {
		std::vector<size_t> res(this->ndims);
		for (size_t i = 0; i < this->ndims; ++i) {
			res[i] = this->lo[i] + this->bigIndex / this->dimProd[i] % this->dims[i];
		}
		return res;
	}

/**
 * Set the iterator to the beginning
 */
	void begin() {
		this->bigIndex = 0;
	}

/**
 * Increment the iterator
 */
	void next() {
		this->bigIndex++;
	}

/**
 * Move the iterator to the given big index
 * @param bi big (flat) index
 */
	void move(size_t bi) {
		this->bigIndex = bi;
	}

/**
 * Compute the big index gor the given index set
 * @param inds index set
 * @return big index
 */
  template<class T>
	T computeBigIndex(const std::vector<T>& inds) const {
		size_t bi = 0;
		for (size_t i = 0; i < this->ndims; ++i) {
			bi += this->dimProd[i] * (inds[i] - this->lo[i]);
		}
		return bi;
	}

private:

// lo/high end index set
	std::vector<size_t> lo;
	std::vector<size_t> hi;

// dimensions along each axis
	std::vector<size_t> dims;

// used to compute the big index
	std::vector<size_t> dimProd;

// flag set to true for ro major indexing
	bool rowMajor;

// total number of terms
	size_t ntot;

// number of axes
	size_t ndims;

// big index
	size_t bigIndex;

};

#endif // MULTI_ARRAY_ITER_H