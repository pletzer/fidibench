
#if (NDIMS > 0)
inds[0] = i / this->dimProd[0] % this->numCells[0];
#endif

#if (NDIMS > 1)
inds[1] = i / this->dimProd[1] % this->numCells[1];
#endif

#if (NDIMS > 2)
inds[2] = i / this->dimProd[2] % this->numCells[2];
#endif