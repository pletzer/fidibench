
if (NDIMS > 0)
	inds[0] = i / dimProdPtr[0] % numCellsPtr[0];

if (NDIMS > 1)
	inds[1] = i / dimProdPtr[1] % numCellsPtr[1];

if (NDIMS > 2)
	inds[2] = i / dimProdPtr[2] % numCellsPtr[2];
