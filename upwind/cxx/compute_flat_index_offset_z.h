
if (NDIMS > 2) {
	upI = dimProdPtr[0] * inds[0];
	upI += dimProdPtr[1] * inds[1];
	upI += dimProdPtr[2] * ((inds[2] + upDirectionPtr[2] + numCellsPtr[2]) % numCellsPtr[2]);
}

