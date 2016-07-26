
if (NDIMS > 1) {
	upI = dimProdPtr[0] * inds[0];
	upI += dimProdPtr[1] * ((inds[1] + upDirectionPtr[1] + numCellsPtr[1]) % numCellsPtr[1]);
	upI += dimProdPtr[2] * inds[2];
}

