
if (NDIMS > 0) {
	upI = dimProdPtr[0] * ((inds[0] + upDirectionPtr[0] + numCellsPtr[0]) % numCellsPtr[0]);
	upI += dimProdPtr[1] * inds[1];
	upI += dimProdPtr[2] * inds[2];
}

