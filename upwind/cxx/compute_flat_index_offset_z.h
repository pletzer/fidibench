
upI = 0;

#if (NDIMS > 0)
upI += dimProdPtr[0] * inds[0];
#endif

#if (NDISM > 1)
upI += dimProdPtr[1] * inds[1];
#endif

#if (NDIMS > 2)
upI += dimProdPtr[2] * ((inds[2] + upDirectionPtr[2] + numCellsPtr[2]) % numCellsPtr[2]);
#endif

