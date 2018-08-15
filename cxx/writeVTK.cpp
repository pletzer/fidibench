/**
 * @brief Write data on a uniform grid to a VTK file
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include "writeVTK.h"

void writeVTK(const std::string& filename, 
	const std::vector<size_t>& globalDims,
	const std::vector<double>& xmins, 
	const std::vector<double>& xmaxs, 
	const std::vector<double>& field,
	const std::string& name) {

	size_t ndims = globalDims.size();
	if (ndims > 3) {
		std::cerr << "WARNING: writeVTK does not support more than 3 dimensions\n";
		return;
	}

  size_t nx = 1;
  size_t ny = 1;
  size_t nz = 1;
  size_t nx1 = 2;
  size_t ny1 = 2;
  size_t nz1 = 2;
  std::vector<double> xLo(3, 0.0);
  std::vector<double> xHi(3, 1.0);
  std::vector<size_t> zeros(3, 0);
  std::vector<size_t> cellDims(3, 1);
  std::vector<size_t> nodeDims(3, 2);
  if (ndims > 0) {
    nx = globalDims[0];
    nx1 = nx + 1;
    xLo[0] = xmins[0];
    xHi[0] = xmaxs[0];
    nodeDims[0] = nx1;
    cellDims[0] = nx;
  }
  if (ndims > 1) {
    ny = globalDims[1];
    ny1 = ny + 1;
    xLo[1] = xmins[1];
    xHi[1] = xmaxs[1];
    nodeDims[1] = ny1;
    cellDims[1] = ny;
  }
  if (ndims > 2) {
    nz = globalDims[2];
    nz1 = nz + 1;
    xLo[2] = xmins[2];
    xHi[2] = xmaxs[2];
    nodeDims[2] = nz1;
    cellDims[2] = nz;
  }
  size_t numPoints = nx1 * ny1 * nz1;
  size_t numCells = nx * ny * nz;
  
  std::ofstream file;
  file.open(filename.c_str());
  file << "# vtk DataFile Version 2.0\n";
  file << "produced by laplacian\n";
  file << "ASCII\n";
  file << "DATASET STRUCTURED_GRID\n";
  file << "DIMENSIONS " << nx1 << ' ' << ny1 << ' ' << nz1 << '\n';
  file << "POINTS " << numPoints << " float\n";

  // write the vertices
  MultiArrayIter pit(zeros, nodeDims, false); // column major
  for (size_t i = 0; i < pit.getNumberOfTerms(); ++i) {
  	std::vector<size_t> inds = pit.getIndices();
  	double x = xLo[0] + (xHi[0] - xLo[0])*inds[0]/double(cellDims[0]);
  	double y = xLo[1] + (xHi[1] - xLo[1])*inds[1]/double(cellDims[1]);
  	double z = xLo[2] + (xHi[2] - xLo[2])*inds[2]/double(cellDims[2]);
  	file << x << ' ' << y << ' ' << z << '\n';  
  	pit.next();
  }

  // write the cell centered field
  file << "CELL_DATA " << numCells << '\n';
  file << "SCALARS " << name << " float\n";
  file << "LOOKUP_TABLE default\n";
  MultiArrayIter fit(zeros, cellDims, false); // column major
  for (size_t i = 0; i < fit.getNumberOfTerms(); ++i) {
  	file << field[i] << '\n';
  	fit.next();
  }

  // clean up
  file.close();
}
