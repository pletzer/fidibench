  void saveVTK(const std::string& filename) {
    std::fstream file;
    file.open(filename.c_str(), std::ios_base::out);
    file << "# vtk DataFile Version 2.0\n";
    file << "upwind.cxx\n";
    file << "ASCII\n";
    file << "DATASET RECTILINEAR_GRID\n";
    file << "DIMENSIONS";
    // in VTK the first dimension varies fastest so need 
    // to invert the order of the dimensions
    if (NDIMS > 2) {
      file << ' ' << this->numCells[2] + 1;
    }
    else {
      file << " 1";
    }
    if (NDIMS > 1) {
      file << ' ' << this->numCells[1] + 1;
    }
    else {
      file << " 1";
    }
    file << ' ' << this->numCells[0] + 1;
    file << "\nX_COORDINATES ";
    if (NDIMS > 2) {
      file << this->numCells[2] + 1 << " double\n";
      for (size_t i = 0; i < this->numCells[2] + 1; ++i) {
        file << ' ' << 0.0 + this->deltas[2] * i;
      }      
    }
    else {
      file << "1 double\n";
      file << "0.0\n";
    }
    file << "\nY_COORDINATES ";
    if (NDIMS > 1) {
      file << this->numCells[1] + 1 << " double\n";
      for (size_t i = 0; i < this->numCells[1] + 1; ++i) {
        file << ' ' << 0.0 + this->deltas[1] * i;
      }      
    }
    else {
      file << "1 double\n";
      file << "0.0\n";
    }
    file << "\nZ_COORDINATES ";
    file << this->numCells[0] + 1 << " double\n";
    for (size_t i = 0; i < this->numCells[0] + 1; ++i) {
      file << ' ' << 0.0 + this->deltas[0] * i;
    }
    file << "\nCELL_DATA " << this->ntot << '\n';
    file << "SCALARS f double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < this->ntot; ++i) {
      file << this->f[i] << " ";
      if ((i + 1) % 10 == 0) file << '\n';
    }
    file << '\n';
    file.close();
  }
