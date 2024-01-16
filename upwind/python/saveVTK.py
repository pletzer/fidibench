#!/usr/bin/env python

def rectilinear(fname, xAxis, yAxis, zAxis, data):
  """
  Save data in reclilinear VTK file format
  @param fname file name
  @param xAxis x axis
  @param yAxis y axis
  @param zAxis z axis
  @param data 
  """

  numCells = data.shape
  ndims = len(numCells)

  with open(fname, 'w') as f:
    f.write("# vtk DataFile Version 2.0\n")
    f.write("upwind\n")
    f.write("ASCII\n")
    f.write("DATASET RECTILINEAR_GRID\n")
    f.write("DIMENSIONS")
    # in VTK the first dimension varies fastest so need 
    # to invert the order of the dimensions
    if ndims > 2:
      f.write(f' {numCells[2] + 1}')
    else:
      f.write(" 1")
    if ndims > 1:
      f.write(f' {numCells[1] + 1}')
    else:
      f.write(" 1")
    f.write(f' {numCells[0] + 1}')
    f.write("\nX_COORDINATES ")
    if ndims > 2:
      f.write(f"{numCells[2] + 1} double")
      for i in range(numCells[2] + 1):
        f.write(f' {xAxis[i]}')
    else:
      f.write("1 double 0.0")
    f.write( "\nY_COORDINATES ")
    if ndims > 1:
      f.write(f"{numCells[1] + 1} double")
      for i in range(numCells[1] + 1): 
        f.write(f' {yAxis[i]}') 
    else:
      f.write("1 double 0.0")
    f.write("\nZ_COORDINATES ")
    f.write(f"{numCells[0] + 1} double")
    for i in range(numCells[0] + 1):
      f.write(f' {zAxis[i]}')

    flatData = data.flat
    ntot = len(flatData)

    f.write(f"\nCELL_DATA {ntot}\n")
    f.write("SCALARS f double 1\n")
    f.write("LOOKUP_TABLE default\n")
      
    for i in range(len(flatData)):
      f.write(f'{flatData[i]}')
      if (i + 1) % 10 == 0:
          f.write('\n')
