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

  f = open(fname, 'w')
  print("# vtk DataFile Version 2.0", file=f)
  print("upwind", file=f)
  print("ASCII", file=f)
  print("DATASET RECTILINEAR_GRID", file=f)
  print("DIMENSIONS", end=' ', file=f)
  # in VTK the first dimension varies fastest so need 
  # to invert the order of the dimensions
  if ndims > 2:
    print(' %d' % (numCells[2] + 1), end=' ', file=f)
  else:
    print(" 1", end=' ', file=f)
  if ndims > 1:
    print(' %d' % (numCells[1] + 1), end=' ', file=f)
  else:
    print(" 1", end=' ', file=f)
  print(' %d' % (numCells[0] + 1), end=' ', file=f)
  print("\nX_COORDINATES ", end=' ', file=f)
  if ndims > 2:
    print(numCells[2] + 1,  " double", file=f)
    for i in range(numCells[2] + 1):
      print(' %f' % xAxis[i], end=' ', file=f) 
  else:
    print("1 double", file=f)
    print("0.0", end=' ', file=f)
  print("\nY_COORDINATES ", end=' ', file=f)
  if ndims > 1:
    print(numCells[1] + 1,  " double", file=f)
    for i in range(numCells[1] + 1): 
      print(' %f' % yAxis[i], end=' ', file=f) 
  else:
    print("1 double", file=f)
    print("0.0", end=' ', file=f)
  print("\nZ_COORDINATES ", end=' ', file=f)
  print(numCells[0] + 1, " double", file=f)
  for i in range(numCells[0] + 1):
    print(' %f' % zAxis[i], end=' ', file=f)

  flatData = data.flat
  ntot = len(flatData)

  print("\nCELL_DATA %d" % ntot, file=f)
  print("SCALARS f double 1", file=f)
  print("LOOKUP_TABLE default", file=f)
    
  for i in range(len(flatData)):
    print(flatData[i], end=' ', file=f)
    if (i + 1) % 10 == 0:
        print(file=f)
  f.close()
