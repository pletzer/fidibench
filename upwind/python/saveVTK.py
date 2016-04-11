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
  print >> f, "# vtk DataFile Version 2.0"
  print >> f, "upwind"
  print >> f, "ASCII"
  print >> f, "DATASET RECTILINEAR_GRID"
  print >> f, "DIMENSIONS",
  # in VTK the first dimension varies fastest so need 
  # to invert the order of the dimensions
  if ndims > 2:
    print >> f, ' %d' % (numCells[2] + 1),
  else:
    print >> f, " 1",
  if ndims > 1:
    print >> f, ' %d' % (numCells[1] + 1),
  else:
    print >> f, " 1",
  print >> f, ' %d' % (numCells[0] + 1),
  print >> f, "\nX_COORDINATES ",
  if ndims > 2:
    print >> f, numCells[2] + 1,  " double"
    for i in range(numCells[2] + 1):
      print >> f, ' %f' % xAxis[i], 
  else:
    print >> f, "1 double"
    print >> f, "0.0",
  print >> f, "\nY_COORDINATES ",
  if ndims > 1:
    print >> f, numCells[1] + 1,  " double"
    for i in range(numCells[1] + 1): 
      print >> f, ' %f' % yAxis[i], 
  else:
    print >> f, "1 double"
    print >> f, "0.0",
  print >> f, "\nZ_COORDINATES ",
  print >> f, numCells[0] + 1, " double"
  for i in range(numCells[0] + 1):
    print >> f, ' %f' % zAxis[i],

  flatData = data.flat
  ntot = len(flatData)

  print >> f, "\nCELL_DATA %d" % ntot
  print >> f, "SCALARS f double 1"
  print >> f, "LOOKUP_TABLE default"
    
  for i in range(len(flatData)):
    print >> f, flatData[i],
    if (i + 1) % 10 == 0:
        print >> f
  f.close()
