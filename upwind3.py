#/usr/bin/env python

import numpy
import copy

class Upwind: 

  def __init__(self, velocity, lengths, numCells):
    self.numCells = numCells
    self.ndims = len(velocity)
    self.deltas = numpy.zeros( (self.ndims,), numpy.float64 )
    self.upDirection = numpy.zeros( (self.ndims,), numpy.float64 )
    self.v = velocity
    self.lengths = lengths
    self.ntot = 1
    self.offsetMat = numpy.identity(self.ndims, numpy.int)
    self.numCellsExt = numpy.outer(self.numCells, numpy.ones((self.ndims,), numpy.int))
    for j in range(self.ndims):
      self.upDirection[j] = -1
      if velocity[j] < 0.: self.upDirection[j] = +1
      self.deltas[j] = lengths[j] / numCells[j]
      self.offsetMat[j, j] = self.upDirection[j]
      self.ntot *= numCells[j]

    self.dimProd = numpy.ones( (self.ndims,), numpy.int )
    for i in range(self.ndims - 2, -1, -1):
        # last index varies fastest
        self.dimProd[i] =  self.dimProd[i + 1] * self.numCells[i + 1]

    self.coeff = self.v * self.upDirection / self.deltas

    # initializing the field
    self.f = numpy.zeros( (self.ntot,), numpy.float64 )
    # initialize lower corner to one
    self.f[0] = 1

  def advect(self, deltaTime):

    # copy 
    oldF = self.f.copy()

    # array of index sets for each cell
    inds = numpy.zeros( (self.ndims, self.ntot), numpy.int )
    for j in range(self.ndims):
      inds[j, :] = numpy.arange(self.ntot)
      inds[j, :] /= self.dimProd[j]
      inds[j, :] %= self.numCells[j]

    indsUp = inds.copy()

    # update the field in each spatial direction
    for j in range(self.ndims):

      # apply offset
      indsUp[j, :] += self.upDirection[j]
      indsUp[j, :] %= self.numCells[j]

      # compute flat indices corresponding to the offset index sets
      flatIndsUp = numpy.dot(self.dimProd, indsUp)

      # update
      self.f -= (deltaTime * self.coeff[j]) * (oldF[flatIndsUp] - oldF)

      # reset
      indsUp[j, :] = inds[j, :]

  def saveVTK(self, fname):
    f = open(fname, 'w')
    print >> f, "# vtk Data Version 2.0"
    print >> f, "upwind3.py"
    print >> f, "ASCII"
    print >> f, "DATASET RECTILINEAR_GRID"
    print >> f, "DIMENSIONS", 
    # in VTK the first dimension varies fastest so need 
    # to invert the order of the dimensions
    if self.ndims > 2:
      print >> f, '%d' % (self.numCells[2] + 1),
    else:
      print >> f, " 1",
    if self.ndims > 1:
      print >> f, '%d' % (self.numCells[1] + 1),
    else:
      print >> f, " 1",
    print >> f, '%d' % (self.numCells[0] + 1)
    print >> f, "X_COORDINATES",
    if self.ndims > 2:
      print >> f, "%d double" % (self.numCells[2] + 1)
      for i in range(self.numCells[2] + 1):
        print >> f, " %f" % (0.0 + self.deltas[2] * i),   
    else:
      print >> f, "1 double"
      print >> f, "0.0"
    print >> f, "\nY_COORDINATES",
    if self.ndims > 1:
      print >> f, "%d double" % (self.numCells[1] + 1)
      for i in range(self.numCells[1] + 1): 
        print >> f, " %f" % (0.0 + self.deltas[1] * i),  
    else:
      print >> f, "1 double"
      print >> f, "0.0"
    print >> f, "\nZ_COORDINATES",
    print >> f, "%d double" % (self.numCells[0] + 1)
    for i in range(self.numCells[0] + 1):
      print >> f, " %f" % (0.0 + self.deltas[0] * i),
    print >> f, "\nCELL_DATA %d" % self.ntot
    print >> f, "SCALARS f double 1"
    print >> f, "LOOKUP_TABLE default"
    for i in range(self.ntot):
      print >> f, self.f[i],
      if (i + 1) % 10 == 0: 
        print >> f
    print >> f
    f.close()

  def checksum(self):
    return numpy.sum(self.f)

  def printOut(self):
    for i in range(len(self.f)):
      print i, ' ', self.f[i]

############################################################################################################
def main():
  import sys

  if len(sys.argv) <= 1:
    print "must specify number of cells in each direction."
    return sys.exit(1)

  ndims = 3
  numCells = [int(sys.argv[1])] * 3
  print "number of cells: ", numCells

  numTimeSteps = 100
  if len(sys.argv) > 2:
    numTimeSteps = int(sys.argv[2])


  velocity = numpy.ones( (ndims,), numpy.float64 )
  lengths = numpy.ones( (ndims,), numpy.float64 )

  # compute dt 
  courant = 0.1
  dt = float('inf')
  for j in range(ndims):
    dx = lengths[j]/ float(numCells[j])
    dt = min(courant * dx / velocity[j], dt)

  up = Upwind(velocity, lengths, numCells)
  #up.saveVTK("up0.vtk")
  for i in range(numTimeSteps):
    up.advect(dt)
    #if i % 10 == 0:
    #up.saveVTK("up" + str(i) + ".vtk")

  #up.printOut()
  print "check sum: ", up.checksum()
  up.saveVTK("up.vtk")

if __name__ == '__main__': main()
