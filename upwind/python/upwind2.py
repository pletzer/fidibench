#/usr/bin/env python

from __future__ import print_function
import numpy
import copy
import saveVTK

class Upwind: 

  def __init__(self, velocity, lengths, numCells):
    self.numCells = numCells
    self.ndims = len(velocity)
    self.deltas = numpy.zeros( (self.ndims,), numpy.float64 )
    self.upDirection = numpy.zeros( (self.ndims,), numpy.int )
    self.v = velocity
    self.lengths = lengths
    self.ntot = 1
    for j in range(self.ndims):
      self.upDirection[j] = -1
      if velocity[j] < 0.: self.upDirection[j] = +1
      self.deltas[j] = lengths[j] / numCells[j]
      self.ntot *= numCells[j]

    self.dimProd = numpy.ones( (self.ndims,), numpy.int )
    for i in range(self.ndims - 2, -1, -1):
        # last index varies fastest
        self.dimProd[i] =  self.dimProd[i + 1] * self.numCells[i + 1]

    self.f = numpy.zeros( (self.ntot,), numpy.float64 )
    # initialize lower corner to one
    self.f[0] = 1

  def advect(self, deltaTime):

    oldF = copy.deepcopy(self.f)

    for j in range(self.ndims):
      c = deltaTime * self.v[j] * self.upDirection[j] / self.deltas[j]

      for i in range(self.ntot):

        inds = self.getIndexSet(i)
        
        oldIndex = inds[j]
        # periodic BCs
        inds[j] += self.upDirection[j] + self.numCells[j]
        inds[j] %= self.numCells[j]

        upI = self.getFlatIndex(inds)
        self.f[i] -= c * (oldF[upI] - oldF[i])
        inds[j] = oldIndex

  def saveVTK(self, fname):
    xAxis = [0.0]
    yAxis = [0.0]
    zAxis = [0.0]
    if self.ndims > 2:
      xAxis = [0.0 + self.deltas[2] * i for i in range(self.numCells[2] + 1)]
    if self.ndims > 1:
      yAxis = [0.0 + self.deltas[1] * i for i in range(self.numCells[1] + 1)]
    if self.ndims > 0:
      zAxis = [0.0 + self.deltas[0] * i for i in range(self.numCells[0] + 1)]
    saveVTK.rectilinear(fname, xAxis, yAxis, zAxis, self.f.reshape(self.numCells))

  def checksum(self):
    return numpy.sum(self.f)

  def printOut(self):
    for i in range(len(self.f)):
      print(i, ' ', self.f[i])

  def getIndexSet(self, flatIndex):
    res = numpy.zeros( (self.ndims,), numpy.int )
    for i in range(self.ndims):
      res[i] = flatIndex / self.dimProd[i] % self.numCells[i]
    return res

  def getFlatIndex(self, inds):
    return numpy.dot(self.dimProd, inds)

#################################################################################################
def main():
  import sys

  if len(sys.argv) <= 1:
    print("must specify number of cells in each direction.")
    return sys.exit(1)

  ndims = 3
  numCells = [int(sys.argv[1])] * 3
  print("number of cells: ", numCells)

  numTimeSteps = 100
  if len(sys.argv) > 2:
    numTimeSteps = int(sys.argv[2])
  print("number of time steps: ", numTimeSteps)

  doVtk = False
  if len(sys.argv) > 3 and sys.argv[3] == 'vtk':
    doVtk = True

  velocity = numpy.ones( (ndims,), numpy.float64 )
  lengths = numpy.ones( (ndims,), numpy.float64 )

  # compute dt 
  courant = 0.1
  dt = float('inf')
  for j in range(ndims):
    dx = lengths[j]/ float(numCells[j])
    dt = min(courant * dx / velocity[j], dt)

  up = Upwind(velocity, lengths, numCells)
  for i in range(numTimeSteps):
    up.advect(dt)

  print("check sum: ", up.checksum())
  if doVtk:
    up.saveVTK("up.vtk")

if __name__ == '__main__': main()
