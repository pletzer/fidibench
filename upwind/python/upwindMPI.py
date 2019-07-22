#!/usr/bin/env python

"""
Upwind discretization of advection equation
@author Alexander Pletzer
"""

import pnumpy
import numpy
import sys
from mpi4py import MPI
import sys
import operator
import saveVTK

class Upwind: 

  def __init__(self, velocity, lengths, numCells):

    self.rk = MPI.COMM_WORLD.Get_rank()
    self.sz = MPI.COMM_WORLD.Get_size()

    # decomposition
    self.dc = pnumpy.CubeDecomp(self.sz, numCells)
    if not self.dc.getDecomp():
      print('*** No uniform decomposition could be found for {0} processes'.format(self.sz))
      print('*** Please adjust the number of cells {0}'.format(numCells))
      sys.exit(1)

    # begin/end indices of local sub-domain
    self.localSlices = self.dc.getSlab(self.rk)
    self.iBeg = numpy.array([s.start for s in self.localSlices])
    self.iEnd = numpy.array([s.stop for s in self.localSlices])
    self.nsLocal = numpy.array([s.stop - s.start for s in self.localSlices])
    print('[{0}] local number of cells: {1}'.format(self.rk, self.nsLocal))

    # global number of cells
    self.numCells = numCells

    self.ndims = len(velocity)
    self.deltas = numpy.zeros( (self.ndims,), numpy.float64 )
    self.upDirection = numpy.zeros( (self.ndims,), numpy.int )
    self.v = velocity
    self.lengths = lengths

    # number of local field values
    self.ntot = 1
    for j in range(self.ndims):
      self.upDirection[j] = -1
      if velocity[j] < 0.: self.upDirection[j] = +1
      self.deltas[j] = lengths[j] / numCells[j]
      self.ntot *= self.nsLocal[j]

    self.coeff = self.v * self.upDirection / self.deltas

    # initializing the field
    self.f = pnumpy.gdaZeros( self.nsLocal, numpy.float64, numGhosts=1 )
    self.fOld = pnumpy.gdaZeros( self.nsLocal, numpy.float64, numGhosts=1 )

    # initialize lower corner to one
    if self.rk == 0:
      self.f[0, 0, 0] = 1

    # get the neighboring ranks
    self.neighbSide = [[] for i in range(self.ndims)]
    direction = numpy.array([0] * self.ndims)
    self.neighbRk = numpy.array([0] * self.ndims)
    periodic = [True for i in range(self.ndims)]
    for i in range(self.ndims):
      direction[i] = self.upDirection[i]
      self.neighbRk[i] = self.dc.getNeighborProc(self.rk, direction, periodic=periodic)
      self.neighbSide[i] = tuple(-direction)
      direction[i] = 0

  def advect(self, deltaTime):
    """
    Advance the field by one time step
    """

    self.fOld[:] = self.f
    c = deltaTime * numpy.sum(self.coeff)

    # handle all local computations first
    self.f += c*self.fOld

    self.f[1:, :, :] -= deltaTime*self.coeff[0]*self.fOld[:-1, :, :]
    self.f[:, 1:, :] -= deltaTime*self.coeff[1]*self.fOld[:, :-1, :]
    self.f[:, :, 1:] -= deltaTime*self.coeff[2]*self.fOld[:, :, :-1]

    # fetch neighboring data. This is where communication takes place
    self.f[:1, :, :] -= deltaTime*self.coeff[0]* \
                        self.fOld.getData(self.neighbRk[0], self.neighbSide[0])
    self.f[:, :1, :] -= deltaTime*self.coeff[1]* \
                        self.fOld.getData(self.neighbRk[1], self.neighbSide[1])
    self.f[:, :, :1] -= deltaTime*self.coeff[2]* \
                        self.fOld.getData(self.neighbRk[2], self.neighbSide[2])

  def checksum(self):
    return self.f.reduce(operator.add, 0.0, rootPe=0)

  def printOut(self):
    for i in range(len(self.f)):
      print('{0} {1}'.format(i, self.f[i]))

  def __del__(self):
    self.f.free()
    self.fOld.free()

  def gatherRoot(self):
    """
    Gather the data on process root
    @return array on rank 0, None on other ranks
    """
    res = None
    if self.rk == 0:
      res = numpy.zeros(self.numCells, numpy.float64)

    fRoot = MPI.COMM_WORLD.gather(self.f, root=0)

    if self.rk == 0:
      for rk in range(self.sz):
        slab = self.dc.getSlab(rk)
        res[slab] = fRoot[rk]

    return res

############################################################################################################
def main():
  import sys

  if len(sys.argv) <= 1:
    print("must specify number of cells in each direction.")
    return sys.exit(1)

  ndims = 3
  numCells = [int(sys.argv[1])] * 3

  numTimeSteps = 100
  if len(sys.argv) > 2:
    numTimeSteps = int(sys.argv[2])

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
  if up.rk == 0:
    print("number of cells: {0}".format(numCells))

  tic, tac = 0, 0
  if up.rk == 0: 
    tic = MPI.Wtime()

  # time iterations
  for i in range(numTimeSteps):
    up.advect(dt)

  if up.rk == 0:
    toc = MPI.Wtime()
    print('Wall clock time spent in advection loop: {0} [sec]'.format(toc - tic))

  chksum = up.checksum()
  if up.rk == 0: 
    print("check sum: {0}".format(chksum))

  if doVtk:
    data = up.gatherRoot()
    if up.rk == 0:
      xAxis = numpy.array([0.0 + i*up.deltas[0] for i in range(numCells[0] + 1)])
      yAxis = numpy.array([0.0 + j*up.deltas[1] for j in range(numCells[1] + 1)])
      zAxis = numpy.array([0.0 + k*up.deltas[1] for k in range(numCells[2] + 1)])
      saveVTK.rectilinear('upMPI.vtk', xAxis, yAxis, zAxis, data)

if __name__ == '__main__': main()
