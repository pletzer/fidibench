#/usr/bin/env python

import numpy
cimport numpy as np
cimport cython

ITYPE = numpy.int
DTYPE = numpy.double
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.double_t DTYPE_t
ctypedef np.int_t ITYPE_t

cdef class Upwind: 

  cdef np.ndarray numCells
  cdef np.ndarray deltas
  cdef np.ndarray upDirection
  cdef np.ndarray v
  cdef np.ndarray lengths
  cdef np.ndarray offsetMat
  cdef np.ndarray numCellsExt
  cdef np.ndarray dimProd
  cdef np.ndarray coeff
  cdef np.ndarray f
  cdef np.ndarray oldF
  cdef int ndims
  cdef int ntot
  cdef np.ndarray inds
  cdef np.ndarray indsUp
  cdef np.ndarray flatIndsUp

  def __init__(self, np.ndarray velocity, 
                     np.ndarray lengths, 
                     np.ndarray numCells):
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
    self.oldF = numpy.zeros( (self.ntot,), numpy.float64 )
    # initialize lower corner to one
    self.f[0] = 1

    # array of index sets for each cell
    self.inds = numpy.zeros( (self.ndims, self.ntot), numpy.int )
    for j in range(self.ndims):
      self.inds[j, :] = numpy.arange(self.ntot)
      self.inds[j, :] /= self.dimProd[j]
      self.inds[j, :] %= self.numCells[j]

    # array of index sets with offsets
    self.indsUp = self.inds.copy()

    # flat indices
    self.flatIndsUp = numpy.zeros( (self.ntot,), numpy.int )

  @cython.boundscheck(False)
  cpdef advect(self, double deltaTime):

    # copy 
    self.oldF[:] = self.f

    # update the field in each spatial direction
    for j in range(self.ndims):

      # apply offset
      self.indsUp[j, :] += self.upDirection[j]
      self.indsUp[j, :] %= self.numCells[j]

      # compute flat indices corresponding to the offset index sets
      self.flatIndsUp = numpy.dot(self.dimProd, self.indsUp)

      # update
      self.f -= (deltaTime * self.coeff[j]) * (self.oldF[self.flatIndsUp] - self.oldF)

      # reset
      self.indsUp[j, :] = self.inds[j, :]


  cdef saveVTK(self, fname):
    f = open(fname, 'w')
    print >> f, "# vtk Data Version 2.0"
    print >> f, "upwind.cxx"
    print >> f, "ASCII"
    print >> f, "DATASET RECTILINEAR_GRID"
    print >> f, "DIMENSIONS"
    # in VTK the first dimension varies fastest so need 
    # to invert the order of the dimensions
    if self.ndims > 2:
      print >> f, ' %d' % (self.numCells[2] + 1),
    else:
      print >> f, " 1",
    if self.ndims > 1:
      print >> f, ' %d' % (self.numCells[1] + 1),
    else:
      print >> f, " 1",
    print >> f, ' %d\n' % (self.numCells[0] + 1)
    print >> f, "X_COORDINATES "
    if self.ndims > 2:
      print >> f, self.numCells[2] + 1,  " double"
      for i in range(self.numCells[2] + 1):
        print >> f, ' %f' % (0.0 + self.deltas[2] * i)    
    else:
      print >> f, "1 double"
      print >> f, "0.0"
    print >> f, "Y_COORDINATES "
    if self.ndims > 1:
      print >> f, self.numCells[1] + 1,  " double"
      for i in range(self.numCells[1] + 1): 
        print >> f, ' %f' % (0.0 + self.deltas[1] * i)  
    else:
      print >> f, "1 double"
      print >> f, "0.0"
    print >> f, "Z_COORDINATES "
    print >> f, self.numCells[0] + 1, " double"
    for i in range(self.numCells[0] + 1):
      print >> f, ' %f' % (0.0 + self.deltas[0] * i)
    print >> f, "CELL_DATA %d" % self.ntot
    print >> f, "SCALARS f double 1"
    print >> f, "LOOKUP_TABLE default"
    for i in range(self.ntot):
      print >> f, self.f[i]
    f.close()

  cpdef checksum(self):
    return self.f.sum()

  cpdef printOut(self):
    for i in range(self.ntot):
      print i, ' ', self.f[i]
