################################################################################
# Upwind advection
# 
# Original code from Christopher Blanton (cjb47@psu.edu)
# Modified by Alex Pletzer (alexander@gokliya.net)


using LinearAlgebra
include("saveVTK.jl")

mutable struct Upwind

    # Number of space dimensions
    ndims::Int

    # Constant advection velocity
    velocity::Array

    # Domain length
    lengths::Array

    # Number of cells along each direction
    numCells::Array

    #
    # Internal variables (derived from the above)
    #

    # Grid size along each direction
    deltas::Array

    # Direction of upwind
    upDirection::Array

    # Total number of cells
    ntot::Int

    # Product of the dimensions, used to go from flat index to index set
    dimProd::Array

    # Field
    f::Array

    # Index set
    inds::Array

    # Array of ones
    eins::Array

    function Upwind(ndims::Int, velocity::Array, lengths::Array, numCells::Array)
        """
        Internal constructor
        @param ndims number of space dimensions
        @param velocity velocity array of size ndims
        @param lengths domain lengths array of size ndims
        @param numCells number of cells, array of size ndims
        """
        this = new(ndims, velocity, lengths, numCells)

        this.deltas = zeros(Float64, ndims::Integer)

        this.upDirection = zeros(Integer, ndims::Integer)
        this.ntot = 1
        for i in 1:ndims
            this.upDirection[i] = -1
            if this.velocity[i] < 0.0
                this.upDirection[i] = +1
            end
            this.deltas[i] = lengths[i]/numCells[i]
            this.ntot *= numCells[i]
        end
        this.dimProd = ones(Integer, ndims::Integer)
        for i = 2:ndims
            this.dimProd[i] = this.dimProd[i-1] * this.numCells[i-1]
        end

        # 
        # Set the initial field
        #
        this.f = zeros(this.ntot)
        this.f[1] = 1.0

        # Build the index set
        this.inds = zeros(Integer, (this.ntot, this.ndims))
        for j = 1:ndims
            this.inds[:, j] = mod.(div.(collect(0:this.ntot - 1), this.dimProd[j]), this.numCells[j])
        end

        this.eins = ones(this.ntot)

        return this
    end
end

function advect!(this::Upwind, deltaTime::Float64)
    """
    Advance the solution by one time step
    @param this instance
    @param deltaTime time increment
    """
    oldF = copy(this.f)

    indsUp = copy(this.inds)

    c = deltaTime * this.velocity .* this.upDirection ./ this.deltas

    for j = 1:this.ndims

        indsUp[:, j] = mod.(this.inds[:, j] + this.upDirection[j]*this.eins, this.numCells[j])

        # compute the flat indices corresponding to the offset index sets
        flatIndsUp = *(indsUp, this.dimProd) .+ 1

        # update
        this.f -= c[j] * oldF[flatIndsUp]
        this.f += c[j] * oldF

        # reset
        indsUp[:, j] = this.inds[:, j]

    end

end

function getFlatIndex(this::Upwind, inds::Array)
    """
    Get the flat index from an index set
    @param this instance
    @param inds index set, eg [i, j, k]
    @return integer index of field array
    """
    return dot(this.dimProd, inds - 1) + 1
end
  
function getIndexSet(this::Upwind, flatIndex::Integer)
    """
    Get the index set of a flat index
    @param this instance
    @param flatIndex flat index 
    @return index set, eg [i, j, k]
    """
    res = zeros(Integer, this.ndims)
    for i = 1:this.ndims
        res[i] = mod(div((flatIndex - 1), this.dimProd[i]), this.numCells[i]) + 1
    end
    return res
end

function save2VTK(this::Upwind, fname::AbstractString)
    """
    Save data in VTK file
    @param this instance
    @param fname file name
    """
    xAxis = linspace(0.0, this.lengths[1], this.numCells[1] + 1)
    yAxis = linspace(0.0, this.lengths[2], this.numCells[2] + 1)
    zAxis = linspace(0.0, this.lengths[3], this.numCells[3] + 1)
    saveVTK.rectilinear(fname, xAxis, yAxis, zAxis, this.f)

end

function  main(ncells, numTimeSteps, savevtk)

    ndims = 3

    # Same resolution in each direction.
    numCells = [ncells,ncells,ncells]
    println("number of cells: ", numCells)
    println("number of time steps: ", numTimeSteps)


    # Initalizing vector velocity and lengths
    velocity = ones(Float64, ndims)
    lengths   = ones(Float64, ndims)

    # Create Upwind instance
    up = Upwind(ndims, velocity, lengths, numCells)

    # Compute time step
    courant = 0.1
    dt = Inf
    for i = 1:ndims
        dx = lengths[i]/float(numCells[i])
        dt = min((courant*dx)/velocity[i], dt)
    end

    # Advect
    for i = 1:numTimeSteps
       advect!(up, dt)
    end 

    # Do the checksum
    println("check sum: ", sum(up.f))

    if savevtk
       save2VTK(up, "upwind.vtk")
    end

end

# Beginning of main program
###########################

ncells = 128
numTimeSteps = 2

# Parse command line arguments
if length(ARGS) >= 1
    ncells = parse(Int, ARGS[1])
    if length(ARGS) >= 2
        numTimeSteps = parse(Int, ARGS[2])
    end
end

savevtk = false
if length(ARGS) >= 3 && ARGS[3] == "vtk"
    savevtk = true
end

main(ncells, numTimeSteps, savevtk)
