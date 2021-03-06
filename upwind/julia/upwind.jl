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

        return this
    end
end

function advect!(this::Upwind, deltaTime::Float64)
    """
    Advance the solution by one time step
    @param this instance
    @param deltaTime time increment
    """
    oldF = deepcopy(this.f)
    for j = 1:this.ndims
        c = deltaTime * this.velocity[j] * this.upDirection[j] / this.deltas[j]
        for i = 1:this.ntot
            inds = getIndexSet(up, i)
            oldIndex = inds[j]
            #periodic BCs
            inds[j] += this.upDirection[j]
            inds[j] = mod1(inds[j], this.numCells[j])
            upI = getFlatIndex(up, inds)
            this.f[i] -= c * (oldF[upI] - oldF[i])
            inds[j] = oldIndex
        end
    end 
end

function getFlatIndex(this::Upwind, inds::Array)
    """
    Get the flat index from an index set
    @param this instance
    @param inds index set, eg [i, j, k]
    @return integer index of field array
    """
    return dot(this.dimProd, inds .- 1) + 1
end
  
function getIndexSet(this::Upwind, flatIndex::Integer)
    """
    Get the index set of a falt index
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

# Beginning of main program

ndims = 3

# Parse command line arguments
if length(ARGS) < 1
    println("Must specify number of cells in each direction.")
    println("Usage: ", basename(Base.source_path()), " numCells [numTimeSteps]")
    exit(1)
end
ncells = parse(Int, ARGS[1])

numTimeSteps = 100
if length(ARGS) > 1
    numTimeSteps = parse(Int, ARGS[2])
end

# Same resolution in each direction.
numCells = [ncells,ncells,ncells]
println("number of cells: ", numCells)
println("number of time steps: ", numTimeSteps)


# Initalizing vector velocity and lengths
velocity = ones(Float64, ndims)
lengths   = ones(Float64, ndims)


# Compute time step
courant = 0.1
dt = Inf
for i = 1:ndims
    dx = lengths[i]/float(numCells[i])
    global dt = min((courant*dx)/velocity[i], dt)
end

# Create Upwind instance
up = Upwind(ndims,velocity,lengths,numCells)

# Advect
for i = 1:numTimeSteps
   advect!(up, dt)
end 

# Do the checksum
println("check sum: ", sum(up.f))

if length(ARGS) > 2 && ARGS[3] == "vtk"
   save2VTK(up, "upwind.vtk")
end
