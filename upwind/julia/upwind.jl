################################################################################
# Upwind advection
# 
# Original code from Christopher Blanton (cjb47@psu.edu)
# Modified by Alex Pletzer (alexander@gokliya.net)

include("saveVTK.jl")
import saveVTK

type Upwind
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

    #
    # Internal constructor
    #
    function Upwind(ndims::Int, velocity::Array, lengths::Array, numCells::Array)

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
        this.dimProd = ones(Integer,ndims::Integer)
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

function advect!(up::Upwind, deltaTime::Float64)
    oldF = deepcopy(up.f)
    for j = 1:up.ndims
        c = deltaTime * up.velocity[j] * up.upDirection[j] / up.deltas[j]
        for i = 1:up.ntot
            inds = getIndexSet(up, i)
            oldIndex = inds[j]
            #periodic BCs
            inds[j] += up.upDirection[j]
            inds[j] = mod1(inds[j], up.numCells[j])
            upI = getFlatIndex(up, inds)
            up.f[i] -= c * (oldF[upI] - oldF[i])
            inds[j] = oldIndex
        end
    end 
end

function getFlatIndex(up::Upwind, inds::Array)
    return dot(up.dimProd, inds - 1) + 1
end
  
function getIndexSet(up::Upwind, flatIndex::Integer)
    res = zeros(Integer,up.ndims)
    for i = 1:up.ndims
        res[i] = mod(div((flatIndex - 1), up.dimProd[i]), up.numCells[i]) + 1
    end
    return res
end

function checksum(up::Upwind)
    return sum(up.f)
end

function save2VTK(up::Upwind, fname::AbstractString)

    xAxis = linspace(0.0, up.lengths[1], up.numCells[1] + 1)
    yAxis = linspace(0.0, up.lengths[2], up.numCells[2] + 1)
    zAxis = linspace(0.0, up.lengths[3], up.numCells[3] + 1)
    saveVTK.rectilinear(fname, xAxis, yAxis, zAxis, up.f)

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
dx = 0.0
for i = 1:ndims
    dx = lengths[i]/float(numCells[i])
    dt = min((courant*dx)/velocity[i], dt)
end

# Create Upwind instance
up = Upwind(ndims,velocity,lengths,numCells)

# Advect
for i = 1:numTimeSteps
   advect!(up, dt)
end 

# Do the checksum
println("check sum: ",checksum(up))

if length(ARGS) > 2 && ARGS[3] == "vtk"
   save2VTK(up, "upwind.vtk")
end
