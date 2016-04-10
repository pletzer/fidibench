################################################################################
# Upwind advection
# 
# Original code from Christopher Blanton (cjb47@psu.edu)
# Modified by Alex Pletzer (alexander@gokliya.net)
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
            inds[j] += up.upDirection[j] + up.numCells[j] 
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

function saveVTK(up::Upwind, fname::AbstractString)

    f = open(fname, "w")
    write(f, "# vtk DataFile Version 2.0\n")
    write(f, "upwind2.jl\n")
    write(f, "ASCII\n")
    write(f,  "DATASET RECTILINEAR_GRID\n")
    write(f,  "DIMENSIONS")

    # Number of nodes in each direction
    if up.ndims >= 1
        write(f, @sprintf(" %d", up.numCells[1] + 1))
        if up.ndims >= 2
            write(f, @sprintf(" %d", up.numCells[2] + 1))
            if up.ndims >= 3
                write(f, @sprintf(" %d\n", up.numCells[3] + 1))
            else
                write(f, " 1\n")
            end
        else
            write(f, " 1")
        end
    else
        write(f, " 1")
    end

    write(f,  "X_COORDINATES ")
    if up.ndims >= 1
      write(f,  @sprintf("%d double\n", up.numCells[1] + 1))
      for i in 1:up.numCells[1] + 1
        write(f,  @sprintf("%f ", 0.0 + up.deltas[1] * (i-1)))
      end 
    else
      write(f,  "1 double\n")
      write(f,  "0.0\n")
    end

    write(f,  "\nY_COORDINATES ")
    if up.ndims >= 2
      write(f, @sprintf("%d double\n", up.numCells[2] + 1))
      for i in 1:up.numCells[2] + 1 
        write(f, @sprintf("%f ", 0.0 + up.deltas[2] * (i-1)))
      end
    else
        write(f,  "1 double\n")
        write(f,  "0.0\n")
    end

    write(f,  "\nZ_COORDINATES ")
    if up.ndims >= 3
        write(f, @sprintf("%d double\n", up.numCells[3] + 1))
        for i in 1:up.numCells[3] + 1
            write(f, @sprintf("%f ", 0.0 + up.deltas[3] * (i-1)))
        end
    else
        write(f,  "1 double\n")
        write(f,  "0.0\n")
    end

    write(f, @sprintf("\nCELL_DATA %d\n", up.ntot))
    write(f, "SCALARS f double 1\n")
    write(f, "LOOKUP_TABLE default\n")
    for i in eachindex(up.f)
        write(f, @sprintf("%f ", up.f[i]))
        if i % 10 == 0
            write(f, "\n")
        end
    end
    write(f, "\n")

    close(f)
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
   saveVTK(up, "upwind.vtk")
