type Upwind
    ndims::Int
    velocity::Array
    lengths::Array
    numCells::Array
    deltas::Array
    upDirection::Array
    ntot::Int
    dimProd::Array
    f::Array
end

function advect!(up::Upwind,deltaTime)
    oldF = deepcopy(up.f)
    for i = 0:up.ntot-1
        inds = getIndexSet(up,i)
        for j = 1:up.ndims
            oldIndex = inds[j]
            #periodic BCs
            inds[j] += up.upDirection[j] + up.numCells[j] 
            inds[j] %= up.numCells[j]
            upI = getFlatIndex(up,inds)
            up.f[i+1] -= (deltaTime * up.velocity[j] * up.upDirection[j]) * (oldF[upI+1]-oldF[i+1])/up.deltas[j]
            inds[j] = oldIndex
        end
    end 
end

function getFlatIndex(up::Upwind,inds)
    ans = Int(dot(up.dimProd,inds))
    return ans
end

function getIndexSet(up::Upwind,flatIndex)
    res = zeros(Int, up.ndims)
    for i = 0:up.ndims-1
        res[i+1] =  mod(Int(floor(flatIndex/up.dimProd[i+1])), up.numCells[i+1])
    end
    return res
end

function  saveVTK(up::Upwind)

end

function checksum(up::Upwind)
    return sum(up.f)
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

#Same resolution in each direction.
numCells = [ncells,ncells,ncells]
println(numCells)


#Initalizing vector velocity and lengths
velocity = ones(Float64,ndims)
lengths   = ones(Float64,ndims)


#computer dt
courant = 0.1
dt = Inf
dx = 0.00e0
for i = 1:ndims
    dx = lengths[i]/float(numCells[i])
    dt = min((courant*dx)/velocity[i],dt)
end


v = velocity
deltas = zeros(Float64,ndims)
upDirection = zeros(Integer,ndims)
ntot = 1
for i in 1:ndims
    upDirection[i] = -1
    if velocity[i] < 0.00
        upDirection[i] = +1
    end
    deltas[i] = lengths[i]/numCells[i]
    ntot *= numCells[i]
end

dimProd= ones(Integer, ndims)
for i = ndims-2:-1:0
    dimProd[i+1] = dimProd[i+2]*numCells[i+2]
end


f = zeros(ntot)
f[1] = 1.00e0
up = Upwind(ndims,velocity,lengths,numCells,deltas,upDirection,ntot,dimProd,f)


for i = 1:numTimeSteps
   advect!(up,dt) 
end 
println("Time evolution complete")


println("Check_sum= ",checksum(up))
