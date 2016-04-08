
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

function advect!(up::Upwind,deltaTime::Float64)
    oldF = deepcopy(up.f)
    for i = 1:up.ntot
        inds = getIndexSet(up,i)
        for j = 1:up.ndims
            oldIndex = inds[j]
            #periodic BCs
            inds[j] += up.upDirection[j] + up.numCells[j] 
            inds[j] = mod1(inds[j],up.numCells[j])
            upI = getFlatIndex(up,inds)
            up.f[i] -= (deltaTime * up.velocity[j] * up.upDirection[j]) * (oldF[upI]-oldF[i])/up.deltas[j]
            inds[j] = oldIndex
        end

    end 
end


function getFlatIndex(up::Upwind,inds::Array)
    ans = mod1(dot(up.dimProd,inds),up.ntot)
    return ans
end

  
function getIndexSet(up::Upwind,flatIndex::Integer)
    res = zeros(Integer,up.ndims)
    for i = 1:up.ndims
        res[i] = mod1(div((flatIndex-1),up.dimProd[i]), up.numCells[i])
    end
    return res
end
 



function checksum(up::Upwind)
    return sum(up.f)
end

#Beginning of main program
ndims = 3
ncells = parse(Int, ARGS[1])
numTimeSteps = parse(Int, ARGS[2])

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
#println(dx)
#println(dt)

#This is where upwind class is constructed in C++ and python version of this code
v = velocity
deltas = zeros(Float64,ndims::Integer)
upDirection = zeros(Integer,ndims::Integer)
ntot = 1
for i in 1:ndims
    upDirection[i] = -1
    if velocity[i] < 0.00
        upDirection[i] = +1
    end
    deltas[i] = lengths[i]/numCells[i]
    ntot *= numCells[i]
end
dimProd= ones(Integer,ndims::Integer)
for i = 2:ndims
   dimProd[i] = dimProd[i-1]*numCells[i-1]
end

f = zeros(ntot)
f[1] = 1.00e0
up = Upwind(ndims,velocity,lengths,numCells,deltas,upDirection,ntot,dimProd,f)


for i = 1:numTimeSteps
   advect!(up,dt)
end 
println("Time evolution complete")

#Do the checksum
println("Check_sum= ",checksum(up))
