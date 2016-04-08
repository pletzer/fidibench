
type Upwind
    ndims::Int
    velocity::Array
    lengths::Array
    numCells::Array
    deltas::Array
    upDirection::Array
    offsetMat::Array
    numCellsExt::Array
    coeff::Array
    ntot::Int
    dimProd::Array
    f::Array
end

function advect!(up::Upwind,deltaTime)

    # copy
    oldF = deepcopy(up.f)

    # index set for each cell
    inds = zeros(Integer, up.ndims, up.ntot)
    
    for j = 1:up.ndims
        inds[j,:] = collect(1:up.ntot)
        inds[j,:] /= up.dimProd[j]
        inds[j,:] %= up.numCells[j]
    end

    indsUp = deepcopy(inds)
    
    # Update the field in each spatial direction
    for j = 1:up.ndims

        # Apply offset
        indsUp[j,:] += up.upDirection[j]
        indsUp[j,:] %= up.numCells[j]

        #compute flat indices corresponding to the offset index sets
        flatIndsUp = dot(up.dimProd, indsUp)

        #update 
        up.f -= (deltaTime * up.coeff[j])*(oldF[flatIndsUp] - oldF)

        # Reset
        indsUp[j,:] = deepcopy(inds[j,:])
    end
end

function getFlatIndex(up::Upwind,inds)
    ans = Int(dot(up.dimProd,inds))
    return ans
end

function getIndexSet(up::Upwind,flatIndex)
    res = zeros(Int,up.ndims)
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



#Beginning of main program
ndims = 3
#For quick development, I will hard code now
#Later I will take these as parameters.
#ncells = 1
ncells = parse(Int, ARGS[1])
#numTimeSteps = 1
numTimeSteps = parse(Int, ARGS[2])

#Same resolution in each direction.
numCells = [ncells,ncells,ncells]
floatingnumCells = [convert(AbstractFloat,ncells),convert(AbstractFloat,ncells),convert(AbstractFloat,ncells)]
println(numCells)



#Initalizing vector velocity and lengths
velocity = ones(Float64,ndims)
#println(velocity)
lengths   = ones(Float64,ndims)
#println(lengths)


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
deltas = zeros(Float64,ndims)
upDirection = zeros(Integer,ndims)
floatupDirection = zeros(Float64,ndims)
offsetMat = eye(ndims,ndims)
numCellsExt = zeros(Float64,ndims,ndims)
BLAS.ger!(1.0,floatingnumCells,ones(Float64,ndims),numCellsExt)

ntot = 1
for i in 1:ndims
    upDirection[i] = -1
    floatupDirection[i] = -1.00e0
    if velocity[i] < 0.00
        upDirection[i] = +1
        floatupDirection[i] = +1.00e0
    end
    deltas[i] = lengths[i]/numCells[i]
    offsetMat[i,i] = upDirection[1]
    ntot *= numCells[i]
end
#println("upDirection= ",upDirection)
#println("ntot= ",ntot)
#println("delta= ",deltas)
dimProd= ones(Integer,ndims)
for i = ndims-2:-1:0
    dimProd[i+1] = dimProd[i+2]*numCells[i+2]
end
#println("dimProd= ",dimProd)

coeff = zeros(Float64,ndims)
for i in 1:ndims
    coeff[i] = velocity[i] * floatupDirection[i] / deltas[i]
end 

f = zeros(ntot)
f[1] = 1.00e0
up = Upwind(ndims,velocity,lengths,numCells,deltas,upDirection,offsetMat,numCellsExt,coeff,ntot,dimProd,f)
#println(up)


#Do the time evolution of the upwind class
for i = 1:numTimeSteps
   advect!(up,dt) 
end 
println("Time evolution complete")

#Do the checksum
#println("Last f= ",up.f)
println("Check_sum= ",checksum(up))
