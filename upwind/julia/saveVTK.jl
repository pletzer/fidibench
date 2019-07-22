module saveVTK

export rectilinear

function rectilinear(fname::AbstractString,
                     xAxis::LinSpace, yAxis::LinSpace, zAxis::LinSpace,
                     data::Array)
    """
    Save data on rectilinear grid in VTK file
    @param fname file name
    @param xAxis x axis (nodes)
    @param yAxis y axis (nodes)
    @param zAxis z axis (nodes)
    @param data field data
    """
    f = open(fname, "w")
    write(f, "# vtk DataFile Version 2.0\n")
    write(f, "upwind.jl\n")
    write(f, "ASCII\n")
    write(f, "DATASET RECTILINEAR_GRID\n")
    write(f, "DIMENSIONS")

    numCells = [length(xAxis) - 1, length(yAxis) - 1, length(zAxis) - 1]
    ntot = length(data)

    # Number of nodes in each direction
    write(f, @sprintf(" %d", numCells[1] + 1))
    write(f, @sprintf(" %d", numCells[2] + 1))
    write(f, @sprintf(" %d\n", numCells[3] + 1))

    write(f,  "X_COORDINATES ")
    write(f,  @sprintf("%d double\n", numCells[1] + 1))
    for i in 1:numCells[1] + 1
        write(f,  @sprintf("%f ", xAxis[i]))
    end 

    write(f,  "\nY_COORDINATES ")
    write(f, @sprintf("%d double\n", numCells[2] + 1))
    for i in 1:numCells[2] + 1 
        write(f, @sprintf("%f ", yAxis[i]))
    end

    write(f,  "\nZ_COORDINATES ")
    write(f, @sprintf("%d double\n", numCells[3] + 1))
    for i in 1:numCells[3] + 1
        write(f, @sprintf("%f ", zAxis[i]))
    end

    write(f, @sprintf("\nCELL_DATA %d\n", ntot))
    write(f, "SCALARS f double 1\n")
    write(f, "LOOKUP_TABLE default\n")
    for i in eachindex(data)
        write(f, @sprintf("%f ", data[i]))
        if i % 10 == 0
            write(f, "\n")
        end
    end
    write(f, "\n")

    close(f)


end

end
