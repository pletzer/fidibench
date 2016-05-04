module upwind_mod

    implicit none
    
    integer, parameter :: r8 = selected_real_kind(12, 100)
    
    type upwind_type
        ! number of space dimensions
        integer :: ndims
        
        ! total number of cells in the domain
        integer :: ntot
        
        ! number of cells in the ndims directions
        integer, allocatable :: numCells(:)
        
        ! upwind direction
        integer, allocatable :: upDirection(:)
        
        ! cell sizes in the ndims directions
        real(r8), allocatable :: deltas(:)
        
        ! velocity
        real(r8), allocatable :: v(:)
        
        ! domain lengths
        real(r8), allocatable :: lengths(:)
        
        ! field as a flat array
        real(r8), allocatable :: f(:)
        
        ! product of the dimensions, used to switch back and forth 
        ! between the flat index and the multi-index representations
        integer, allocatable :: dimProd(:)
    end type

contains

    ! Constructor
    ! @param velocity velocity field (constant)
    ! @param lengths domain lengths
    ! @param numCells number of cells in the x, y, ... directions
    subroutine upwind_new(obj, velocity, lengths, numCells)
    
        type(upwind_type) :: obj
        real(r8), intent(in) :: velocity(:)
        real(r8), intent(in) :: lengths(:)
        integer, intent(in) :: numCells(:)
        
        integer :: j
        
        ! f95 will allocate and copy
        obj % ndims = size(velocity)
        obj % v = velocity
        obj % lengths = lengths
        obj % numCells = numCells
        
        allocate(obj % upDirection(obj % ndims))
        allocate(obj % deltas(obj % ndims))
        
        ! compute the otal number of cells and other stuff
        obj % ntot = 1       
        do j = 1, obj % ndims
            obj % upDirection(j) = -1
            if (velocity(j) < 0.) then
                obj % upDirection(j) = 1
            endif
            obj % deltas(j) = lengths(j) / numCells(j)
            obj % ntot = obj % ntot * numCells(j)
        enddo
        
        allocate(obj % dimProd(obj % ndims))
        
        obj % dimProd(obj % ndims) = 1
        do j = obj % ndims - 1, 1, -1
            obj % dimProd(j) =  obj % dimProd(j + 1) * obj % numCells(j + 1)
        enddo
        
        allocate(obj % f(obj % ntot))
        
        ! initialize the field, zero everywhere except for the 
        ! low corner cell where the field is one
        obj % f = 0
        obj % f(1) = 1

    end subroutine 
    
    ! Destructor
    subroutine upwind_del(obj)
    
        type(upwind_type) :: obj
        
        deallocate(obj % v)
        deallocate(obj % lengths)
        deallocate(obj % numCells)    
        deallocate(obj % upDirection)
        deallocate(obj % deltas)
        deallocate(obj % dimProd)
        deallocate(obj % f)
        
    end subroutine
    
    ! Advance by one time step
    ! @param deltaTime time step
    subroutine upwind_advect(obj, deltaTime)
    
        type(upwind_type) :: obj
        real(r8), intent(in) :: deltaTime

        real(r8), allocatable :: oldF(:)
        integer :: i, j, oldIndex, upI
        integer :: inds(obj % ndims)
        
        ! allocate and copy the field (f95)
        oldF = obj % f

        ! iterate over the cells
        do i = 1, obj % ntot

            ! compute the index set of this cell
            call upwind_getIndexSet(obj, i, inds)

            ! iterate over the directions
            do j = 1, obj % ndims
                
                ! cache the cell index
                oldIndex = inds(j)
                
                ! increment the cell index
                inds(j) = inds(j) + obj % upDirection(j)
                
                ! apply periodic BCs
                inds(j) = mod(inds(j) + obj % numCells(j) - 1, obj % numCells(j)) + 1
                  
                ! compute the new flat index 
                call upwind_getFlatIndex(obj, inds, upI)
                    
                ! update the field
                obj % f(upI) = obj % f(upI) - &
              &               deltaTime*obj % v(j)*obj % upDirection(j)*(oldF(upI) - oldF(i))/obj % deltas(j)
                    
                ! reset the index
                inds(j) = oldIndex
           enddo
        enddo

    end subroutine

    subroutine upwind_saveVTK(obj, filename)
        type(upwind_type) :: obj
        character(len=*), intent(in) :: filename
        
        integer iunit, i
        
        ! f2008
        open(newunit = iunit, file = filename, status = 'new')
        write(iunit, *) '# vtk DataFile Version 2.0'
        write(iunit, *) 'upwind.f90'
        write(iunit, *) 'ASCII'
        write(iunit, *) 'DATASET RECTILINEAR_GRID'

        write(iunit, *) 'DIMENSIONS'
        ! in VTK the first dimension varies fastest so need 
        ! to invert the order of the dimensions
        if (obj % ndims > 2) then
            write(iunit, *) ' ', obj % numCells(3) + 1
        else
            write(iunit, *) ' 1'
        endif 
        
        if (obj % ndims > 1) then
            write(iunit, *) ' ', obj % numCells(2) + 1
        else
            write(iunit, *) ' 1'
        endif
        
        write(iunit, *) ' ', obj % numCells(1) + 1

        write(iunit, *) 'X_COORDINATES '
        if (obj % ndims > 2) then
            write(iunit, *) obj % numCells(3) + 1, ' double'
            do i = 1, obj % numCells(3)
                write(iunit, *) ' ', 0.0 + obj % deltas(3) * i
            enddo     
  
        else 
            write(iunit, *) '1 double'
            write(iunit, *) '0.0'
        endif
    
        write(iunit, *) 'Y_COORDINATES '
        if (obj % ndims > 1) then
            write(iunit, *) obj % numCells(2) + 1, ' double'
            do i = 1, obj % numCells(2)
                write(iunit, *) ' ', 0.0 + obj % deltas(2) * i
            enddo      
        else
            write(iunit, *) '1 double'
            write(iunit, *) '0.0'
        endif
      
        write(iunit, *) 'Z_COORDINATES '
        write(iunit, *) obj % numCells(1) + 1, ' double'
        do i = 1, obj % numCells(1)
            write(iunit, *) ' ', 0.0 + obj % deltas(1) * i
        enddo

        write(iunit, *) 'CELL_DATA ', obj % ntot
        write(iunit, *) 'SCALARS f double 1'
        write(iunit, *) 'LOOKUP_TABLE default'
        do i = 1, obj % ntot
            write(iunit, *) obj % f(i), ' '
        enddo
    
        close(iunit)
    
    end subroutine

    subroutine upwind_print(obj)
            type(upwind_type) :: obj
        
        integer :: i
        
        do i = 1, obj % ntot
            print *, 'i = ', i, ' f = ',  obj % f(i)
        enddo
 
    end subroutine

    subroutine upwind_getIndexSet(obj, flatIndex, res) 
        type(upwind_type) :: obj
        integer, intent(in) :: flatIndex
        integer, intent(out) :: res(:)
    
        integer :: i
        do i = 1, obj % ndims
            res(i) = mod((flatIndex - 1)/obj % dimProd(i), obj % numCells(i)) + 1
        enddo
    end subroutine
    
    subroutine upwind_getFlatIndex(obj, inds, res) 
        type(upwind_type) :: obj
        integer, intent(in) :: inds(:)
        integer, intent(out) :: res
        
        res = dot_product(obj % dimProd, inds - 1) + 1
        
    end subroutine

end module

! Converts a string into an integer
subroutine str2int(str, i, stat)
    implicit none
    ! Arguments
    character(len=*), intent(in) :: str
    integer, intent(out)         :: i
    integer, intent(out)         :: stat

    read(str,*,iostat=stat) i
end subroutine str2int

program main

    use upwind_mod
    
    implicit none

    integer, parameter :: ndims = 3
    
    integer :: argc, numCells(ndims), n, ier, numTimeSteps, i, j
    character(len=32) :: argv
    real(r8) :: velocity(ndims)
    real(r8) :: lengths(ndims)
    real(r8) :: courant, dt, dx, val, chksum
    type(upwind_type) :: up

    numCells = -1
    
    ! default number of steps
    numTimeSteps = 100  
    argc = 0
    do 
        call get_command_argument(argc, argv)
        if (len_trim(argv) == 0) exit
        call str2int(argv, n, ier)
        if (argc == 1) then
            numCells = n
        else if (argc == 2) then
            numTimeSteps = n
        endif
        argc = argc + 1
    enddo

    if (argc < 2) then
        print *, 'must specify number of cells in each direction.'
        stop(1)
    endif

    print *, 'number of cells: ', numCells
    print *, 'number of time steps: ', numTimeSteps
    
    ! velocity field
    velocity = 1
    
    ! domain lengths
    lengths = 1

    ! compute time step from Courant's condition 
    courant = 0.1_r8
    dt = huge(1.0_r8)
    do j = 1, ndims
        dx = lengths(j) / numCells(j)
        val = courant * dx / velocity(j)
        dt = min(val, dt)
    enddo

    ! instantiate upwind object 
    call upwind_new(up, velocity, lengths, numCells)
    
    ! call upwind_saveVTK(up, 'up0.vtk')
    
    ! advance 
    do i = 1,numTimeSteps
        call upwind_advect(up, dt)
    enddo
    
    ! call upwind_print(up)
    print *, 'check sum: ', sum(up % f)
    call upwind_saveVTK(up, 'up.vtk')
    
    ! clean up
    call upwind_del(up)

end program
