#ifndef _THREADS_PER_TEAM
#define _THREADS_PER_TEAM 32
#endif
#ifndef _TEAMS
#define _TEAMS 64
#endif
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
        real(r8), pointer :: f(:)
        
        ! product of the dimensions, used to switch back and forth 
        ! between the flat index and the multi-index representations
        integer, allocatable :: dimProd(:)

        contains 

            procedure :: new => upwind_new
            procedure :: del => upwind_del
            procedure :: advect => upwind_advect
            procedure :: saveVTK => upwind_saveVTK
            procedure :: getFlatIndex => upwind_getFlatIndex

    end type

contains

    subroutine getIndexSet(i, ndims, numCells, dimProd, inds)
        implicit none
        integer, intent(in) :: i
        integer, intent(in) :: ndims
        integer, intent(in) :: numCells(:)
        integer, intent(in) :: dimProd(:)
        integer, intent(out) :: inds(:)

#if defined(HAVE_OPENACC)
        !$acc routine
#else
        !$omp declare target
#endif

        integer :: j

        do j = 1, ndims
            inds(j) = mod((i - 1)/dimProd(j), numCells(j)) + 1
        enddo
    end subroutine getIndexSet

    ! Constructor
    ! @param velocity velocity field (constant)
    ! @param lengths domain lengths
    ! @param numCells number of cells in the x, y, ... directions
    subroutine upwind_new(this, velocity, lengths, numCells)
    
        class(upwind_type) :: this
        real(r8), intent(in) :: velocity(:)
        real(r8), intent(in) :: lengths(:)
        integer, intent(in) :: numCells(:)
        
        integer :: j
        
        this % ndims = size(velocity)
        
        allocate(this % upDirection(this % ndims))
        allocate(this % deltas(this % ndims))
        allocate(this % v(this % ndims))
        allocate(this % lengths(this % ndims))
        allocate(this % numCells(this % ndims))
        allocate(this % dimProd(this % ndims))
        
        this % v = velocity
        this % lengths = lengths
        this % numCells = numCells

        ! compute the total number of cells and other stuff
        this % ntot = 1       
        do j = 1, this % ndims
            this % upDirection(j) = -1
            if (velocity(j) < 0.) then
                this % upDirection(j) = 1
            endif
            this % deltas(j) = lengths(j) / numCells(j)
            this % ntot = this % ntot * numCells(j)
        enddo
        
        this % dimProd(this % ndims) = 1
        do j = this % ndims - 1, 1, -1
            this % dimProd(j) =  this % dimProd(j + 1) * this % numCells(j + 1)
        enddo
        
        allocate(this % f(this % ntot))
        
        ! initialize the field, zero everywhere except for the 
        ! low corner cell where the field is one
        this % f = 0
        this % f(1) = 1

    end subroutine 
    
    ! Destructor
    subroutine upwind_del(this)

        class(upwind_type) :: this
        
        deallocate(this % v)
        deallocate(this % lengths)
        deallocate(this % numCells)    
        deallocate(this % upDirection)
        deallocate(this % deltas)
        deallocate(this % dimProd)
        deallocate(this % f)
        
    end subroutine
    
    ! Advance by one time step
    ! @param deltaTime time step
    subroutine upwind_advect(this, deltaTime)
    
        class(upwind_type) :: this
        real(r8), intent(in) :: deltaTime

        real(r8), allocatable :: oldF(:)
        integer :: i, j, oldIndex, upI
        integer :: inds(this % ndims)
        real(r8), pointer :: fptr(:)
        real(r8) :: v(this % ndims), deltas(this % ndims)
        integer :: ntot, ndims, numCells(this % ndims), dimProd(this % ndims)
        integer :: upDirection(this % ndims)
        
        ! allocate and copy the field
        allocate(oldF(this % ntot))
        
        fptr => this % f
        ntot = this % ntot
        ndims = this % ndims
        do j = 1, ndims
            numCells(j) = this % numCells(j)
            upDirection(j) = this % upDirection(j)
            dimProd(j) = this % dimProd(j)
            v(j) = this % v(j)
            deltas(j) = this % deltas(j)
        enddo

#if defined(HAVE_OPENACC)
    #if defined(HAVE_OPENACC_KERNELS)
        #warning "Using OpenACC kernels directive"
        !$acc kernels
    #else
        !$acc parallel loop 
    #endif
#elif defined(HAVE_OPENMP_OFFLOAD)
        !$omp target
        !$omp teams distribute parallel do thread_limit(_THREADS_PER_TEAM)
#elif defined(HAVE_OPENMP)
        !$omp parallel do
#endif
        do i = 1, ntot
            oldF(i) = fptr(i)
        enddo
#if defined(HAVE_OPENACC)
    #if defined(HAVE_OPENACC_KERNELS)
        !$acc end kernels
    #else
        !$acc end parallel loop
    #endif
#elif defined(HAVE_OPENMP_OFFLOAD)
        !$omp end teams distribute parallel do
        !$omp end target
#elif defined(HAVE_OPENMP)
        !$omp end parallel do
#endif

        ! iterate over the cells
#if defined(HAVE_OPENACC)
    #if defined(HAVE_OPENACC_KERNELS)
        #warning "Using OpenACC kernels directive"
        !$acc kernels
        !$acc loop independent, private(inds, oldIndex, upI, i, j)
    #else
        !$acc parallel loop private(inds, j, oldIndex, upI)
    #endif
#elif defined(HAVE_OPENMP_OFFLOAD)
        !$omp target
        !$omp teams distribute parallel do private(inds, j, oldIndex, upI) thread_limit(_THREADS_PER_TEAM)
#elif defined(HAVE_OPENMP)
        !$omp parallel do private(inds, j, oldIndex, upI)
#endif
        do i = 1, ntot

            ! compute the index set of this cell
            call getIndexSet(i, ndims, numCells, dimProd, inds)
            !do j = 1, ndims
            !    inds(j) = mod((i - 1)/dimProd(j), numCells(j)) + 1
            !enddo

            do j = 1, ndims
                    
                ! cache the cell index
                oldIndex = inds(j)
                
                ! increment the cell index
                inds(j) = inds(j) + upDirection(j)
                
                ! apply periodic BCs
                inds(j) = modulo(inds(j) + numCells(j) - 1, numCells(j)) + 1
                  
                ! compute the new flat index 
                upI = dot_product(dimProd, inds - 1) + 1
                
                ! update the field
                fptr(i) = fptr(i) - &
                  &   deltaTime*v(j)*upDirection(j)*(oldF(upI) - oldF(i))/deltas(j)
                    
                ! reset the index
                inds(j) = oldIndex
            enddo
        enddo
#if defined(HAVE_OPENACC)
    #if defined(HAVE_OPENACC_KERNELS)
        !$acc end kernels
    #else
        !$acc end parallel loop
    #endif
#elif defined(HAVE_OPENMP_OFFLOAD)  
        !$omp end teams distribute parallel do
        !$omp end target
#elif defined(HAVE_OPENMP)
        !$omp end parallel do
#endif
        deallocate(oldF)

    end subroutine

    subroutine upwind_saveVTK(this, filename)
        class(upwind_type) :: this
        character(len=*), intent(in) :: filename
        
        integer iunit, i
        
        ! f2008
        !open(newunit = iunit, file = filename, status = 'unknown')
        iunit = 10
        ! f95
        open(unit = iunit, file = filename, status = 'unknown')
        write(iunit, '(a)') '# vtk DataFile Version 2.0'
        write(iunit, '(a)') 'upwind.f90'
        write(iunit, '(a)') 'ASCII'
        write(iunit, '(a)') 'DATASET RECTILINEAR_GRID'

        ! in VTK the first dimension varies fastest so need 
        ! to invert the order of the dimensions
        if (this % ndims > 2) then
            write(iunit, '(a, i10, i10, i10)') 'DIMENSIONS ', &
             & this % numCells(3) + 1, this % numCells(2) + 1, this % numCells(1) + 1
        else
            if (this % ndims > 1) then
                write(iunit, '(a, i10, i10)') 'DIMENSIONS 1', &
             & this % numCells(2) + 1, this % numCells(1) + 1
            else
                write(iunit, '(a, i10)') 'DIMENSIONS 1 1', this % numCells(1) + 1
            endif 
        endif 
        
        write(iunit, '(a, i10, a)') 'X_COORDINATES ', this % numCells(1) + 1, ' double'
        do i = 1, this % numCells(1) + 1
            write(iunit, '(e20.7)') 0.0 + this % deltas(1) * (i - 1)
        enddo

        write(iunit, *) 
        if (this % ndims > 1) then
            write(iunit, '(a, i10, a)') 'Y_COORDINATES ', this % numCells(2) + 1, ' double'
            do i = 1, this % numCells(2) + 1
                write(iunit, '(e20.7)') 0.0 + this % deltas(2) * (i - 1)
            enddo      
        else
            write(iunit, '(a)') 'Y_COORDINATES 1 double'
            write(iunit, '(a)') '0.0'
        endif

        write(iunit, *) 
        if (this % ndims > 2) then
            write(iunit, '(a, i10, a)') 'Z_COORDINATES ', this % numCells(3) + 1, ' double'
            do i = 1, this % numCells(3) + 1
                write(iunit, '(e20.7)') 0.0 + this % deltas(3) * (i - 1)
            enddo      
        else
            write(iunit, '(a)') 'Z_COORDINATES 1 double'
            write(iunit, '(a)') '0.0'
        endif

        write(iunit, '(a, i20)') 'CELL_DATA ', this % ntot
        write(iunit, '(a)') 'SCALARS f double 1'
        write(iunit, '(a)') 'LOOKUP_TABLE default'
        do i = 1, this % ntot
            write(iunit, '(e20.7)') this % f(i)
        enddo
    
        close(iunit)
    
    end subroutine

    subroutine upwind_print(this)
        class(upwind_type), intent(in) :: this
        integer :: i
        
        do i = 1, this % ntot
            write(*, '(a, i10, a, e20.13)') 'i = ', i, ' f = ',  this % f(i)
        enddo
 
    end subroutine
    
    pure function upwind_getFlatIndex(this, inds) result(res) 
        class(upwind_type), intent(in) :: this
        integer, intent(in) :: inds(:)
        integer :: res
        
        res = dot_product(this % dimProd, inds - 1) + 1
        
    end function

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
#ifdef HAVE_OPENMP
    use omp_lib
#endif
    implicit none

    integer, parameter :: ndims = 3
    
    integer :: argc, numCells(ndims), n, ier, numTimeSteps, i, j, & 
   &           numThreads, maxNumThreads, threadId
    logical :: doVtk
    character(len=32) :: argv
    real(r8) :: velocity(ndims)
    real(r8) :: lengths(ndims)
    real(r8) :: courant, dt, dx, val, chksum
    type(upwind_type) :: up
#ifdef HAVE_OPENMP
    !$omp parallel
    numThreads = omp_get_num_threads()
    !$omp end parallel
    maxNumThreads = omp_get_max_threads()
    write(*, '(A,I4,A,I4)') 'Number of OpenMP threads: ', numThreads, ' out of a max number: ', maxNumThreads
#endif
 
    numCells = -1
    doVtk = .FALSE.
    
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
        else if (argc == 3 .and. argv == 'vtk') then
            doVtk = .TRUE.
        endif
        argc = argc + 1
    enddo

    if (argc < 2) then
        stop 'must specify number of cells in each direction.'
    endif

    write(*, '(a)', advance='no') 'number of cells: '
    do i = 1, ndims
        write(*, '(i10, a)', advance='no') numCells(i), ' '
    enddo
    write(*, *) ' ' ! new line
    write(*, '(a,i10)') 'number of time steps: ', numTimeSteps
    
    ! velocity field
    velocity = 1
    
    ! domain lengths
    lengths = 1

    ! compute time step from Courant's condition 
    courant = 0.1_r8
    dt = huge(1.0_r8)
    do j = 1, ndims
        dx = lengths(j) / real(numCells(j), r8)
        val = courant * dx / velocity(j)
        dt = min(val, dt)
    enddo

    ! instantiate up
    call up % new(velocity, lengths, numCells)
    
    ! call up % saveVTK('up0.vtk')
    
    ! advance 
    do i = 1, numTimeSteps
        call up % advect(dt)
    enddo
    
    write(*,'(a, f15.9)') 'check sum: ', sum(up % f)

    if (doVtk) then 
        call up % saveVTK('up.vtk')
    endif
    
    ! clean up
    call up % del()

end program
