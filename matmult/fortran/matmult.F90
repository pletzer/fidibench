module matmult_mod

    contains

    ! matrix initialization
    subroutine initmat(a, b)

        implicit none

        real(8), intent(out) :: a(:, :)
        real(8), intent(out) :: b(:, :)
        real(8) :: count
        integer :: i, j, k, n, m, el
        n = size(a, 1)
        m = size(b, 1)
        el = size(b, 2)
        count = 1
        do j = 1, m
            do i = 1, n
                a(i, j) = count
                count = count + 1
            enddo
        enddo
        do k  = 1, el
            do j = 1, m
                b(j, k) = count
                count = count + 1
            enddo
        enddo
    
    end subroutine initmat

    ! naive implementation
    subroutine matmult0(a, b, c)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(in) :: b(:, :)
        real(8), intent(out) :: c(:, :)
        integer :: n, m, el
        integer :: i, j, k
        n = size(a, 1)
        m = size(b, 1)
        el = size(c, 2)
        do k = 1, el
            do i = 1, n
                c(i, k) = 0
                do j = 1, m
                    c(i, k) = c(i, k) + a(i, j)*b(j, k)
                enddo
            enddo
        enddo
    end subroutine matmult0

    ! parallelized
    subroutine matmult1(a, b, c)
        implicit none
        real(8), intent(in) :: a(:, :)
        real(8), intent(in) :: b(:, :)
        real(8), intent(out) :: c(:, :)
        integer :: n, m, el
        integer :: i, j, k
        n = size(a, 1)
        m = size(b, 1)
        el = size(c, 2)
        do concurrent (k = 1:el)
            do i = 1, n
                c(i, k) = 0
                do j = 1, m
                    c(i, k) = c(i, k) + a(i, j)*b(j, k)
                enddo
            enddo
        enddo
    end subroutine matmult1


end module matmult_mod

program test
    use matmult_mod
    implicit none
    real(8), allocatable :: a(:, :), b(:, :)
    real(8), allocatable :: c0(:, :), c1(:, :)
    integer :: n, m, el
    integer :: argc, ier, val
    character(len=32) :: argv
    real(8) :: count, diff

    ! default matrix sizes
    n = 3
    m = 4
    el = 5

    ! parse the command line arguments
    argc = 1
    do 
        call get_command_argument(argc, argv)
        if (len_trim(argv) == 0) exit

        ! covnert string to integer
        read(argv, *) val

        if (ier == 0) then

            ! is an integer
            if (argc == 1) then 
                n = val
            else if (argc == 2) then 
                m = val
            else if (argc == 3) then 
                el = val
            endif
        endif

        argc = argc + 1
    enddo

    print *, 'matrix sizes: ', n, ' x', m, ' x', el
    allocate (a(n, m))
    allocate (b(m, el))
    allocate (c0(n, el))
    allocate (c1(n, el))

    call initmat(a, b)

    call matmult0(a, b, c0)
    call matmult1(a, b, c1)

    diff = sum(abs(c0 - c1))
    if (diff > 1.e-10) then
        print *,'ERROR: signifiant difference diff=', diff
        stop
    endif
    print *,'SUCCESS.'

end program test