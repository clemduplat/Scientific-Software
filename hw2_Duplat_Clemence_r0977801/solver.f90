module solver
    ! Write your own solver module using sgesv and dgesv to solve in single an double precision
    ! Being able to handle both single and double precision without recompilation -> generic interface
    !----generic interface----!
    ! automatically select the appropriate procedure based on the arguments
    implicit none
    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)
    interface solve
        module procedure solve_single, solve_double
    end interface solve
contains
    !-----Explanation Lapack Library used----!
    !SGESV computes the solution to a real system of linear equations
    !    A * X = B,
    !where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
    !  N is INTEGER
    !      The number of linear equations, i.e., the order of the
    !      matrix A.  N >= 0.
    !NRHS is INTEGER
    !        The number of right hand sides, i.e., the number of columns
    !        of the matrix B.  NRHS >= 0.
    !IPIV is INTEGER array, dimension (N)
    !      The pivot indices that define the permutation matrix P;
    !      row i of the matrix was interchanged with row IPIV(i).
    ! INFO is INTEGER
    !      = 0:  successful exit
    !      < 0:  if INFO = -i, the i-th argument had an illegal value
    !      > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
    !            has been completed, but the factor U is exactly
    !            singular, so the solution could not be computed.
    subroutine solve_single(A,bx)
        real(sp),dimension(:,:), intent(inout) :: A  
        real(sp), dimension(:), intent(inout) :: bx   
        integer, dimension(:), allocatable :: ipiv 
        integer :: n, nrhs, info
        n = size(A, 1)       
        nrhs = 1  ! ex: number of colums of matrix b-> 1 in our case           
        allocate(ipiv(n))   ! integer array dimension N 

        !------Solve the equation using LAPACK routine sgesv-----!
        call sgesv(n, nrhs, A, n, ipiv, bx, n, info)

        !-------Check for a successful solution-----!
        if (info /= 0) then
            if (info < 0) then
                print *, 'Argument ', -info, ' had an illegal value.'
            else
                print *, 'The ', info, 'th diagonal element of the triangular factor is zero.'
            endif
        endif

        deallocate(ipiv)
    end subroutine solve_single
    subroutine solve_double(A,bx)
        real(dp),dimension(:,:), intent(inout) :: A  
        real(dp), dimension(:), intent(inout) :: bx   
        integer, dimension(:), allocatable :: ipiv 
        integer :: n, nrhs, info
        n = size(A, 1)       
        nrhs = 1             
        allocate(ipiv(n))   

        !------Solve the equation using LAPACK routine dgesv-----!
        call dgesv(n, nrhs, A, n, ipiv, bx, n, info)

        !-------Check for a successful solution-----!
        if (info /= 0) then
            if (info < 0) then
                print *, 'Argument ', -info, ' had an illegal value.'
            else
                print *, 'The ', info, 'th diagonal element of the triangular factor is zero.'
            endif
        endif

        deallocate(ipiv)
    end subroutine solve_double
end module solver