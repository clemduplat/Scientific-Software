program blas_divide_and_conquer
    integer, parameter :: dp = selected_real_kind(15,307), N = 2000
    integer, dimension(:), allocatable :: seed
    real(kind=dp), dimension(N,N) :: a, b, c

    ! Make sure we use the same pseudo-random numbers each run by initializing
    ! the seed to a certain value.
    call random_seed(size=k)
    allocate(seed(k))
    seed = N
    call random_seed(put=seed)

    call random_number(A)
    call random_number(B)

    call divide_and_conquer(A,B,C)    
    
    deallocate(seed)

contains

subroutine divide_and_conquer(A,B,C)
    ! One step of divide and conquer, uses a blas routine to multiply 
    ! the submatrices. Assumes that the matrix size N divisible by 2.

    ! Do not modify these lines
    real(kind=dp), dimension(:,:), intent(in) :: A, B
    real(kind=dp), dimension(:,:), intent(out):: C
    character,parameter :: transpA = 'N', transpB = 'N'
    integer :: m, demi_m
    m = size(A,1)
    demi_m = m / 2

    ! Modify these lines: choose the appropriate blas routine, and values for M,N,K,alpha,Axx,
    ! LDA,Bxx,LDB,beta,Cxx and LDC.
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(1,1), m, B(1,1), m, 0.0_dp, C(1,1), m) ! C11 = A11*B11
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(1:demi_m,demi_m+1:m), m, B(demi_m+1:m,1:demi_m), m, 1.0_dp, C(1:demi_m,1:demi_m), m) ! C11 += A12*B21
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(1:demi_m,1:demi_m), m, B(1:demi_m,demi_m+1:m), m, 0.0_dp, C(1:demi_m,demi_m+1:m), m) ! C12 = A11*B12
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(1:demi_m,demi_m+1:m), m, B(demi_m+1:m,demi_m+1:m), m, 1.0_dp, C(1:demi_m,demi_m+1:m), m) ! C12 += A12*B22
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(demi_m+1:m,1:demi_m), m, B(1:demi_m,1:demi_m), m, 0.0_dp, C(demi_m+1:m,1:demi_m), m) ! C21 = A21*B11
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(demi_m+1:m,demi_m+1:m), m, B(demi_m+1:m,1:demi_m), m, 1.0_dp, C(demi_m+1:m,1:demi_m), m) ! C21 += A22*B21
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(demi_m+1:m,1:demi_m), m, B(1:demi_m,demi_m+1:m), m, 0.0_dp, C(demi_m+1:m,demi_m+1:m), m) ! C22 = A21*B12
    call dgemm(transpA, transpB, demi_m, demi_m, demi_m, 1.0_dp, A(demi_m+1:m,demi_m+1:m), m, B(demi_m+1:m,demi_m+1:m), m, 1.0_dp, C(demi_m+1:m,demi_m+1:m), m) ! C22 += A22*B22
end subroutine divide_and_conquer    
    
end program blas_divide_and_conquer