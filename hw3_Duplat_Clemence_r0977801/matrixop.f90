!   Several implementations to compute the matrix-matrix product of 
!	two square (NxN), double precision matrices A and B.
!   For the indices i, j and k we use the following notation:
!   	C(i,j) = sum_{k=1}^{N} A(i,k)*B(k,j) for i,j = 1,..,N

module matrixop
    implicit none
    save
    integer, parameter :: dp = selected_real_kind(15,307) ! Get at least double precision
    integer :: threshold = 128
    interface
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            DOUBLE PRECISION ALPHA,BETA
            INTEGER K,LDA,LDB,LDC,M,N
            CHARACTER TRANSA,TRANSB
            DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
        end subroutine dgemm
    end interface
contains
    !--------------------------------------------------------------------------
    ! 1. Three nested loops
    !
    ! NOTE: use the following notation for the indices
    !    C(i,j) = sum_{k=1}^{N} A(i,k)*B(k,j)
    !       i = Row index of A and C
    !       j = Column index of B and C
    !       k = Column index of A and row index of B
    !--------------------------------------------------------------------------
    subroutine mm_ijk(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer i,j,k
        C = 0.0_dp 
        do i = 1,size(A,1)
            do j = 1,size(B,2)
                do k = 1,size(A,2)
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
            enddo     
        enddo
    end subroutine mm_ijk

    subroutine mm_ikj(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer i,k,j
        C = 0.0_dp 
        
        do i = 1, size(A,1)
            do k = 1, size(A,2)
                do j = 1, size(B,2)
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
            enddo     
        enddo
        ! TODO: complete this subroutine
    end subroutine mm_ikj

    subroutine mm_jik(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer i,j,k
        C = 0.0_dp
        do j = 1, size(B,2)
            do i = 1, size(A,1)
                do k = 1, size(A,2)
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
            enddo     
        enddo
    end subroutine mm_jik

    subroutine mm_jki(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer i,j,k
        C = 0.0_dp
        do j = 1, size(B,2)
            do k = 1, size(B,1)
                do i = 1, size(A,1)
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
            enddo     
        enddo
    end subroutine mm_jki
    

    subroutine mm_kij(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer i,j,k
        C = 0.0_dp
        do k = 1, size(A,2)
            do i = 1, size(A,1)
                do j = 1, size(B,2)
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
            enddo     
        enddo
    end subroutine mm_kij
    
    subroutine mm_kji(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer i,j,k
        C = 0.0_dp
        do k = 1, size(A,2)
            do j = 1, size(B,2)
                do i = 1, size(A,1)
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
            enddo     
        enddo
    end subroutine mm_kji
    !--------------------------------------------------------------------------
    ! 2. Two nested loops with vector operations
    !--------------------------------------------------------------------------
    subroutine mm_ikj_vect(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer i,k
        C = 0.0_dp 
        do i = 1, size(A,1)
            do k = 1, size(A,2)
                C(i,:) = C(i,:) + A(i,k)*B(k,:)
            enddo
        enddo
    end subroutine mm_ikj_vect

    subroutine mm_jki_vect(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B 
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer j,k
        C = 0.0_dp
        do j = 1, size(B,2)
            do k = 1, size(B,1)
                C(:,j) = C(:,j) + A(:,k)*B(k,j)
            enddo
        enddo
    end subroutine mm_jki_vect

    subroutine mm_kij_vect(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer k,i
        C = 0.0_dp
        do k = 1, size(A,2)
            do i = 1, size(A,1)
                C(i,:) = C(i,:) + A(i,k)*B(k,:)
            enddo
        enddo
    end subroutine mm_kij_vect

    subroutine mm_kji_vect(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer k,j
        C = 0.0_dp
        do k = 1, size(B,1)
            do j = 1, size(B,2)
                C(:,j) = C(:,j) + A(:,k)*B(k,j)
            enddo
        enddo
    end subroutine mm_kji_vect
    !--------------------------------------------------------------------------
    ! 3. Two nested loops with dot_product
    !--------------------------------------------------------------------------
    !---We should know RESULT = DOT_PRODUCT(VECTOR_A, VECTOR_B)----!
    subroutine mm_ijk_dot_product(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer i,j
        C = 0.0_dp
        do i = 1, size(A,1)
            do j = 1, size(B,2)
                C(i,j) = dot_product(A(i,:), B(:,j))
            enddo
        enddo
    end subroutine mm_ijk_dot_product

    subroutine mm_jik_dot_product(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
         ! TODO: complete this subroutine
        integer i,j
        C = 0.0_dp 
        do j = 1, size(B,2)
            do i = 1, size(A,1)
                C(i,j) = dot_product(A(i,:), B(:,j))
            enddo
        enddo
    end subroutine mm_jik_dot_product
    !--------------------------------------------------------------------------
    ! 4. Two nested loops with dot_product and explicit transpose of matrix A
    !--------------------------------------------------------------------------
    subroutine mm_transp_ijk_dot_product(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        !---use an a additional variable to store the transpose of A, inverse dimensions of A---!
        real(kind=dp), dimension(size(A,2),size(A,1)):: trans
         ! TODO: complete this subroutine
        integer i,j
        trans = transpose(A)
        C = 0.0_dp 
        do i = 1, size(trans,2)
            do j = 1, size(B,2)
                C(i,j) = dot_product(trans(:,i), B(:,j))
            enddo
        enddo
    end subroutine mm_transp_ijk_dot_product

    subroutine mm_transp_jik_dot_product(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        !---use an a additional variable to store the transpose of A, inverse dimensions of A---!
        real(kind=dp), dimension(size(A,2),size(A,1)):: trans
         ! TODO: complete this subroutine
        integer i,j 
        C = 0.0_dp
        trans=transpose(A)
        do j = 1, size(B,2)
            do i = 1, size(trans,2)
                C(i,j) = dot_product(trans(:,i), B(:,j))
            enddo
        enddo
    end subroutine mm_transp_jik_dot_product 
    !--------------------------------------------------------------------------
    ! 5. In blocks
    !--------------------------------------------------------------------------
    subroutine mm_blocks_a(A, B, C, blocksize)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer, intent(in) :: blocksize
        
        !---blocked variant of the fastest method with three nested loops---!
        !---assume that block size is a divisor of N---!
        !---use 6 nested loops, order of outer and inner loop should be the same---!
        !---mm_kij should iterate over k-i-j-k-i-j-----!
        integer i, j, k, ii, jj, kk
        ! TODO: Complete this subroutine
         
        !We will divide the input matrices into smaller blocks
        !Then multiply each block individually
        !----The choice of the loop order is important----!
        !----fastest: column by column----! -> in my case jki
        C = 0.0_dp 
        do jj = 1, size(A,2), blocksize
            do kk = 1, size(B,1), blocksize
                do ii = 1, size(A,1), blocksize
                    do j = jj, min(jj+blocksize-1, size(A,2))
                        do k = kk, min(kk+blocksize-1, size(B,1))
                            do i = ii, min(ii+blocksize-1, size(A,1))
                                C(i,j) = C(i,j) + A(i,k)*B(k,j)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine mm_blocks_a
    subroutine mm_blocks_b(A, B, C, blocksize)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer, intent(in) :: blocksize
        integer i, j, k, ii, jj, kk
        C = 0.0_dp ! TODO: Complete this subroutine
        !---blocked variant of the slowest method with three nested loops---!
        !---assume that block size is a divisor of N---!
        !---use 6 nested loops, order of outer and inner loop should be the same---!
        !----Slowest method----! -> in my case kij
        do kk = 1, size(A,2), blocksize
            do ii = 1, size(A,1), blocksize
                do jj = 1, size(B,2), blocksize
                    do k = kk, min(kk+blocksize-1, size(A,2))
                        do i = ii, min(ii+blocksize-1, size(A,1))
                            do j = jj, min(jj+blocksize-1, size(B,2))
                                C(i,j) = C(i,j) + A(i,k)*B(k,j)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine mm_blocks_b
    !--------------------------------------------------------------------------
    ! 6. Intrinsic matmul function
    !--------------------------------------------------------------------------
    subroutine mm_matmul(A, B, C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        C = matmul( A, B ) 
    end subroutine mm_matmul
    !--------------------------------------------------------------------------
    ! 7. Using BLAS
    !--------------------------------------------------------------------------
    subroutine mm_blas(A,B,C)
        real(kind=dp), dimension(:,:), intent(in)  :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer :: m, n, k
        real(kind=dp) :: alpha, beta
        C = 0.0_dp    ! TODO: Complete this subroutine
        !DGEMM  performs one of the matrix-matrix operations
        !C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
        !op( X ) = X   or   op( X ) = X**T,
        !alpha and beta are scalars, and A, B and C are matrices, with op( A )
        !an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
        
        !C=alpha*A*B +beta*C
        ! and here we want C=A*B so :
        alpha = 1
        beta = 0

        m = size(A, 1)  ! Number of rows of A and C
        n = size(B, 2)  ! Number of columns of B and C
        k = size(A, 2)  ! Number of columns of A and rows of B
        
        call dgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C, m)
    end subroutine mm_blas
    !--------------------------------------------------------------------------
    ! 8. Reference implementation to illustrate working principle of
    !    the divide and conquer method (for N = 2^k). Not optimized for 
    !    performance! Only for illustrative purposes. 
    !--------------------------------------------------------------------------
    recursive subroutine mm_divide_and_conquer(A, B, C)
        real(kind=dp), dimension(:,:), intent(in) :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer :: N,j,k
        real(kind=dp), dimension(size(A,1)/2,size(A,1)/2) :: M1,M2,M3,M4,M5,M6,M7,M8
        N = size(A,1)
        if (N > threshold) then
            call mm_divide_and_conquer(A(1:N/2,1:N/2),B(1:N/2,1:N/2),M1)
            call mm_divide_and_conquer(A(1:N/2,N/2+1:N),B(N/2+1:N,1:N/2),M2)
            call mm_divide_and_conquer(A(1:N/2,1:N/2),B(1:N/2,N/2+1:N),M3)
            call mm_divide_and_conquer(A(1:N/2,N/2+1:N),B(N/2+1:N,N/2+1:N),M4)
            call mm_divide_and_conquer(A(N/2+1:N,1:N/2),B(1:N/2,1:N/2),M5)
            call mm_divide_and_conquer(A(N/2+1:N,N/2+1:N),B(N/2+1:N,1:N/2),M6)
            call mm_divide_and_conquer(A(N/2+1:N,1:N/2),B(1:N/2,N/2+1:N),M7)
            call mm_divide_and_conquer(A(N/2+1:N,N/2+1:N),B(N/2+1:N,N/2+1:N),M8)
            C(1:N/2,1:N/2)     = M1 + M2
            C(1:N/2,N/2+1:N)   = M3 + M4
            C(N/2+1:N,1:N/2)   = M5 + M6
            C(N/2+1:N,N/2+1:N) = M7 + M8
        else 
            C = 0._dp 
            do j=1,N
                do k=1,N
                    C(:,j) = C(:,j) + A(:,k)*B(k,j)
                enddo
            enddo 
        endif
    end subroutine mm_divide_and_conquer
    !--------------------------------------------------------------------------
    ! 9. Reference implementation to illustrate working principle of 
    !    Strassen's algorithm (for N = 2^k). Not optimized for performance!
    !    Only for illustrative purposes.
    !--------------------------------------------------------------------------
    recursive subroutine mm_strassen(A, B, C)
        real(kind=dp), dimension(:,:), intent(in) :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer :: N,j,k
        real(kind=dp), dimension(size(A,1)/2,size(A,1)/2) :: M1,M2,M3,M4,M5,M6,M7
        N = size(A,1)
        if (N > threshold) then
            ! Compute M1
            call mm_strassen(A(1:N/2,1:N/2)+A(N/2+1:N,N/2+1:N), & 
                    B(1:N/2,1:N/2)+B(N/2+1:N,N/2+1:N),M1)
            ! Compute M2
            call mm_strassen(A(N/2+1:N,1:N/2)+A(N/2+1:N,N/2+1:N), &
                    B(1:N/2,1:N/2),M2)
            ! Compute M3
            call mm_strassen(A(1:N/2,1:N/2),B(1:N/2,N/2+1:N)- & 
                    B(N/2+1:N,N/2+1:N),M3)
            ! Compute M4
            call mm_strassen(A(N/2+1:N,N/2+1:N),&
                    B(N/2+1:N,1:N/2)-B(1:N/2,1:N/2),M4)
            ! Compute M5
            call mm_strassen(A(1:N/2,1:N/2)+A(1:N/2,N/2+1:N),& 
                    B(N/2+1:N,N/2+1:N),M5)
            ! Compute M6
            call mm_strassen(A(N/2+1:N,1:N/2)-A(1:N/2,1:N/2), & 
                    B(1:N/2,1:N/2)+B(1:N/2,N/2+1:N),M6)
            ! Compute M7
            call mm_strassen(A(1:N/2,N/2+1:N)-A(N/2+1:N,N/2+1:N), &
                 B(N/2+1:N,1:N/2)+B(N/2+1:N,N/2+1:N),M7)
            
            C(1:N/2,1:N/2) = M1 + M4 - M5 + M7
            C(1:N/2,N/2+1:N) = M3 + M5
            C(N/2+1:N,1:N/2) = M2 + M4 
            C(N/2+1:N,N/2+1:N) = M1 - M2 + M3 + M6
        else 
            C = 0._dp 
            do j=1,N
                do k=1,N
                    C(:,j) = C(:,j) + A(:,k)*B(k,j)
                enddo
            enddo 
        endif
    end subroutine mm_strassen
end module matrixop