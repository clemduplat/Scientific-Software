! Compares several implementations of the square matrix-matrix product in Fortran

program mm_driver
    use matrixop
    use timings
    implicit none
    
  
    !--------------------------------------------------------------------------
    ! Abstract interfaces
    !
    ! NOTE: this simplifies the timings.
    !--------------------------------------------------------------------------
    abstract interface
        subroutine mm_interface(a, b, c)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: a, b
            real(kind=dp), dimension(:,:), intent(out) :: c
        end subroutine mm_interface
        subroutine mm_blocks_interface(a,b,c,blocksize)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: a, b
            real(kind=dp), dimension(:,:), intent(out) :: c
            integer, intent(in) :: blocksize
        end subroutine mm_blocks_interface

    end interface

    !--------------------------------------------------------------------------
    ! Main timing program
    !--------------------------------------------------------------------------
    integer :: k
    integer, parameter :: N = 1600, blocksize = 100
    !integer, parameter :: N = 2000, blocksize = 500!100!20
    integer, dimension(:), allocatable :: seed
    real(kind=dp), dimension(N,N) :: a, b, c, c_matmul

    ! Make sure we use the same pseudo-random numbers each run by initializing
    ! the seed to a certain value.
    call random_seed(size=k)
    allocate(seed(k))
    seed = N
    call random_seed(put=seed)
   

    call random_number(a)
    call random_number(b)
    call mm_matmul(a,b,c_matmul) ! Reference value

    ! Start the timings
    print *, ""
    write(unit=*, fmt="(A)") "TIMING RESULTS:"
    
    ! 1. Three nested loops
    call do_timing( "IJK", mm_ijk )
    call do_timing( "IKJ", mm_ikj )
    call do_timing( "JIK", mm_jik )
    call do_timing( "JKI", mm_jki )
    call do_timing( "KIJ", mm_kij )
    call do_timing( "KJI", mm_kji )
    
   
    ! 2. Two nested loops with vector operations
    call do_timing( "IKJ, J VECT", mm_ikj_vect )
    call do_timing( "JKI, I VECT", mm_jki_vect )
    call do_timing( "KIJ, J VECT", mm_kij_vect )
    call do_timing( "KJI, I VECT", mm_kji_vect )
    
    ! 3. Two nested loops with dot_product
    call do_timing( "IJ DOT_PRODUCT", mm_ijk_dot_product )
    call do_timing( "JI DOT_PRODUCT", mm_jik_dot_product )
    
    ! 4. Two nested loops with dot_product and explicit transpose of matrix A
    call do_timing( "IJ TP DOT_PRODUCT", mm_transp_ijk_dot_product )
    call do_timing( "JI TP DOT_PRODUCT", mm_transp_jik_dot_product ) 
    
    ! 5. In blocks
    call do_timing( "IN BLOCKS A", method_blocks=mm_blocks_a )
    call do_timing( "IN BLOCKS B", method_blocks=mm_blocks_b )

    ! 6. Intrinsic matmul function
    call do_timing( "MATMUL", mm_matmul )

    ! 7. Divide and conquer
    call do_timing( "DIVIDE AND CONQUER", mm_divide_and_conquer )

    ! 8. Strassen algorithm
    call do_timing( "STRASSEN ALGORITHM", mm_strassen )

    ! 9. Using BLAS
    call do_timing( "BLAS", mm_blas )
    


    ! Clean up
    deallocate(seed)
contains

    real(dp) function relError( )
        ! Compute the Frobenius norm of the error with respect to the 
        ! result returned by matmul and divide by the Frobenius norm 
        ! of the result from matmul to get some kind of relative error.
        relError = sqrt(sum((C_matmul-C)**2))/sqrt(sum(C_matmul**2))
    end function relError

    subroutine do_timing(funcName,method, method_blocks)
        character(len=*), intent(in) :: funcName
        procedure(mm_interface), optional :: method
        procedure(mm_blocks_interface), optional :: method_blocks
        real(kind(0.e0)) :: elapsedTime
        ! Do the timing
        call tic()
        if( present(method) ) then
            call method( a, b, c )
        else
            call method_blocks( a, b, c, blocksize)
        end if
        call toc(elapsedTime)
        print "(A18, A, F10.4, A, ES9.2)", funcName, ": CPU time", elapsedTime, " sec, relative error = ", relError()        

    end subroutine do_timing

end program mm_driver