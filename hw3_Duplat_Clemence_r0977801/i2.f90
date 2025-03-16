program i2
    use matrixop
    use timings
    implicit none
    abstract interface
        subroutine mm_interface(A, B, C)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: A, B
            real(kind=dp), dimension(:,:), intent(out) :: C
        end subroutine mm_interface
        subroutine mm_blocks_interface(A,B,C,blocksize)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: A, B
            real(kind=dp), dimension(:,:), intent(out) :: C
            integer, intent(in) :: blocksize
        end subroutine mm_blocks_interface

    end interface
    !-----Pointer to the slowest and fastest method---!
    !-> make it easier to change of method
    procedure(mm_interface), pointer :: slowest_method_ptr, fastest_method_ptr
    !---Defining slowest and fastest method-----!
    
    procedure(mm_blocks_interface), pointer :: blocked_slowest_method_ptr, blocked_fastest_method_ptr
    

    integer, parameter :: N_warmup = 2000 
    integer, parameter :: block_size = 100
    real(kind=dp), allocatable, dimension(:,:) :: A_warmup,&
         B_warmup,C_warmup
    real(kind=dp), allocatable, dimension(:,:) :: A, B, C
    integer :: statusAllocate
    integer :: N
    real(kind(0.e0)) :: elapsedTimeSlowest, elapsedTimeFastest
    real(kind=dp) ::  mflopsSlowest, mflopsFastest
    ! Compute the matrix matrix multiplication of two large matrices to avoid warmup effects
    allocate(A_warmup(N_warmup,N_warmup),B_warmup(N_warmup,N_warmup),C_warmup(N_warmup,N_warmup),STAT=statusAllocate)
    if (statusAllocate>0) then
        print *, "Failed to allocate warmup matrices. Terminating program."
        call exit;
    endif
    call random_number(A_warmup)
    call random_number(B_warmup)
    C_warmup = matmul(A_warmup,B_warmup)
    C_warmup = matmul(A_warmup,B_warmup)
    deallocate(A_warmup,B_warmup,C_warmup)
    slowest_method_ptr => mm_kij !30sec
    fastest_method_ptr => mm_jki !2sec
    blocked_slowest_method_ptr => mm_blocks_b
    blocked_fastest_method_ptr => mm_blocks_a
    
    ! Add your implementation here. 
    ! Do not forget to call random_number for each matrix size.
    OPEN(10,file="i2.dat")
    do N = 10,100,10
        !----Allocating matrices and verify the status----!
        allocate(A(N, N), B(N, N), C(N, N),STAT=statusAllocate)
        if (statusAllocate>0) then
            print *, "Failed to allocate A,B,C matrices. Terminating program."
            call exit;
        endif

        !----Initialize with random numbers-----!
        call random_number(A)
        call random_number(B)
        C=0.0_dp
        !--------slowest method------!
        call time_method(slowest_method_ptr, A, B, C, elapsedTimeSlowest)
        !--------fastest method--------!
        call time_method(fastest_method_ptr, A, B, C, elapsedTimeFastest)
        !-----Flop calculation------!
        !----Avoid having negative values-----!
        if ( mod(N,10) == 0) then
            mflopsSlowest = 2*(N/10)**3 /(1000._dp*elapsedTimeSlowest)
            mflopsFastest = 2*(N/10)**3 /(1000._dp*elapsedTimeFastest)
        else
            mflopsSlowest = 2*(real(N,kind=dp)/100._dp)**3 /elapsedTimeSlowest
            mflopsFastest = 2*(real(N,kind=dp)/100._dp)**3 /elapsedTimeFastest
        end if
        !mflopsSlowest = (2.0_dp * N**3) / (elapsedTimeSlowest * 1e6_dp)
        !mflopsFastest = (2.0_dp * N**3) / (elapsedTimeFastest * 1e6_dp)

        !----------Output format---------------!
        !write(10,'(i0,2("  ",f12.4))') N, mflopsSlowest, mflopsFastest
        

        deallocate(A, B, C)
    end do
   
    do N = 100,1600,100
        
        allocate(A(N, N), B(N, N), C(N, N),STAT=statusAllocate)
        if (statusAllocate>0) then
            print *, "Failed to allocate A,B,C matrices. Terminating program."
            call exit;
        endif
        
        ! Initialize matrices with random numbers
        call random_number(A)
        call random_number(B)
        !C=0.0_dp
        !--------slowest method------!

        !call time_method(slowest_method_ptr, A, B, C, elapsedTimeSlowest)
        !mflopsSlowest = (2.0_dp * N**3) / (elapsedTimeSlowest * 1e6_dp)
        !--------fastest method--------!
        !call time_method(fastest_method_ptr, A, B, C, elapsedTimeFastest)
        !-------------Using blocks-------------------!
        call time_method_block(blocked_slowest_method_ptr,A,B,C,block_size,elapsedTimeSlowest)
        call time_method_block(blocked_fastest_method_ptr,A,B,C,block_size,elapsedTimeFastest)
        !-----Flop calculation------!
        !----Avoid having negative values-----!
        if ( mod(N,10) == 0) then
            mflopsFastest = 2*(N/100)**3 /elapsedTimeFastest
            mflopsSlowest = 2*(N/100)**3 /elapsedTimeSlowest
        else
            mflopsFastest = 2*(real(N,kind=dp)/100._dp)**3 /elapsedTimeFastest
            mflopsSlowest = 2*(real(N,kind=dp)/100._dp)**3 /elapsedTimeSlowest
        end if
        !mflopsFastest = (2.0_dp * N**3) / (elapsedTimeFastest * 1e6_dp)

        !----------Output format---------------!
        
        !write(10,'(i0,2("  ",f12.4))') N, mflopsSlowest, mflopsFastest
        !---For block size----!
        write(10,'(i0,3("  ",f12.4))') N, mflopsSlowest, mflopsFastest
        deallocate(A, B, C)
    end do

    
   
    

    CLOSE(10)
contains
    subroutine time_method(mm_method, A, B, C, elapsedTime)
        procedure(mm_interface), pointer :: mm_method
        real(kind=dp), dimension(:,:), intent(in) :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        real(kind(0.e0)), intent(out) :: elapsedTime

        call tic()
        call mm_method(A, B, C)
        call toc(elapsedTime)
    end subroutine time_method
    subroutine time_method_block(mm_method, A, B, C, blocksize, elapsedTime)
        procedure(mm_blocks_interface), pointer :: mm_method
        real(kind=dp), dimension(:,:), intent(in) :: A, B
        real(kind=dp), dimension(:,:), intent(out) :: C
        integer, intent(in) :: blocksize
        real(kind(0.e0)), intent(out) :: elapsedTime
    
        call tic()
        call mm_method(A, B, C, blocksize)
        call toc(elapsedTime)
    end subroutine time_method_block
    

end program i2