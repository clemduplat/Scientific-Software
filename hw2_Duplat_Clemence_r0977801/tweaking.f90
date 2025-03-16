
program tweaking
    use timing
    use solvers
    implicit none
    ! Fortran program wich expects N,T and delta(beta) and target maximal number ofinfections as command line argument
    ! Read a given initial state model parameter as last assignement
    ! Solve problem F(beta)=target wich returns maximum F_betaue of I_k+Q_k over all k, using Newton's method
    ! Compute F(beta) using a finite difference approximation
    ! Pick a suitable stop criterion for the iteration and stop your choice
    ! write code in a modular way, easy to switch between simulation methods
    !------Declare the variables------! 
    
    real(KREAL):: tol
    real(KREAL) :: optimal_beta,initial_beta
    integer  :: max_iter = 1000
    real(KREAL) :: F_beta_prev
    real(KREAL) ::  update_step
    real(KREAL) :: max_step
    integer :: a
    integer :: b
    

    
    call get_parameters_tweaking
    write(*,*) "Give range of value i for Delta(beta)=10^i:"
    read(*,*) b
    dt = T / real(N, KREAL)

    !call cpu_time(start_time)

    call get_method
    initial_beta=beta
    
    do a=1,b
        call tic()
        F_beta_prev = 1.0e30
        deltaBeta = 10.0_KREAL**(-a)
        beta = initial_beta
        tol=1.0E-6_KREAL
        max_step=0.001_KREAL
        do iteration = 1, max_iter
            !-------Compute F_beta-------!
            F_beta = compute_F_beta(beta)
            
            
            !------Implementing my stopping criterium----!
            if (abs(F_beta - target) < abs(target) * tol) then
                !print*,"Condition1"
                exit
            endif

            
            !----Derivative using difference approximation----!
            lastBeta = beta
            beta = beta + deltaBeta
            F_beta_plus_delta = compute_F_beta(beta)
            dF = (F_beta_plus_delta - F_beta) / deltaBeta
            
            !------Checking if dF is not too small for Newton method-------!
            
            F_beta_prev=F_beta
            
            
            update_step=-(F_beta-target)/dF
            if(abs(dF)<tol) then
                !-----Manually correct the update_step----!
                update_step=min(max(-max_step,update_step),max_step)
                
            endif
            
            !----------Newton Method--------!
            beta=lastBeta+update_step
            !if(abs(beta-lastBeta)<tol*lastBeta) exit
            
            
            
        end do
        
        optimal_beta=beta
        
        !call cpu_time(end_time)
        !total_time = end_time - start_time
        !write(*, '(F20.15)') deltaBeta
        !write(*, '(F20.15)') optimal_beta
        
        !write(*, '(F20.15)') toc()
        !write(*, *) iteration
        write(*, '(I10,1X,F20.15,1X,F20.15,1X,I10)') a, beta, toc(), iteration
        
    end do
    
    
    !------- Dont forget to deallocate arrays-----------!
    
contains
    function compute_F_beta(beta) result(result_value)
        real(KREAL), intent(in) :: beta
        real(KREAL) :: result_value
        real(KREAL) :: total_population
        allocate( S(1:N+1), I(1:N+1), Q(1:N+1), R(1:N+1), D(1:N+1))
        total_population=11e6
        S(1) =  11e6-5.0_KREAL
        I(1) = 5.0_KREAL
        Q(1) = 0.0_KREAL
        R(1) = 0.0_KREAL
        D(1)  =0.0_KREAL
        
        result_value=0.0_KREAL
        do k=1, N 
            SELECT CASE (method_number)
            CASE (1)
                call solver_forward(S(k), I(k), Q(k), R(k), D(k), S(k+1), I(k+1), Q(k+1), R(k+1), D(k+1))
            CASE (2)
                call solver_heun(S(k), I(k), Q(k), R(k), D(k), S(k+1), I(k+1), Q(k+1), R(k+1), D(k+1))
            CASE (3)
                call solver_backward(S(k), I(k), Q(k), R(k), D(k), S(k+1), I(k+1), Q(k+1), R(k+1), D(k+1))
            CASE DEFAULT
                ! Handle the case when method_number doesn't match any of the above cases
                ! You can add error handling or other logic here if needed
                print*,"Wrong method"
            END SELECT
            
            result_value=max(result_value, I(k)+Q(k))
            
        end do
        !print*, result_value
        deallocate( S, I, Q, R, D)
    end function compute_F_beta
    



end program tweaking