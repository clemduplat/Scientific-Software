program callSolver
    use solvers
    
    implicit none
    
    
    real(KREAL) ::  somme
    
    
    
    !--------- Get the parameters--------!
    call get_parameters_file
    !call get_parameters_question
    
    !-------------Allocate arrays and put initial conditions-----------!
    dt = T / real(N, KREAL)
    !allocate(t_arr(1:N+1), S(1:N+1), I(1:N+1), Q(1:N+1), R(1:N+1), D(1:N+1))
    allocate(t_arr(1:N+1), S(1:N+1), I(1:N+1), Q(1:N+1), R(1:N+1), D(1:N+1))
    call get_method
    t_arr(1) = 0.0
    S(1) = S0
    I(1) = I0
    Q(1) = 0.0
    R(1) = 0.0
    D(1) = 0.0
    do k = 1, N
        t_arr(k+1) = t_arr(k) + dt
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
    end do

    
    
    
    somme=S(k)+I(k)+Q(k)+R(k)+D(k)
    
    print*,"La somme doit être égale à 105 : f ", somme
    !----------------Output the results in a file so that I have : t S I Q R D at each iteration---------!
    !---Forward---!
    open(unit=20, file='simulation_results_method.data', status='replace', action='write')
    do k = 1, N+1
        write(20, '(6(1X, E22.5, 1X))', advance='yes') t_arr(k), S(k), I(k), Q(k), & 
                                                    R(k), D(k)
    end do
    close(20)

    

    !------- Dont forget to deallocate arrays-----------!
    deallocate(t_arr, S, I, Q, R, D)

end program callSolver