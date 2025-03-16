program back_euler 
    use solver_gfortran
    use Derivatives
    
    
    
    implicit none

    
    integer :: k
    
    real(KREAL), allocatable :: t_arr(:)
    real(KREAL), allocatable :: S(:), I(:), Q(:), R(:), D(:)
    real(KREAL), dimension(5) :: x_k, x,f
    integer, parameter :: max_iterations = 1000
    real(KREAL), parameter :: tolerance = 1.0E-6_KREAL
    real(KREAL) :: somme
    ! ------Open file for parameters or put them yourself--------!
    !open(unit=10, file='parameters.in', status='old', action='read')
    !read(10, *) beta, mu, gamma, alpha, delta, S0, I0
    !close(10)
    !T=10.0
    !N=100
    beta=10.00
    gamma=10.00
    alpha=1.00
    delta=0.00
    mu=0.00
    S0=100.00
    I0=5.00
    !T=30.0
    !N=150
    T=30.0
    N=150
    !-------------Allocate arrays and put initial conditions-----------!
    dt = T / real(N, KREAL)
    allocate(t_arr(1:N+1), S(1:N+1), I(1:N+1), Q(1:N+1), R(1:N+1), D(1:N+1))

    t_arr(1) = 0.0
    S(1) = S0
    I(1) = I0
    Q(1) = 0.0
    R(1) = 0.0
    D(1) = 0.0
    
    do k = 1, N
        t_arr(k+1) = t_arr(k) + dt
        
        s_fk=S_der(S(k),I(k),R(k))
    
        i_fk=I_der(S(k),I(k),R(k)) 
        
        q_fk=Q_der(I(k),Q(k))
        
        r_fk=R_der(I(k),Q(k),R(k))
        
        d_fk=D_der(I(k),Q(k))
        
        !!-------------Call Backward Euler Method---------------!
        !---- Put SIQRD in a vector so that I can use Jacobian and identity as asked in the assignment-----!
        f=[s_fk,i_fk,q_fk,r_fk,d_fk]
        x_k = [S(k), I(k), Q(k), R(k), D(k)]
        call backward_euler(dt, x_k, f, x, max_iterations, tolerance) ! remove the ! if you wan't to use backward
        S(k+1) = x(1) ! remove the ! if you wan't to use backward
        I(k+1) = x(2) ! remove the ! if you wan't to use backward
        Q(k+1) = x(3) ! remove the ! if you wan't to use backward
        R(k+1) = x(4) ! remove the ! if you wan't to use backward
        D(k+1) = x(5) ! remove the ! if you wan't to use backward
        
        
    end do
    somme=S(k)+I(k)+Q(k)+R(k)+D(k)
    !print*,"La somme doit ête égale à 105",somme
    !----------------Output the results in a file so that I have : t S I Q R D at each iteration---------!
    open(unit=20, file='simulation_results_backwardeuler.data', status='replace', action='write')
    do k = 1, N+1
        write(20, '(6(1X, E12.5, 1X))', advance='yes') t_arr(k), S(k), I(k), Q(k), R(k), D(k)
        end do
    close(20)

    !------- Don't forget to deallocate arrays-----------!
    deallocate(t_arr, S, I, Q, R, D)

contains
    subroutine backward_euler(dt, x_0, f, x, max_iterations, tolerance)
        !A popular method to solve such systems of non-linear equations
        !is Newton’s method. Starting from an initial guess x(0)k+1, this method katively updates its estimate
        real(KREAL), intent(in) :: dt, tolerance
        integer, intent(in) :: max_iterations
        real(KREAL), dimension(5),intent(in) :: x_0,f
        real(KREAL), dimension(5),intent(out) :: x
        
         

        real(KREAL), dimension(5, 5) :: A
        real(KREAL), dimension(5) :: b
        real(KREAL), dimension(5) :: x_new

        integer :: k
        integer :: nb_iteration_forward_error = 0
        integer :: nb_iteration_backward_error = 0

        real(KREAL) :: forward_error
        real(KREAL) :: backward_error
        !--------------Starting from an initial guess-------------!
        x = x_0 
        ! ------------------iterate untill one of the stopping criteria is met ---------------------!
        do k = 1, max_iterations
            ! dt* \frac{\delta f}{\delta x}(x^{s}_k+1)-I_{5}
            A = calculate_jacobian(x_0(1), x_0(2), x_0(4)) * dt - identity_matrix(5)

            ! (x_k +dt*f(x^{s}_{k+1})-x^{s}_{k+1}) in  equation (2) from assignment
            b = x_0 + dt * f - x

            ! Function solve is given to us, make us able to "calculate the inverse" that will be stocked in b as bx
            call solve(A, b)

            ! x^{s+1}_k+1=x^{s}_k+1 -b -----> our initial equation
            x_new = x_0 - b

            !------ First stop criterium: Calculate if the forward error is small enough ---------------!
            forward_error = norm2(x_new - x)
            !print*, "forward error is: ", forward_error
            !if (forward_error < tolerance) then
                !print*, "Number of iterations with forward error:", nb_iteration_forward_error
                !print*,"Condition be met"
                !exit
            !end if
            nb_iteration_forward_error = nb_iteration_forward_error + 1

            !------ Second stop criterium: Calculate the backward error ---------------!
            backward_error = norm2(x_0 + dt*f - x_new)
            !print*, "backward error is: ", backward_error

            if (backward_error < tolerance) then
                !print*, "Number of iterations with backward error:", nb_iteration_backward_error
                !print*,"Condition be met"
                exit
            end if
            nb_iteration_backward_error = nb_iteration_backward_error + 1

            ! ----------------Update x for the next iteration----------------!
            x = x_new
            
        end do

        
    end subroutine backward_euler
end program back_euler