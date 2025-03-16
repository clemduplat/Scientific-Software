module solvers
    use Derivatives
    use solver
    !use solver_gfortran
    
    implicit none
    
    contains
    subroutine solver_forward(S,I,Q,R,D,S_1,I_1,Q_1,R_1,D_1)
        real(KREAL), intent(in) :: S,I,Q,R,D
        real(KREAL), intent(out) :: S_1,I_1,Q_1,R_1,D_1
        real(KREAL), dimension(5) :: fk

        call calculate_derivatives(S, I, Q, R, D, fk)

        call forward_euler(dt, S, fk(1), S_1)
        call forward_euler(dt, I, fk(2), I_1)
        call forward_euler(dt, Q, fk(3), Q_1)
        call forward_euler(dt, R, fk(4), R_1)
        call forward_euler(dt, D, fk(5), D_1)
        end subroutine solver_forward
    subroutine forward_euler(dt, x_0,f, x)

        real(KREAL), intent(in) :: dt,x_0,f
        real(KREAL), intent(out) :: x
        
        !----------Forward euler formula to compute approximation-------------!
        x = x_0 + dt * f
        end subroutine forward_euler
    subroutine solver_heun(S, I, Q, R, D, S_1, I_1, Q_1, R_1, D_1)
        real(KREAL), intent(in) :: S, I, Q, R, D
        real(KREAL), intent(out) :: S_1, I_1, Q_1, R_1, D_1
        real(KREAL), dimension(5) :: f_fk, f_heun
        real(KREAL) :: s_2, i_2, q_2, r_2, d_2
        call calculate_derivatives(S, I, Q, R, D, f_fk)
        ! Predict intermediate values
        s_2 = S + dt * f_fk(1)
        i_2 = I + dt * f_fk(2)
        q_2 = Q + dt * f_fk(3)
        r_2 = R + dt * f_fk(4)
        d_2 = D + dt * f_fk(5)
        call calculate_derivatives(s_2, i_2, q_2, r_2, d_2, f_heun)
        call heun(dt, S, f_fk(1), f_heun(1), S_1)
        call heun(dt, I, f_fk(2), f_heun(2), I_1)
        call heun(dt, Q, f_fk(3), f_heun(3), Q_1)
        call heun(dt, R, f_fk(4), f_heun(4), R_1)
        call heun(dt, D, f_fk(5), f_heun(5), D_1)
    end subroutine solver_heun


    subroutine heun(dt,x_0,f,f_k,x)
        real(KREAL), intent(in) :: dt,x_0,f,f_k
        real(KREAL), intent(out) :: x
        !--------------Heun formula to compute approximation------------------!
        x= x_0 + dt*(0.5_KREAL*f+0.5_KREAL*f_k)
        end subroutine heun
    subroutine solver_backward(S, I, Q, R, D, S_1, I_1, Q_1, R_1, D_1)
        real(KREAL), intent(in) :: S, I, Q, R, D
        real(KREAL), intent(out) :: S_1, I_1, Q_1, R_1, D_1
        real(KREAL), dimension(5) :: x, x_k, f
        real(KREAL) :: total_population
        call calculate_derivatives(S, I, Q, R, D, f)

        x_k = [S, I, Q, R, D]
        call backward_euler(dt, x_k, f, x) 
        total_population=11e6
        S_1=x(1)
        I_1=x(2)
        Q_1=x(3)
        R_1=x(4)
        D_1=x(5)
        S_1 = max(min(x(1), total_population), 0.0_KREAL)
        I_1 = max(min(x(2), total_population), 0.0_KREAL)
        Q_1 = max(min(x(3), total_population), 0.0_KREAL)
        R_1 = max(min(x(4), total_population), 0.0_KREAL)
        D_1 = max(min(x(5), total_population), 0.0_KREAL)
    end subroutine solver_backward

    subroutine backward_euler(dt, x_0, f, x)
        !A popular method to solve such systems of non-linear equations
        !is Newton’s method. Starting from an initial guess x(0)k+1, this method katively updates its estimate
        real(KREAL),intent(in) :: dt
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
        real(KREAL) :: absolute_change, custom_tolerance

        real(KREAL) :: reduced_dt
        logical :: nan_detected
        custom_tolerance = 1.0E-6_KREAL  
        reduced_dt = dt
        nan_detected = .false.
        !--------------Starting from an initial guess-------------!
        x = x_0 
        ! ------------------iterate untill one of the stopping criteria is met ---------------------!
        do k = 1, max_iterations
            ! dt* \frac{\delta f}{\delta x}(x^{s}_k+1)-I_{5}
            
            A = calculate_jacobian(x_0(1), x_0(2), x_0(4)) * dt - identity_matrix(5)

            ! (x_k +dt*f(x^{s}_{k+1})-x^{s}_{k+1}) in  equation (2) from assignment
            b = x_0 + dt * f - x
            if(norm2(b)<tolerance*norm2(x)) then
                !print*,"C1"
                exit 
            endif
            ! Function solve is given to us, make us able to "calculate the inverse" that will be stocked in b as bx
            call solve(A, b)

            ! x^{s+1}_k+1=x^{s}_k+1 -b -----> our initial equation
            x_new = x_0 - b
            nan_detected = any(isnan(x_new))

            if (nan_detected) then
                !------Handle NaN, last known good state and reduce step size-----!
                x = x_0
                reduced_dt = reduced_dt / 2  
                cycle  !Skip rest of the loop and start the next iteration!
            endif
            !------ First stop criterium: Calculate if the forward error is small enough ---------------!
            forward_error = norm2(x_new - x)
            !print*, "forward error is: ", forward_error
            if (forward_error < custom_tolerance) then
                !print*, "Number of iterations with forward error:", nb_iteration_forward_error
                !print*,"Condition be met"
                exit
            end if
            nb_iteration_forward_error = nb_iteration_forward_error + 1

            !------ Second stop criterium: Calculate the backward error ---------------!
            backward_error = norm2(x_0 + dt*f - x_new)
            !print*, "backward error is: ", backward_error

            !if (backward_error < tolerance) then
                !print*, "Number of iterations with backward error:", nb_iteration_backward_error
                !print*,"Condition be met"
                !exit
            !end if
            nb_iteration_backward_error = nb_iteration_backward_error + 1
            !---------- Implement a more suitable stopping criteria----!
            
            
            
            ! ----------------Update x for the next iteration----------------!
            x = x_new
            
        end do

        
        end subroutine backward_euler

    end module solvers