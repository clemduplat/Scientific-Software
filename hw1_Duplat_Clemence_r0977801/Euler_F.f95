program Euler_F
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
    open(unit=10, file='parameters.in', status='old', action='read')
    read(10, *) beta, mu, gamma, alpha, delta, S0, I0
    close(10)
    !T=10.0
    !N=100
    !beta=10.00
    !gamma=10.00
    !alpha=1.00
    !delta=0.00
    !mu=0.00
    !S0=100.00
    !I0=5.00
    !T=30.0
    !N=150
    T=10.0
    N=100
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
        
        !!-------------Call forward Euler Method-----------!
        s_fk=S_der(S(k),I(k),R(k))
        call forward_euler(dt,S(k),s_fk,S(k+1)) 
        i_fk=I_der(S(k),I(k),R(k)) 
        call forward_euler(dt,I(k),i_fk,I(k+1)) 
        q_fk=Q_der(I(k),Q(k))
        call forward_euler(dt,Q(k),q_fk,Q(k+1)) 
        r_fk=R_der(I(k),Q(k),R(k))
        call forward_euler(dt,R(k),r_fk,R(k+1)) 
        d_fk=D_der(I(k),Q(k))
        call forward_euler(dt,D(k),d_fk,D(k+1)) 
 
        
    end do
    somme=S(k)+I(k)+Q(k)+R(k)+D(k)
    print*,"somme 105=",somme
    !----------------Output the results in a file so that I have : t S I Q R D at each iteration---------!
    open(unit=20, file='simulation_results_forwardeuler.txt', status='replace', action='write')
    do k = 1, N+1
        write(20, '(6(1X, E12.5, 1X))', advance='yes') t_arr(k), S(k), I(k), Q(k), R(k), D(k)
        end do
    close(20)

    !------- Dont forget to deallocate arrays-----------!
    deallocate(t_arr, S, I, Q, R, D)

contains
    subroutine forward_euler(dt, x_0,f, x)
        real(KREAL), intent(in) :: dt,x_0,f
        real(KREAL), intent(out) :: x
        
        !----------Forward euler formula to compute approximation-------------!
        x = x_0 + dt * f
        end subroutine forward_euler
end program Euler_F