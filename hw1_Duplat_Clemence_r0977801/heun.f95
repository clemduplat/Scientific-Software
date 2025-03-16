program heun_meth
    use Derivatives
    
    
    implicit none

    
    integer :: k
    
    real(KREAL), allocatable :: t_arr(:)
    real(KREAL), allocatable :: S(:), I(:), Q(:), R(:), D(:)
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
        
        !!-------------Call forward Euler Method-----------!
        s_fk=S_der(S(k),I(k),R(k))
    
        i_fk=I_der(S(k),I(k),R(k)) 
        
        q_fk=Q_der(I(k),Q(k))
        
        r_fk=R_der(I(k),Q(k),R(k))
        
        d_fk=D_der(I(k),Q(k))
        


        
        

        !!-------------Call Heun Method---------------!
        s_2= S(k) + dt*s_fk
        i_2= I(k) + dt*i_fk
        q_2= Q(k) + dt*q_fk
        r_2= R(k) + dt*r_fk
        d_2= D(k) + dt*d_fk

        s_heun=S_der(s_2,i_2,r_2)
        call heun(dt,S(k),s_fk,s_heun,S(k+1)) ! remove the ! if you wan't to use heun
        
        i_heun=I_der(s_2,i_2,r_2)
        call heun(dt,I(k),i_fk,i_heun,I(k+1)) ! remove the ! if you wan't to use heun
        
        q_heun=Q_der(i_2,q_2) 
        call heun(dt,Q(k),q_fk,q_heun,Q(k+1)) ! remove the ! if you wan't to use heun
        
        r_heun=R_der(i_2,q_2,r_2)
        call heun(dt,R(k),r_fk,r_heun,R(k+1)) ! remove the ! if you wan't to use heun
        
        d_heun=d_fk
        call heun(dt,D(k),d_fk,d_heun,D(k+1))  ! remove the ! if you wan't to use heun
            
        
    end do
    somme=S(k)+I(k)+Q(k)+R(k)+D(k)
    print*,"La somme doit ête égale à 105",somme
    !----------------Output the results in a file so that I have : t S I Q R D at each iteration---------!
    open(unit=20, file='simulation_results_heun.data', status='replace', action='write')
    do k = 1, N+1
        write(20, '(6(1X, E12.5, 1X))', advance='yes') t_arr(k), S(k), I(k), Q(k), R(k), D(k)
        end do
    close(20)

    !------- Don't forget to deallocate arrays-----------!
    deallocate(t_arr, S, I, Q, R, D)

contains
    subroutine heun(dt,x_0,f,f_k,x)
        real(KREAL), intent(in) :: dt,x_0,f,f_k
        real(KREAL), intent(out) :: x
        !--------------Heun formula to compute approximation------------------!
        x= x_0 + dt*(0.5_KREAL*f+0.5_KREAL*f_k)
        end subroutine heun
end program heun_meth