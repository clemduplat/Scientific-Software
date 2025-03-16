module Parametres
    implicit none
    
    !-------------Module with my parametsr so I don't need to write them systematically over-----------!
    !integer, parameter :: KREAL = selected_real_kind(6, 37)
    integer, parameter :: KREAL = selected_real_kind(15, 307)
    real(KREAL) :: beta, mu, gamma, alpha, delta, S0, I0, &
                   s_fk,i_fk,q_fk,r_fk,d_fk,s_2,i_2,q_2,  &
                   r_2,d_2,s_heun,i_heun,q_heun,r_heun,d_heun
    !---------------------------------!
    real(KREAL) :: dt
    integer, parameter :: max_iterations = 100
    real(KREAL), parameter :: tolerance = 1.0E-6_KREAL
    real(KREAL) :: T 
    integer :: N 
    !---------------!
    integer :: k
    real(KREAL), allocatable :: t_arr(:)
    real(KREAL), allocatable :: S(:), I(:), Q(:), R(:), D(:)
    integer :: method_number

    real(KREAL), dimension(5) :: x_k, x,f

    !-----Added for Tweaking module---!
    real(KREAL) :: deltaBeta, target
    real(KREAL) :: optimalBeta, currentMaxInfection
    real(KREAL) :: F_beta, F_prime_beta, F_beta_plus_delta
    integer :: iteration
    real(KREAL) :: lastBeta, parameterStep
    real(KREAL) ::dF
    
    !--------Added for Stability-------!
    real(KREAL) :: R0
    real(KREAL), dimension(5,5):: Jacobian
    real(KREAL), dimension(5) :: wr, wi  ! Real and imaginary parts of eigenvalues
    real(KREAL), dimension(:), allocatable :: work  ! Workspace
    real(KREAL) :: opti_workspace(1)
    integer :: lwork,z,size_jacobian
    integer :: lda, ldvl, ldvr, lwmax
    real(KREAL), dimension(5,5) :: v  ! Dummy array for VL and VR
contains
    subroutine get_parameters_file
        open(unit=10, file='parameters.in', status='old', action='read')
        read(10, *) beta, mu, gamma, alpha, delta, S0, I0
        close(10)
        !--------- N ant T are now read as first and second argument--------!
        write(*,*) "Give N (integer):"
        read(*,*) N
        write(*,*) "Give the T (real):"
        read(*,*) T
    end subroutine get_parameters_file
    subroutine get_parameters_tweaking
        beta= 0.4_KREAL
        gamma=0.3_KREAL
        alpha=0.005_KREAL
        delta=0.00_KREAL
        mu=0.00_KREAL
        !S0=0.00
        !I0=5.00
        T=100.0_KREAL
        N=1000
        
        write(*,*) "Give a target"
        read(*,*) target
        !target= 20000
        
    end subroutine get_parameters_tweaking
    subroutine get_parameters_tweaking_into_table
        beta= 0.4_KREAL
        gamma=0.3_KREAL
        alpha=0.005_KREAL
        delta=0.00_KREAL
        mu=0.00_KREAL
        !S0=0.00
        !I0=5.00
        T=100.0_KREAL
        N=1000
        
        
        target=20000
       
    end subroutine get_parameters_tweaking_into_table
    subroutine get_parameters_stability
        size_jacobian=5
        lda = size_jacobian
        ldvl = 1  ! Not computing left eigenvectors
        ldvr = 1  ! Not computing right eigenvectors
        lwmax = 1000  ! Max limit for workspace size
        !beta= 0.4_KREAL
        !gamma=0.3_KREAL
        !alpha=0.005_KREAL
        !delta=0.00_KREAL
        !mu=0.00_KREAL
        !S0=0.00
        !I0=5.00
        T=100.0_KREAL
        N=1000_KREAL
        !------Compute Jacobian around equilibium point---!
        !S0= 1.0_KREAL
        !I0=0.0_KREAL
        open(unit=10, file='eigenvalues2.in', status='old', action='read')
        read(10, *) beta, mu, gamma, alpha, delta, S0, I0
        close(10)
        R0=0.0_KREAL
    end subroutine get_parameters_stability
    subroutine get_parameters_open
        beta=10.00
        gamma=10.00
        alpha=1.00
        delta=0.00
        mu=0.00
        S0=100.00
        I0=5.00
        !--------- N ant T are now read as first and second argument--------!
        write(*,*) "Give N (integer):"
        read(*,*) N
        write(*,*) "Give the T (real):"
        read(*,*) T
    end subroutine get_parameters_open
    subroutine get_method
        write(*,*) "Wich method do you want : Forward Euler(1), Heun (2), Backward Euler(3)"
        read(*,*) method_number

    end subroutine get_method
    
end module Parametres
