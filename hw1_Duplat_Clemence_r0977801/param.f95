module Parametres
    implicit none
    
    !-------------Module with my parametsr so I don't need to write them systematically over-----------!
    integer, parameter :: KREAL = selected_real_kind(6, 37)
    real(KREAL) :: beta, mu, gamma, alpha, delta, S0, I0, &
                   s_fk,i_fk,q_fk,r_fk,d_fk,s_2,i_2,q_2,  &
                   r_2,d_2,s_heun,i_heun,q_heun,r_heun,d_heun
    real(KREAL) :: dt
    real(KREAL) :: T 
    integer :: N 
    
end module Parametres
