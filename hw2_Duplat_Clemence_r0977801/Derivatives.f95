module Derivatives
    use Parametres
    implicit none
    contains
    subroutine calculate_derivatives(S, I, Q, R, D, derivatives)
        real(KREAL), intent(in) :: S, I, Q, R, D
        real(KREAL), dimension(5), intent(out) :: derivatives

        derivatives(1) = S_der(S, I, R)
        derivatives(2) = I_der(S, I, R)
        derivatives(3) = Q_der(I, Q)
        derivatives(4) = R_der(I, Q, R)
        derivatives(5) = D_der(I, Q)
    end subroutine calculate_derivatives

    function S_der(S,I,R) result(result_value)
        !-----Function to calculate S(t)-------!
        
        real(KREAL), intent(in):: I,S,R
        real(KREAL) :: result_value
        result_value= - beta * (I / (S + I + R)) * S + mu * R
        end function S_der

    function I_der(S,I,R) result(result_value)
        !-----Function to calculate I(t)-------!
        
        real(KREAL), intent(in):: I,S,R
        real(KREAL) :: result_value
        result_value=(beta * (S / (S + I + R)) - gamma - delta - alpha) * I
        end function I_der
    function Q_der(I,Q) result(result_value)
        !-----Function to calculate Q(t)-------!
        
        real(KREAL), intent(in):: I,Q
        real(KREAL):: result_value
        result_value=delta * I - (gamma + alpha) * Q
        end function Q_der
    function R_der(I,Q,R) result(result_value)
        !-----Function to calculate R(t)-------!
        
        real(KREAL), intent(in):: I,Q,R
        real(KREAL) :: result_value
        result_value= gamma * (I + Q) - mu * R
        end function R_der
    function D_der(I,Q) result(result_value)
        !-----Function to calculate D(t)-------!
        
        real(KREAL), intent(in):: I,Q
        real(KREAL):: result_value
        result_value=alpha * (I + Q)
        end function D_der
    function calculate_jacobian(S, I, R) result(Jacobian)
        !-------Function to calculate the Jacobian of my SIQRD model---------!
        
        real(KREAL):: Jacobian(5, 5)
        real(KREAL), intent(in) :: S, I, R
        
        real(KREAL) dSdS, dSdI, dSdR, dSdQ, dSdD
        real(KREAL) dIdS, dIdI, dIdR, dIdQ, dIdD
        real(KREAL) dRdS, dRdI, dRdR, dRdQ, dRdD
        real(KREAL) dQdS, dQdI, dQdR, dQdQ, dQdD
        real(KREAL) dDdS, dDdI, dDdR, dDdQ, dDdD

        !----------DS-------!
        !dSdS = -beta * (I*(I+R)/ ((S + I + R)**2) )
        dSdS = (beta * I /(S+I+R))*((S/(S+I+R)) -1 )
        !dSdI = -beta * S *(R+S)/ ((S + I + R)**2)
        dSdI = (beta * S /(S+I+R) )*((I/(S+I+R)) -1) 
        dSdQ = 0.0_KREAL 
        dSdR = beta * S *I/((S + I + R)**2) + mu
        dSdD = 0.0_KREAL

        !----------DI-------!
        !dIdS = beta * I*(I+R) / ((S + I + R)**2)
        dIdS =(beta * I/(S+I+R)) * (1-S/(S+I+R)) 
        !dIdI = - beta * S *I/ ((S + I + R)**2) + beta*S/(I+R+S) - ( gamma + delta + alpha)
        dIdI= (beta*S /(S+I+R))*(1-I/(S+I+R)) - gamma - delta - alpha
        dIdQ = 0.0_KREAL
        dIdR = - (beta * I* S) / ((S + I + R)**2) 
        dIdD = 0.0_KREAL

        !----------DQ-------!
        dQdS = 0.0_KREAL
        dQdI = delta
        dQdQ = -(gamma + alpha)
        dQdR = 0.0_KREAL
        dQdD = 0.0_KREAL

        !----------DR-------!
        dRdS = 0.0_KREAL
        dRdI = gamma
        dRdQ = gamma
        dRdR = - mu
        dRdD = 0.0_KREAL

        !----------DD-------!
        dDdS = 0.0_KREAL
        dDdI = alpha
        dDdQ = alpha
        dDdR = 0.0_KREAL
        dDdD = 0.0_KREAL

        Jacobian(1, 1) = dSdS
        Jacobian(1, 2) = dSdI
        Jacobian(1, 3) = dSdQ
        Jacobian(1, 4) = dSdR
        Jacobian(1, 5) = dSdD

        Jacobian(2, 1) = dIdS
        Jacobian(2, 2) = dIdI
        Jacobian(2, 3) = dIdQ
        Jacobian(2, 4) = dIdR
        Jacobian(2, 5) = dIdD

        Jacobian(3, 1) = dQdS
        Jacobian(3, 2) = dQdI
        Jacobian(3, 3) = dQdQ
        Jacobian(3, 4) = dQdR
        Jacobian(3, 5) = dQdD

        Jacobian(4, 1) = dRdS
        Jacobian(4, 2) = dRdI
        Jacobian(4, 3) = dRdQ
        Jacobian(4, 4) = dRdR
        Jacobian(4, 5) = dRdD

        Jacobian(5, 1) = dDdS
        Jacobian(5, 2) = dDdI
        Jacobian(5, 3) = dDdQ
        Jacobian(5, 4) = dDdR
        Jacobian(5, 5) = dDdD
        end function calculate_jacobian
    function identity_matrix(n) result(mat)
        !--------Function to calculate an identity matrix of size n------!
        integer, intent(in) :: n
        real(KREAL), dimension(n, n) :: mat
        integer :: i
        mat = 0.0_KREAL
        do i = 1, n
            mat(i, i) = 1.0_KREAL
        end do
        end function identity_matrix

end module Derivatives