program exponential_growth
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307) ! Double precision given format in the programs session 1
    integer, parameter :: sp = selected_real_kind(6, 37)   ! Single precision given format in the programs session 1
    integer, parameter :: precision = dp ! to make it easy to switch between the 2 precision format
    ! Change the line 5 of code to switch: precision=sp (single precision) or precision=dp (double precsion)
    real(kind=precision) :: Id,r,d,Io
    
    write(*,*) "Give the number of initial infections :"
    read(*,*) Io
    write(*,*) "Give the number of days : "
    read(*,*) d
    write(*,*) "Give the daily infection rate: "
    read(*,*) r

    Id=((1.0_precision+r)**d)*Io
    write(*,*) "The number of infected individauls after a certain number of days is: ",Id
    

    end program exponential_growth