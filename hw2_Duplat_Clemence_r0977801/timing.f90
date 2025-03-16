module timing
    use Parametres
    implicit none
    real(KREAL) :: start_time

contains
    subroutine tic()
        call cpu_time(start_time)
    end subroutine tic

    function toc() result(elapsed_time)
        real(KREAL) :: elapsed_time, end_time
        call cpu_time(end_time)
        elapsed_time = end_time - start_time
    end function toc
end module timing
