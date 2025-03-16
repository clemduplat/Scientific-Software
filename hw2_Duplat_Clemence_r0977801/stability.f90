program stability
    use Derivatives
    implicit none
    

    call get_parameters_stability
    Jacobian= calculate_jacobian(S0, I0, R0)
    !----dgeev lapack librarie----!
    !---Recommend by lapack documentation---!
    !---First call to dgeev-----!
    !If LWORK = -1, then a workspace query is assumed; the routine
    !      only calculates the optimal size of the WORK array, returns
    !      this value as the first entry of the WORK array, and no error
    !      message related to LWORK is issued by XERBLA.
    ! In other words we will call dgeev with lwork=-1 to find the optimal workspace size required to compute its eigenvalues 
    lwork = -1
    call dgeev('N', 'N', size_jacobian, Jacobian, lda, wr, wi, v, ldvl, v, ldvr, opti_workspace, lwork, info)
    if (info /= 0) then
        print*, "Error"
        stop
    end if
    ! round the number to the nearest integer to be sure and have enough size for the workspace
    lwork = min(lwmax, ceiling(opti_workspace(1)))
    allocate(work(lwork))

    !------Second call to DGEEV-----!
    ! this time to compute eigenvalues, with correct workspace !
    call dgeev('N', 'N', size_jacobian, Jacobian, lda, wr, wi, v, ldvl, v, ldvr, work, lwork, info)
    if (info /= 0) then
        print*, "Error"
        stop
    end if

    
    !----Write real part and imaginary part----!
    print *, 'Eigenvalues of the Jacobian:'
    do z = 1, size_jacobian
        print *, 'Eigenvalue ', z, ': ', wr(z), ' + z*', wi(z)
        write(10, *) wr(z), wi(z)
    end do
  

    !----- don't forget to deallocate-----!
    if (allocated(work)) then
        deallocate(work)
    end if
   
end program stability