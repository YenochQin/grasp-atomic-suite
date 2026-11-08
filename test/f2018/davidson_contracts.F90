! Production initial estimates, sparse multiplication and GDVD callback contract.
program davidson_contracts
    use vast_kind_param, only: double
    use dinit_I, only: dinit
    use iniest_I, only: iniest
    use iniest2_I, only: iniest2
    use gdvd_I, only: gdvd
    use spicmv2_I, only: spicmv2
    use hmat_C, only: emt, irow, iendc
#ifdef DAVIDSON_MPI
    use mpi_C, only: myid, nprocs, host, lenhost, startmpi, MPI_Finalize
#endif
    implicit none
    integer, parameter :: n=3, nd=6, lim=6, lwork=2*nd*lim+2*lim*lim+11*lim+1, liwork=6*lim+1
    integer :: niv, ilow, ihigh, mblock, maxiter, loops, nmv, ierr, i
    integer :: jcol(0:3), rows(3), blocks(3), selected(lim), iwork(liwork)
    real(double) :: values(3), basis(8), work(lwork), diag(nd), h(n,n), big(nd,nd)
    real(double) :: result(n,2), vectors(n,2), crite, critc, critr, ortho, guard(7)
    logical :: hiend
#ifdef DAVIDSON_MPI
    call startmpi(myid,nprocs,host,lenhost)
#endif
    guard = -99
    call dinit(3, 2.0_double, guard(2:6), 2)
    if (any(guard /= [-99.0_double,2.0_double,-99.0_double,2.0_double, &
                     -99.0_double,2.0_double,-99.0_double])) error stop 'DINIT stride/guards'
    call dinit(3, 4.0_double, guard(2:6), -2)
    if (any(guard(2:6:2) /= 4.0_double)) error stop 'DINIT negative stride'
    call dinit(0, 7.0_double, guard, 1)
    if (guard(1) /= -99.0_double) error stop 'DINIT zero count'

    values = [1.0_double,2.0_double,4.0_double]
    rows = [1,2,3]
    jcol = [0,1,2,3]
    blocks = 1
    niv = 2
    basis = -99
    call iniest(n,n,niv,values,jcol,rows,basis,blocks,1)
    call check_estimates()
    basis = -99
    call iniest2(n,n,niv,basis,values,jcol,rows)
    call check_estimates()

    h = reshape([2.0_double,1.0_double,0.0_double, &
                 1.0_double,2.0_double,0.0_double, &
                 0.0_double,0.0_double,4.0_double], [n,n])
    allocate(emt(4),irow(4),iendc(0:n))
    emt = [2.0_double,1.0_double,2.0_double,4.0_double]
    irow = [1,1,2,3]
    iendc = [0,1,3,4]
    vectors = reshape([1.0_double,2.0_double,3.0_double,4.0_double,5.0_double,6.0_double],[n,2])
    result = -99
    call spicmv2(n,2,vectors,result)
    if (maxval(abs(result-matmul(h,vectors))) > 1e-12_double) error stop 'SPICMV2 sparse result'
    deallocate(emt,irow,iendc)

    ! GDVD reserves three restart slots, so use LIM=6 even for this small block.
    ! Use a fully coupled matrix: the MPI solver requires more than two
    ! iterations before accepting a residual (its existing stopping rule).
    big=0
    do i=1,nd
        big(i,i)=2.0_double
    end do
    do i=1,nd-1
        big(i,i+1)=-1.0_double
        big(i+1,i)=-1.0_double
    end do
    work=0
    work(1)=1
    work(nd+1)=big(1,1)
    do i=1,nd
        diag(i)=big(i,i)
    end do
    niv=1
    ilow=1
    ihigh=1
    mblock=1
    maxiter=30
    selected=0
    iwork=0
    loops=0
    hiend=.false.
    crite=1e-14_double
    critc=1e-12_double
    critr=1e-12_double
    ortho=1e-12_double
    call gdvd(multiply,nd,lim,diag,ilow,ihigh,selected,niv,mblock,crite,critc,critr,ortho, &
              maxiter,work,lwork,iwork,liwork,hiend,loops,nmv,ierr)
    if (ierr /= 0 .or. nmv < 1) error stop 'GDVD callback execution'
    if (abs(work(nd+1)-(2.0_double-2.0_double*cos(acos(-1.0_double)/real(nd+1,double)))) > 1e-12_double) &
        error stop 'GDVD eigenvalue'
    if (maxval(abs(matmul(big,work(:nd))-work(nd+1)*work(:nd))) > 1e-12_double) error stop 'GDVD residual'
    print *, 'PASS Davidson interfaces and numerical contracts'
#ifdef DAVIDSON_MPI
    call MPI_Finalize(ierr)
#endif
contains
    subroutine check_estimates()
        if (any(abs(basis(7:8)-values(:2)) > 1e-12_double)) error stop 'INIEST eigenvalues'
        if (abs(abs(basis(1))-1.0_double) > 1e-12_double) error stop 'INIEST first vector'
        if (abs(abs(basis(5))-1.0_double) > 1e-12_double) error stop 'INIEST second vector'
        if (any(abs(basis([2,3,4,6])) > 1e-12_double)) error stop 'INIEST zero components'
    end subroutine
    subroutine multiply(size_n, count, b, c)
        integer, intent(in) :: size_n, count
        real(double), intent(in) :: b(size_n,count)
        real(double), intent(inout) :: c(size_n,count)
        c = matmul(big,b)
    end subroutine
end program
