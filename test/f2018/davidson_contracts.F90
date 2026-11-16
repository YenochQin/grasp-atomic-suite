! Production initial estimates, sparse multiplication and GDVD callback contract.
program davidson_contracts
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64
    use dinit_I, only: dinit
    use iniest_I, only: iniest
    use iniest2_I, only: iniest2
    use gdvd_I, only: gdvd
    use spicmv2_I, only: spicmv2
    use adds_I, only: adds
    use newvec_I, only: newvec
    use hmat_c, only: emt => hamiltonian_values, irow => matrix_partner_indices, iendc => stored_segment_end_offsets
    implicit none (type, external)
    integer, parameter :: n=3, nd=6, lim=6, lwork=2*nd*lim+2*lim*lim+11*lim+1, liwork=6*lim+1
    integer :: niv, ilow, ihigh, mblock, maxiter, loops, nmv, ierr, i
    integer :: jcol(0:3), rows(3), blocks(3), selected(lim), iwork(liwork)
    real(real64) :: values(3), basis(8), work(lwork), diag(nd), h(n,n), big(nd,nd)
    real(real64) :: result(n,2), vectors(n,2), crite, critc, critr, ortho, guard(7)
    logical :: hiend
    character(16) :: scenario
    real(real64) :: expected
    guard = -99
    call dinit(3, 2.0_real64, guard(2:6), 2)
    if (any(guard /= [-99.0_real64,2.0_real64,-99.0_real64,2.0_real64, &
                     -99.0_real64,2.0_real64,-99.0_real64])) error stop 'DINIT stride/guards'
    call dinit(3, 4.0_real64, guard(2:6), -2)
    if (any(guard(2:6:2) /= 4.0_real64)) error stop 'DINIT negative stride'
    call dinit(0, 7.0_real64, guard, 1)
    if (guard(1) /= -99.0_real64) error stop 'DINIT zero count'

    values = [1.0_real64,2.0_real64,4.0_real64]
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

    h = reshape([2.0_real64,1.0_real64,0.0_real64, &
                 1.0_real64,2.0_real64,0.0_real64, &
                 0.0_real64,0.0_real64,4.0_real64], [n,n])
    allocate(emt(4),irow(4),iendc(0:n))
    emt = [2.0_real64,1.0_real64,2.0_real64,4.0_real64]
    irow = [1,1,2,3]
    iendc = [0,1,3,4]
    vectors = reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64,5.0_real64,6.0_real64],[n,2])
    result = -99
    call spicmv2(n,2,vectors,result)
    if (maxval(abs(result-matmul(h,vectors))) > 1e-12_real64) error stop 'SPICMV2 sparse result'
    deallocate(emt,irow,iendc)
    call check_partial_updates()

    ! GDVD reserves three restart slots, so use LIM=6 even for this small block.
    call get_command_argument(1, scenario)
    big=0
    select case (trim(scenario))
    case ('')
        do i=1,nd
            big(i,i)=2.0_real64
        end do
        do i=1,nd-1
            big(i,i+1)=-1.0_real64
            big(i+1,i)=-1.0_real64
        end do
        expected=2.0_real64-2.0_real64*cos(acos(-1.0_real64)/real(nd+1,real64))
    case ('early')
        ! An invariant two-dimensional block reaches an exact Ritz pair
        ! before a third expansion. It must never normalize a zero residual.
        do i=1,nd
            big(i,i)=real(i+1,real64)
        end do
        big(2,2)=2.0_real64
        big(1,2)=1.0_real64
        big(2,1)=1.0_real64
        expected=1.0_real64
    case default
        error stop 'Unknown Davidson scenario'
    end select
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
    crite=1e-14_real64
    critc=1e-12_real64
    critr=1e-12_real64
    ortho=1e-12_real64
    call gdvd(multiply,nd,lim,diag,ilow,ihigh,selected,niv,mblock,crite,critc,critr,ortho, &
              maxiter,work,lwork,iwork,liwork,hiend,loops,nmv,ierr)
    if (ierr /= 0 .or. nmv < 1) error stop 'GDVD callback execution'
    if (.not. all(ieee_is_finite(work(:nd+1)))) error stop 'GDVD nonfinite eigenpair'
    if (abs(work(nd+1)-expected) > 1e-12_real64) &
        error stop 'GDVD eigenvalue'
    if (maxval(abs(matmul(big,work(:nd))-work(nd+1)*work(:nd))) > 1e-12_real64) error stop 'GDVD residual'
    if (abs(sum(work(:nd)**2)-1.0_real64) > 1e-12_real64) error stop 'GDVD normalization'
    print *, 'PASS Davidson interfaces and numerical contracts'
contains
    subroutine check_partial_updates()
        real(real64) :: b(3,3), ab(3,3), packed(6), coeff(6), eigen(3), old(2), scratch(3), shifts(3)
        integer :: flags(2), pending(2), count_new
        logical :: done
        b=0
        b(1,1)=1; b(2,2)=1; b(3,3)=1
        packed=-99
        call adds(3,3,0,1,b,h,packed)
        if(abs(packed(1)-2.0_real64)>0.or.maxval(abs(packed(2:)+99.0_real64))>0) &
            error stop 'ADDS first-column guards'
        call adds(3,3,1,1,b,h,packed)
        call adds(3,3,2,1,b,h,packed)
        if(.not.all(ieee_is_finite(packed))) error stop 'Nonfinite projected matrix'
        if(maxval(abs(packed-[2.0_real64,1.0_real64,2.0_real64,0.0_real64,0.0_real64,4.0_real64]))>0) &
            error stop 'ADDS must preserve previous projected columns'
        ! Root 2 is already converged. Updating root 1 must retain that flag.
        ab=0; ab(1,1)=1; ab(2,2)=2
        coeff=[1.0_real64,0.0_real64,0.0_real64,1.0_real64,0.0_real64,0.0_real64]
        eigen=[1.0_real64,2.0_real64,3.0_real64];old=[0.0_real64,0.0_real64]
        flags=[0,1];pending=[1,2];count_new=1;scratch=-99;shifts=-99
        call newvec(3,2,3,1,2,1e-12_real64,count_new,pending,coeff,eigen,old,ab,b,flags,scratch,shifts,done)
        if(.not.done.or.any(flags/=1)) error stop 'NEWVEC must retain other roots convergence state'
    end subroutine
    subroutine check_estimates()
        if (any(abs(basis(7:8)-values(:2)) > 1e-12_real64)) error stop 'INIEST eigenvalues'
        if (abs(abs(basis(1))-1.0_real64) > 1e-12_real64) error stop 'INIEST first vector'
        if (abs(abs(basis(5))-1.0_real64) > 1e-12_real64) error stop 'INIEST second vector'
        if (any(abs(basis([2,3,4,6])) > 1e-12_real64)) error stop 'INIEST zero components'
    end subroutine
    subroutine multiply(size_n, count, b, c)
        integer, intent(in) :: size_n, count
        real(real64), intent(in) :: b(size_n,count)
        real(real64), intent(inout) :: c(size_n,count)
        c = matmul(big,b)
    end subroutine
end program
