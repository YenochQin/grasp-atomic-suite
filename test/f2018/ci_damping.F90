! Compile the unmodified production MATRIX body; only assembly/eigensolution
! and memory-record input are controlled seams. BLAS is the linked real library.
module ci_fixture
    use vast_kind_param, only: double, long
    use eigv_C, only: evec, eval
    implicit none
    integer :: solves = 0
contains
    subroutine eigenpairs()
        solves = solves + 1
        eval = [-2.0_double, -1.0_double]
        if (solves == 1) then
            evec(3:6) = [1.0_double, 0.0_double, 0.0_double, 1.0_double]
        else
            evec(3:6) = [0.6_double, 0.8_double, -0.8_double, 0.6_double]
        endif
    end subroutine
end module

program ci_damping
    use ci_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use damp_C, only: cdamp
    use def_C, only: iccmin, ncmin, ncmax
    use debug_C, only: ldbpg
    use hblock_C, only: nblock, ncfblk, nevblk, ncmaxblk
    use hmat_C, only: emt, irow, iendc
    use mcpa_C, only: kmaxf
    use orb_C, only: ncf, nw
    use pos_C, only: ncftot, nvecsiz, ncfpast, ncminpast, nevecpast
    use peav_C, only: eavblk
    use syma_C, only: jpgg
#ifdef CI_MPI
    use mpi_C
    use matrixmpi_I
#else
    use mpi_s
    use matrix_I
#endif
#ifdef CI_MEM
    use rmcdhf_mem_C, only: imin1_30, imax1_30, imin2_30, imax2_30
#endif
    implicit none
    integer :: unit, i, records, local_size
    real(double), parameter :: tol = 1e-12_double
#ifdef CI_MPI
    call startmpi(myid,nprocs,host,lenhost)
#else
    myid = 0
    nprocs = 1
#endif
    nblock=1; ncf=2; ncftot=2; nw=1; ncmin=2; ncmax=2; nvecsiz=8; kmaxf=1
    ldbpg=.false.; jpgg=1
    allocate(ncfblk(1), nevblk(1), ncmaxblk(1), ncfpast(1), ncminpast(1), nevecpast(1), eavblk(1))
    allocate(iccmin(2), cdamp(2), evec(8), eval(2))
    ncfblk=2; nevblk=2; ncmaxblk=2; ncfpast=0; ncminpast=0; nevecpast=2
    iccmin=[1,2]; cdamp=[0.0_double,0.5_double]; evec=99.0_double
    nullify(emt,irow,iendc)
    local_size = count([(mod(i-1,nprocs)==myid, i=1,2)])
#ifdef CI_MEM
    allocate(imin1_30(1),imax1_30(1),imin2_30(1),imax2_30(1))
    imin1_30=0; imax1_30=2; imin2_30=0; imax2_30=local_size
#else
    do unit=30,33
        open(unit,status='scratch',form='unformatted')
        records=3
        if(unit==30) records=6
        do i=1,records
            write(unit) 0
        enddo
        if(unit==30) then
            write(unit) 'MCP',1,2
            write(unit) local_size
            write(unit) [(i,i=1,local_size)], [(i,i=1,local_size)]
        endif
    enddo
#endif
    if(myid==0) then
        open(25,status='scratch',form='unformatted')
        write(25) 'G92MIX'
        write(25) 2,2,1,2,8,1
    endif
    call matrix_step(.true.)
    ! With fcheck=all, the production allocation/zeroing must include IENDC(0).
    if(maxval(abs(evec(3:6)-[1.0_double,0.0_double,0.0_double,1.0_double]))>tol) &
        error stop 'FIRST must bypass damping'
    call matrix_step(.false.)
    if(.not.all(ieee_is_finite(evec))) error stop 'nonfinite damped roots'
    if(maxval(abs(evec(3:6)-[0.6_double,0.8_double,0.8_double,-0.6_double]))>tol) &
        error stop 'DDOT damping branch differs from analytic roots'
    if(abs(dot_product(evec(3:4),evec(5:6)))>tol) error stop 'orthogonality'
    if(abs(sum(evec(3:4)**2)-1.0_double)>tol.or.abs(sum(evec(5:6)**2)-1.0_double)>tol) &
        error stop 'normalization'
    if(any(evec(1:2)/=99.0_double).or.any(evec(7:8)/=99.0_double)) error stop 'block offset guards'
    cdamp=0.0_double
    call matrix_step(.false.)
    if(maxval(abs(evec(3:6)-[0.6_double,0.8_double,-0.8_double,0.6_double]))>tol) &
        error stop 'zero damping control'
    if(solves/=3) error stop 'fixture did not execute three production calls'
    if(myid==0) then
        close(25)
        print *, 'PASS: production FIRST / nonzero multi-root DDOT / zero damping'
    endif
#ifndef CI_MEM
    do unit=30,33
        close(unit)
    enddo
#endif
#ifdef CI_MPI
    call MPI_Finalize(ierr)
#endif
contains
    subroutine matrix_step(first)
        logical,intent(in)::first
#ifdef CI_MPI
        call matrixmpi(first)
#else
        call matrix(first)
#endif
    end subroutine
end program

#ifdef CI_MPI
subroutine maneigmpi(dvdfirst,lprint,jblock,ncfpat,ncminpat,nevecpat,ncftot)
#else
subroutine maneig(dvdfirst,lprint,jblock,ncfpat,ncminpat,nevecpat,ncftot)
#endif
    use ci_fixture
    implicit none
    logical,intent(in)::dvdfirst
    logical::lprint
    integer,intent(in)::jblock,ncfpat,ncminpat,nevecpat
    integer::ncftot
    if(jblock/=1.or.nevecpat/=2.or.ncftot/=2) error stop 'eigensolver seam arguments'
    call eigenpairs()
end subroutine

#ifdef CI_MEM
subroutine setham(jblock,myid,nprocs,ncontr_tot,ncontr_tot2)
#else
subroutine setham(jblock,myid,nprocs)
#endif
    use vast_kind_param, only: double,long
    use hmat_C, only: emt
    implicit none
    integer,intent(in)::jblock,myid,nprocs
#ifdef CI_MEM
    integer(long)::ncontr_tot(20),ncontr_tot2(20)
#endif
    if(jblock/=1) error stop 'Hamiltonian seam block'
    emt=0.0_double
end subroutine

#ifdef CI_MEM
subroutine read03_mem(jblock)
    use hmat_C, only: irow,iendc
#ifdef CI_MPI
    use mpi_C, only: myid,nprocs
#else
    use mpi_s, only: myid,nprocs
#endif
    implicit none
    integer,intent(in)::jblock
    integer::i,k
    if(jblock/=1) error stop 'memory seam block'
    k=0
    do i=myid+1,2,nprocs
        k=k+1
        iendc(i)=k
        irow(k)=i
    enddo
end subroutine
subroutine read1_mem(nfile,ncontr_tot2,lab,ncontr)
    use vast_kind_param, only: long
    implicit none
    integer,intent(in)::nfile
    integer(long),intent(in)::ncontr_tot2
    integer,intent(out)::lab,ncontr
    error stop 'unexpected empty block seam'
end subroutine
#endif
