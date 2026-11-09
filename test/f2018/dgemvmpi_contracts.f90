module error_fixture
    implicit none
    integer :: last_info=0
end module

program dgemvmpi_contracts
    use vast_kind_param, only: double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mpi_C
    use dvd_mpi_I, only: dgemvmpi
    use error_fixture
    implicit none
    real(double) :: a(3,2),x(3),guarded(5)
    integer :: kind
    character(1) :: trans
    call startmpi(myid,nprocs,host,lenhost)
    a=reshape([1.0_double,2.0_double,3.0_double,4.0_double,5.0_double,6.0_double],[3,2])
    x=[2.0_double,-1.0_double,0.5_double]
    do kind=1,3
        trans='N'
        if(kind==2) trans='T'
        if(kind==3) trans='C'
        guarded=99
        call dgemvmpi(trans,3,2,1.0_double,a,3,x,1,0.0_double,guarded(2),1)
        if(last_info/=0.or..not.all(ieee_is_finite(guarded))) error stop 'DGEMVMPI valid call'
        if(kind==1) then
            if(maxval(abs(guarded(2:4)-[-2.0_double,-1.0_double,0.0_double]))>0) error stop 'DGEMVMPI N'
        else
            if(maxval(abs(guarded(2:4)-[1.5_double,6.0_double,99.0_double]))>0) error stop 'DGEMVMPI T/C'
        end if
        if(abs(guarded(1)-99.0_double)+abs(guarded(5)-99.0_double)>0) error stop 'DGEMVMPI guards'
    end do
    guarded=99
    call dgemvmpi('N',3,2,1.0_double,a,3,x,2,0.0_double,guarded(2),1)
    if(last_info/=8.or.maxval(abs(guarded-99.0_double))>0) error stop 'INCX must be rejected before access'
    last_info=0
    call dgemvmpi('T',3,2,1.0_double,a,3,x,1,0.0_double,guarded(2),-1)
    if(last_info/=11.or.maxval(abs(guarded-99.0_double))>0) error stop 'INCY must be rejected before access'
    last_info=0
    call dgemvmpi('N',0,2,1.0_double,a,3,x,1,0.0_double,guarded(2),1)
    if(last_info/=0.or.maxval(abs(guarded-99.0_double))>0) error stop 'Empty product quick return'
    if(myid==0) print *, 'PASS MPI products, empty rank, guards and unsupported stride rejection'
    call MPI_Finalize(ierr)
end program

! Return from the error handler so the test proves the production routine
! itself returns before indexing unsupported strides or entering collectives.
subroutine xerbla(name,info)
    use error_fixture
    implicit none
    character(*),intent(in) :: name
    integer,intent(in) :: info
    if(name/='DGEMV ') error stop 'Unexpected BLAS error'
    last_info=info
end subroutine
