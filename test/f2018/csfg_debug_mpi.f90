program csfg_debug_mpi
    use mpi_C
    use csfg_scf_C, only: ldbpdacon
    use setdbgmpi_I, only: setdbgmpi
    implicit none
    call startmpi(myid,nprocs,host,lenhost)
    ldbpdacon=[.false.,.true.]
    call setdbgmpi('coefficients')
    if(.not.ldbpdacon(1).or.ldbpdacon(2)) error stop 'DACON coefficient flag broadcast'
    call setdbgmpi('radial')
    if(ldbpdacon(1).or..not.ldbpdacon(2)) error stop 'DACON radial flag broadcast'
    if(myid==0) print *, 'PASS CSFG diagnostic broadcasts'
    call MPI_Finalize(ierr)
end program

! Controlled root-only input; the production MPI distributor is tested.
subroutine setdbg(filename)
    use debug_C, only: ldbpa,ldbpg,ldbpr
    use csfg_scf_C, only: ldbpdacon
    implicit none
    character(*),intent(in) :: filename
    ldbpa=.false.;ldbpg=.false.;ldbpr=.false.
    select case(filename)
    case('coefficients')
        ldbpdacon=[.true.,.false.]
    case('radial')
        ldbpdacon=[.false.,.true.]
    case default
        error stop 'Unexpected diagnostic input'
    end select
end subroutine
