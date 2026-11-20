! Shares the Slater integral cache and counts each rank's share of integral evaluations.
module distributed_slater_integrals
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer, pointer :: indteirk(:)
    real(real64), pointer :: valteirk(:)

contains

    ! GENINTRK assigns entry i (starting at 1) to MOD(i, nprocs).
    pure integer function rk_local_count(n, myid, nprocs) result(nlocal)
        ! Arguments
        integer, intent(in) :: n, myid, nprocs
        if (n < 0 .or. nprocs < 1) error stop "RK_LOCAL_COUNT: invalid rank layout"
        if (myid < 0 .or. myid >= nprocs) error stop "RK_LOCAL_COUNT: invalid rank layout"
        nlocal = n/nprocs
        if (myid > 0 .and. myid <= mod(n,nprocs)) nlocal = nlocal + 1
    end function rk_local_count
end module distributed_slater_integrals
