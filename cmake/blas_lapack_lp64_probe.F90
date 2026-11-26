! Use the same Fortran interfaces and imported libraries as the Suite targets.
program blas_lapack_lp64_probe
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use blas_lapack_interfaces, only: dgemv, ddot
#ifdef GRASP_LAPACK_PROBE
    use blas_lapack_interfaces, only: dspevx, slamch, dlamch
#endif
    implicit none (type, external)

    if (storage_size(0) /= 32 .or. storage_size(0.0_real32) /= 32 .or. storage_size(0.0_real64) /= 64) &
        error stop "Suite requires 32-bit default integers and 32/64-bit real kinds"

#ifdef GRASP_LAPACK_PROBE
    call check_lapack()
#endif
    call check_blas()
    print *, "Suite LP64 ABI verified"

contains

    subroutine check_blas()
        ! Each scalar has space for a full 64-bit access, with the second word fixed.
        ! With LP64 this is a 1-by-1 multiply. With ILP64 the row count is negative
        ! on little endian, or greater than the leading dimension on big endian.
        ! DGEMV checks these parameters before using the array arguments.
        integer, volatile, allocatable :: args(:,:)
        real(real64) :: a(1,1), x(1), y(1)

        allocate(args(2,5))
        args(1,:) = 1
        args(2,:) = 0
        args(2,1) = -1
        a = 2.0_real64
        x = 3.0_real64
        y = 0.0_real64
        call dgemv("N", args(1,1), args(1,2), 1.0_real64, a, args(1,3), x, args(1,4), &
            0.0_real64, y, args(1,5))
        if (y(1) /= 6.0_real64) error stop "DGEMV did not accept the LP64 argument layout"
        ! Use the same padded leading words for this scalar return check.
        x = 2.0_real64
        if (ddot(args(1,2), x, args(1,4), y, args(1,5)) /= 12.0_real64) error stop "DDOT return ABI"
    end subroutine check_blas

#ifdef GRASP_LAPACK_PROBE
    subroutine check_lapack()
        integer, volatile, allocatable :: args(:,:)
        integer :: iwork(2), ifail(2)
        real(real64) :: ap(1), w(1), z(1,1), work(1)

        ! N=0 returns without touching the work arrays. M and INFO should both be
        ! written as zero without changing the next word of either integer slot.
        allocate(args(2,6))
        args = 0
        args(1,4) = 1
        args(:,5:6) = -1
        call dspevx("N", "A", "U", args(1,1), ap, 0.0_real64, 0.0_real64, args(1,2), args(1,3), 0.0_real64, &
            args(1,5), w, z, args(1,4), work, iwork, ifail, args(1,6))
        if (any(args(2,5:6) /= -1)) error stop "DSPEVX wrote outside the 32-bit integer outputs"
        if (any(args(1,5:6) /= 0)) error stop "DSPEVX zero-order LP64 call failed"
        if (slamch("B") /= real(radix(0.0_real32),real32)) error stop "SLAMCH return ABI"
        if (dlamch("B") /= real(radix(0.0_real64),real64)) error stop "DLAMCH return ABI"
    end subroutine check_lapack
#endif
end program blas_lapack_lp64_probe
