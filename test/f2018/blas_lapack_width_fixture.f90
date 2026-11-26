! Both fixture archives export the same Fortran symbols. Compile with the
! default integer or -fdefault-integer-8 to vary the ABI without a vendor library.
subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    character(1), intent(in) :: trans
    integer, intent(in) :: m, n, lda, incx, incy
    real(real64), intent(in) :: alpha, a(*), x(*), beta
    real(real64), intent(inout) :: y(*)

    if (trans /= "N" .or. any([m, n, lda, incx, incy] /= 1)) return
    y(1) = alpha*a(1)*x(1) + beta*y(1)
end subroutine dgemv

function ddot(n, x, incx, y, incy) result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer, intent(in) :: n, incx, incy
    real(real64), intent(in) :: x(*), y(*)
    real(real64) :: value

    value = 0.0_real64
    if (all([n, incx, incy] == 1)) value = x(1)*y(1)
end function ddot

subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    character(1), intent(in) :: jobz, range, uplo
    integer, intent(in) :: n, il, iu, ldz
    integer, intent(out) :: m, iwork(*), ifail(*), info
    real(real64), intent(in) :: vl, vu, abstol
    real(real64), intent(inout) :: ap(*)
    real(real64), intent(out) :: w(*), z(*), work(*)

    if (n /= 0) error stop "Width fixture only supports DSPEVX with N=0"
    m = 0
    info = 0
end subroutine dspevx

function slamch(cmach) result(value)
    use, intrinsic :: iso_fortran_env, only: real32
    implicit none (type, external)
    character(1), intent(in) :: cmach
    real(real32) :: value

    if (cmach /= "B") error stop "Width fixture only supplies the floating-point radix"
    value = real(radix(0.0_real32),real32)
end function slamch

function dlamch(cmach) result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    character(1), intent(in) :: cmach
    real(real64) :: value

    if (cmach /= "B") error stop "Width fixture only supplies the floating-point radix"
    value = real(radix(0.0_real64),real64)
end function dlamch
