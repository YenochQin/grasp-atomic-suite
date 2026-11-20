! LP64 BLAS/LAPACK interfaces use default integers and the precision declared by each routine.
module blas_lapack_interfaces
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private
    public :: dcopy, ddot, daxpy, dscal, dgemv, dspevx, dswap, idamax
    public :: slamch, dlamch

    ! Procedure interfaces
    interface
        real(real32) function slamch(cmach)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real32
            implicit none (type, external)

            ! Arguments
            character(1), intent(in) :: cmach
        end function slamch
        function dlamch(cmach) result(value)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            character(1), intent(in) :: cmach
            real(real64) :: value
        end function dlamch
        subroutine dswap(n, x, incx, y, incy)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n, incx, incy
            real(real64), intent(inout) :: x(*), y(*)
        end subroutine dswap
        function idamax(n, x, incx) result(index)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n, incx
            real(real64), intent(in) :: x(*)
            integer :: index
        end function idamax
        subroutine dcopy(n, x, incx, y, incy)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n, incx, incy
            real(real64), intent(in) :: x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine dcopy
        function ddot(n, x, incx, y, incy) result(value)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n, incx, incy
            real(real64), intent(in) :: x(*), y(*)
            real(real64) :: value
        end function ddot
        subroutine daxpy(n, alpha, x, incx, y, incy)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n, incx, incy
            real(real64), intent(in) :: alpha, x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine daxpy
        subroutine dscal(n, alpha, x, incx)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n, incx
            real(real64), intent(in) :: alpha
            real(real64), intent(inout) :: x(*)
        end subroutine dscal
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            character(1), intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine dgemv
        subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info)
            ! Module dependencies
            import real64
            implicit none (type, external)

            ! Arguments
            character(1), intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, il, iu, ldz
            real(real64), intent(in) :: vl, vu, abstol
            real(real64), intent(inout) :: ap(*)
            integer, intent(out) :: m, info
            real(real64), intent(out) :: w(*), z(ldz,*), work(*)
            integer, intent(out) :: iwork(*), ifail(*)
        end subroutine dspevx
    end interface
end module blas_lapack_interfaces
