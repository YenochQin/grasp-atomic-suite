module dvd_mpi_I
    implicit none
    interface
        subroutine ddotmpi(n, x, y, ddot)
            use vast_kind_param, only: double
            integer :: n
            real(double) :: x(n), y(n), ddot
        end subroutine
        subroutine dgemvmpi(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use vast_kind_param, only: double
            character(1) :: trans
            integer :: m, n, lda, incx, incy
            real(double) :: alpha, beta, a(lda,*), x(*), y(*)
        end subroutine
        subroutine multbc_col(n, k, m, c, temp, b)
            use vast_kind_param, only: double
            integer, intent(in) :: n
            integer :: k, m
            real(double) :: c(k*m), temp(m), b(n*k)
        end subroutine
    end interface
end module
