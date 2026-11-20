module reduce_sparse_matrix_product_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine reduce_sparse_matrix_product (n, db, dc, da, dconst, dl)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n
            real(real64), intent(in) :: db(*)
            real(real64), intent(inout) :: dc(*)
            real(real64), intent(in) :: da(n)
            real(real64), intent(in) :: dconst
            real(real64), intent(out) :: dl
        end subroutine reduce_sparse_matrix_product
    end interface
end module reduce_sparse_matrix_product_interface
