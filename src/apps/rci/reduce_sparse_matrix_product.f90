module reduce_sparse_matrix_product_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: reduce_sparse_matrix_product
contains

    !***********************************************************************
    !                                                                      *
    subroutine reduce_sparse_matrix_product(n, db, dc, da, dconst, dl)
        ! Module dependencies
        !
        !  Used by dnimcv
        !  Developed from dmerge. The only diff is: idy not needed here
        !
        !***********************************************************************
        !                     Gediminas Gaigalas  10/057
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: n
        real(real64), intent(in) :: dconst
        real(real64), intent(out) :: dl
        real(real64), intent(in) :: db(*)
        real(real64), intent(inout) :: dc(*)
        real(real64), intent(in) :: da(n)

        ! Local declarations
        real(real64) :: dsum
        !-----------------------------------------------

        dsum = 0.0_real64
        dsum = dot_product(da,db(:n))
        dc(:n) = dc(:n) + dconst*da
        dl = dsum

        return
    end subroutine reduce_sparse_matrix_product
end module reduce_sparse_matrix_product_interface
