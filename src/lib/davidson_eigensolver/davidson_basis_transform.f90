module davidson_basis_transform
    implicit none (type, external)

    ! Module interface
    private
    public :: multbc
contains
    subroutine multbc(row_count, input_column_count, output_column_count, coefficients, row_workspace, vectors)
        ! Module dependencies
        ! Right-multiply the vector matrix by the coefficients, overwriting its leading columns.
        use, intrinsic :: iso_fortran_env, only: real64
        use blas_lapack_interfaces, only: dcopy, dgemv
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: row_count
        integer, intent(in) :: input_column_count
        integer, intent(in) :: output_column_count
        real(real64)  :: coefficients(input_column_count*output_column_count)
        real(real64)  :: row_workspace(output_column_count)
        real(real64)  :: vectors(row_count*input_column_count)

        ! Local declarations
        integer :: row_index
        ! Compute the entire row before overwriting it. Successive columns are row_count elements apart.
        do row_index = 1, row_count
            call dgemv ("Transp", input_column_count, output_column_count, 1.0_real64, coefficients, &
                input_column_count, vectors(row_index), row_count, 0.0_real64, row_workspace, 1)
            call dcopy (output_column_count, row_workspace, 1, vectors(row_index), row_count)
        end do

        return
    end subroutine multbc
end module davidson_basis_transform
