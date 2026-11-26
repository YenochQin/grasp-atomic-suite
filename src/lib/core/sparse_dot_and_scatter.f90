! Accumulate a sparse dot product and scatter the scaled entries into a second vector.
module sparse_dot_and_scatter
    implicit none (type, external)

    ! Module interface
    private
    public :: sparse_dot_scatter
contains
    subroutine sparse_dot_scatter(element_count, source_vector, updated_vector, target_indices, matrix_values, update_scale, &
        & dot_product)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64, int64
        implicit none (type, external)

        ! Arguments
        integer(int64), intent(in) :: element_count
        real(real64), intent(in)    :: update_scale
        real(real64), intent(out)   :: dot_product
        integer, intent(in) :: target_indices(element_count)
        real(real64), intent(in) :: source_vector(*)
        real(real64), intent(inout) :: updated_vector(*)
        real(real64), intent(in) :: matrix_values(element_count)

        ! Local declarations
        integer(int64) :: entry_index
        real(real64) :: dot_product_sum

        dot_product_sum = 0.0
        ! Reuse each stored entry for a local dot product and a scattered symmetric contribution.
        do entry_index = 1, element_count
            dot_product_sum = dot_product_sum + matrix_values(entry_index)*source_vector(target_indices(entry_index))
            updated_vector(target_indices(entry_index)) = updated_vector(target_indices(entry_index)) + &
                & update_scale*matrix_values(entry_index)
        end do
        dot_product = dot_product_sum

        return
    end subroutine sparse_dot_scatter
end module sparse_dot_and_scatter
