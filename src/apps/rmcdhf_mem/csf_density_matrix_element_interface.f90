module csf_density_matrix_element_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function csf_density_matrix_element (extended_optimal_level, bra_csf_index, ket_csf_index, &
            & symmetry_block_index)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            logical, intent(in) :: extended_optimal_level
            integer, intent(in) :: bra_csf_index
            integer, intent(in) :: ket_csf_index
            integer :: symmetry_block_index
        end function csf_density_matrix_element
    end interface
end module csf_density_matrix_element_interface
