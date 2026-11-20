module locate_mcp_restart_position_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine locate_mcp_restart_position (first_bra_csf_index, first_ket_csf_index, matrix_element_position, &
            & one_electron_coefficient_count, &
            & two_electron_counts_by_rank)
            implicit none (type, external)

            ! Arguments
            integer, intent(inout) :: first_bra_csf_index
            integer, intent(inout) :: first_ket_csf_index
            integer, intent(inout) :: matrix_element_position
            integer, intent(out) :: one_electron_coefficient_count
            integer, pointer :: two_electron_counts_by_rank(:)
        end subroutine locate_mcp_restart_position
    end interface
end module locate_mcp_restart_position_interface
