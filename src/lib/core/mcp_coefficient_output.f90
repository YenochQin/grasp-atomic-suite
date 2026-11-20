! Append one labeled angular coefficient to the shared MCP buffer, growing it as needed.
module mcp_coefficient_output
    implicit none (type, external)

    ! Module interface
    private
    public :: speak
contains
    subroutine speak(bra_csf_index, ket_csf_index, bra_orbital_1, bra_orbital_2, ket_orbital_1, ket_orbital_2, &
        & tensor_rank, coefficient)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use angular_coefficient_buffer_state,        only: &
            angular_coefficient_labels, angular_coefficients, angular_buffer_capacity, &
            & angular_coefficient_count
        use debug_output_settings ,        only: ibug1
        use orbital_csf_basis,           only: orbital_principal_n, orbital_symmetry_labels
        use angular_coefficient_buffer, only: alcbuf
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index
        integer, intent(in) :: ket_csf_index
        integer, intent(in) :: bra_orbital_1
        integer, intent(in) :: bra_orbital_2
        integer, intent(in) :: ket_orbital_1
        integer, intent(in) :: ket_orbital_2
        integer, intent(in) :: tensor_rank
        real(real64), intent(in) :: coefficient
        if (ibug1 /= 0) then
            call ensure_debug_output()
            write (angular_debug_unit, '(2(1x,1i2),4(1x,i2,a2),1x,i2,1x,1p,d19.12)') &
                bra_csf_index, ket_csf_index, orbital_principal_n(bra_orbital_1), &
                & orbital_symmetry_labels(bra_orbital_1), &
                & orbital_principal_n(bra_orbital_2), orbital_symmetry_labels(bra_orbital_2&
                ), orbital_principal_n(ket_orbital_1), orbital_symmetry_labels(ket_orbital_1), &
                & orbital_principal_n(ket_orbital_2), orbital_symmetry_labels(ket_orbital_2), tensor_rank, coefficient
        end if
        angular_coefficient_count = angular_coefficient_count + 1
        if (angular_coefficient_count > angular_buffer_capacity) call alcbuf (2)
        angular_coefficient_labels(1,angular_coefficient_count) = bra_orbital_1
        angular_coefficient_labels(2,angular_coefficient_count) = bra_orbital_2
        angular_coefficient_labels(3,angular_coefficient_count) = ket_orbital_1
        angular_coefficient_labels(4,angular_coefficient_count) = ket_orbital_2
        angular_coefficient_labels(5,angular_coefficient_count) = tensor_rank
        angular_coefficients(angular_coefficient_count) = coefficient
        return
        return
    end subroutine speak
end module mcp_coefficient_output
