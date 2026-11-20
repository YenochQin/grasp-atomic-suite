module breit_angular_coefficient_buffer
    implicit none (type, external)

    ! Module interface
    private
    public :: talk, talk_fallback

    ! Procedure interfaces
    ! Derive the external interface from the fallback. RCI supplies talk at link time.
    procedure(talk_fallback) :: talk
contains

    ! Store every input coefficient, including zero, with its radial-integral labels.
    subroutine talk_fallback(bra_csf_index, ket_csf_index, radial_integral_order, bra_orbital_1, bra_orbital_2, &
        ket_orbital_1, ket_orbital_2, integral_type_code, coefficient)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
            angular_coefficient_labels, angular_coefficients
        use debug_output_settings, only: angular_coefficient_debug_flag => ibug1
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use angular_coefficient_buffer, only: alcbuf
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index, ket_csf_index, radial_integral_order
        integer, intent(in) :: bra_orbital_1, bra_orbital_2, ket_orbital_1, ket_orbital_2, integral_type_code
        real(real64), intent(in) :: coefficient

        if (angular_coefficient_debug_flag /= 0) then
            call ensure_debug_output()
            write (angular_debug_unit, '(2(1x,1i2),4(1x,i2,a2),1x,1i2,1x,1i2,1x,1p,d19.12)') &
                bra_csf_index, ket_csf_index, &
                orbital_principal_n(bra_orbital_1), orbital_symmetry_labels(bra_orbital_1), &
                orbital_principal_n(bra_orbital_2), orbital_symmetry_labels(bra_orbital_2), &
                orbital_principal_n(ket_orbital_1), orbital_symmetry_labels(ket_orbital_1), &
                orbital_principal_n(ket_orbital_2), orbital_symmetry_labels(ket_orbital_2), &
                radial_integral_order, integral_type_code, coefficient
        end if

        angular_coefficient_count = angular_coefficient_count + 1
        ! Buffer operation 2 doubles capacity while preserving existing entries.
        if (angular_coefficient_count > angular_buffer_capacity) call alcbuf (2)

        ! Consumers read the four orbital indices, then the order and signed integral type.
        angular_coefficient_labels(1,angular_coefficient_count) = bra_orbital_1
        angular_coefficient_labels(2,angular_coefficient_count) = bra_orbital_2
        angular_coefficient_labels(3,angular_coefficient_count) = ket_orbital_1
        angular_coefficient_labels(4,angular_coefficient_count) = ket_orbital_2
        angular_coefficient_labels(5,angular_coefficient_count) = radial_integral_order
        angular_coefficient_labels(6,angular_coefficient_count) = integral_type_code
        angular_coefficients(angular_coefficient_count) = coefficient

        return
        return
    end subroutine talk_fallback
end module breit_angular_coefficient_buffer
