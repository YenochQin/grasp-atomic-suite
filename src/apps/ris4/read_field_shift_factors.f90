module read_field_shift_factors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_field_shift_factors
contains

    !       SUBROUTINE DENSREAD(DINT1,DINT2,DINT3,                           &
    !                           DINT4,DINT5,DINT6,                           &
    !                           DINT7,DINT1VEC,DENS1VEC)
    !
    !    IF angular coefficients already exist
    !    This routine controls combines the radial and angular parts for the*
    !    calculation of the NMS parameter, the electron density at the
    !    origin and radial expectation values.
    !
    ! Written by C\'edric Naz\'e
    !
    ! Last revision: Feb. 2011
    !
    ! Translated by Gediminas Gaigalas 11/18/19
    subroutine read_field_shift_factors(contact_density_integrals,total_nms_integrals,radius_integrals, &
        & radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals,leading_nms_integrals)
        ! Module dependencies
        use calculation_file_units, only: one_body_cache_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: orbital_index_radix, max_orbitals
        use selected_level_indices, only: active_level_count
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use orbital_csf_basis, only: csf_count, max_orbitals
        use eigenstate_expansions, only: csf_mixing_coefficients
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: contact_density_integrals(max_orbitals, max_orbitals), &
            & total_nms_integrals(max_orbitals, max_orbitals),  &
            radius_integrals(max_orbitals, max_orbitals), radius_squared_integrals(max_orbitals, max_orbitals), &
            inverse_radius_integrals(max_orbitals, max_orbitals), inverse_radius_squared_integrals(max_orbitals, max_orbitals),  &
            leading_nms_integrals(max_orbitals, max_orbitals)

        ! Local constants

        ! Local declarations
        real(real64) :: cached_shell_coefficients(max_orbitals)
        real(real64) :: contact_density_element, total_nms_element, radius_element, radius_squared_element, inverse_radius_element
        real(real64) :: inverse_radius_squared_element, leading_nms_element
        real(real64) :: contact_density_contribution, total_nms_contribution, radius_contribution, &
            & radius_squared_contribution, inverse_radius_contribution
        real(real64) :: inverse_radius_squared_contribution, leading_nms_contribution
        integer :: previous_column_csf, previous_row_csf, io_status, ia, ib, column_csf, row_csf, i, level_index, &
            & level_mixing_offset, packed_orbital_label
        integer :: coefficient_count

        previous_column_csf = 0
        previous_row_csf = 0

        rewind (one_body_cache_unit)
        do
            read (one_body_cache_unit,iostat = io_status) column_csf,row_csf,coefficient_count
            if (io_status /= 0) exit
            ! *      print*, 'ic', IC,IR,NCOUNT
            !
            !    Read successful; decode the labels of I(ab)
            !
            !       Initialise
            if ((column_csf/=previous_column_csf).or.(row_csf/=previous_row_csf)) then
                ! *            ISPARC = ISPAR (IC)
                ! *            ITJPOC = ITJPO (IC)
                ! *            ITJPOR = ITJPO (IR)
                ! *            IDIFF = ITJPOC - ITJPOR
                previous_column_csf = column_csf
                previous_row_csf = row_csf
                ! *        DO I = 1,NELMNT
                contact_density_element = 0.0_real64
                total_nms_element = 0.0_real64
                radius_element = 0.0_real64
                radius_squared_element = 0.0_real64
                inverse_radius_element = 0.0_real64
                inverse_radius_squared_element = 0.0_real64
                leading_nms_element = 0.0_real64
                ! *        ENDDO
            end if
            do i= 1,coefficient_count
                read(one_body_cache_unit) cached_shell_coefficients(i),packed_orbital_label
                ia = packed_orbital_label/orbital_index_radix
                ib = mod(packed_orbital_label,orbital_index_radix)
                ! *           LOC = INDEX(I)
                contact_density_element = contact_density_element + contact_density_integrals(ia,ib)*cached_shell_coefficients(i)
                total_nms_element = total_nms_element + total_nms_integrals(ia,ib)*cached_shell_coefficients(i)
                radius_element = radius_element + radius_integrals(ia,ib)*cached_shell_coefficients(i)
                radius_squared_element = radius_squared_element + radius_squared_integrals(ia,ib)*cached_shell_coefficients(i)
                inverse_radius_element = inverse_radius_element + inverse_radius_integrals(ia,ib)*cached_shell_coefficients(i)
                inverse_radius_squared_element = inverse_radius_squared_element + inverse_radius_squared_integrals(ia,ib)* &
                    & cached_shell_coefficients(i)
                leading_nms_element = leading_nms_element + leading_nms_integrals(ia,ib)*cached_shell_coefficients(i)
            end do
            ! Return to the start of the loop

            ! *      ICI = 0
            ! *      DO 21 I = 1,NELMNT
            ! *        IRI = IROW(I)
            ! *         IF (I .GT. IENDC(ICI)) ICI = ICI+1
            do level_index = 1,active_level_count
                level_mixing_offset = (level_index-1)*csf_count
                contact_density_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)* &
                    & contact_density_element
                total_nms_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)*total_nms_element
                radius_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)*radius_element
                radius_squared_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)*radius_squared_element
                inverse_radius_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)*inverse_radius_element
                inverse_radius_squared_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)* &
                    & inverse_radius_squared_element
                leading_nms_contribution = csf_mixing_coefficients(column_csf+level_mixing_offset)* &
                    & csf_mixing_coefficients(row_csf+level_mixing_offset)*leading_nms_element
                if (row_csf/=column_csf) then
                    contact_density_contribution = 2.0_real64 * contact_density_contribution
                    total_nms_contribution = 2.0_real64 * total_nms_contribution
                    radius_contribution = 2.0_real64 * radius_contribution
                    radius_squared_contribution = 2.0_real64 * radius_squared_contribution
                    inverse_radius_contribution = 2.0_real64 * inverse_radius_contribution
                    inverse_radius_squared_contribution = 2.0_real64 * inverse_radius_squared_contribution
                    leading_nms_contribution = 2.0_real64 * leading_nms_contribution
                end if
                contact_density_by_level(level_index) = contact_density_by_level(level_index) + contact_density_contribution
                total_nms_by_level(level_index) = total_nms_by_level(level_index) + total_nms_contribution
                radius_by_level(level_index) = radius_by_level(level_index) + radius_contribution
                radius_squared_by_level(level_index) = radius_squared_by_level(level_index) + radius_squared_contribution
                inverse_radius_by_level(level_index) = inverse_radius_by_level(level_index) + inverse_radius_contribution
                inverse_radius_squared_by_level(level_index) = inverse_radius_squared_by_level(level_index) + &
                    & inverse_radius_squared_contribution
                leading_nms_by_level(level_index) = leading_nms_by_level(level_index) + leading_nms_contribution
            end do
            ! *   21 CONTINUE

        end do
        return
    end subroutine read_field_shift_factors
end module read_field_shift_factors_interface
