module calculate_seltzer_field_shift_factors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_seltzer_field_shift_factors
contains

    !***********************************************************************
    !                                                                      *
    !JE   SUBROUTINE DENSNEW(DOIT,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
    subroutine calculate_seltzer_field_shift_factors(angular_cache_write_mode,contact_density_integrals, &
        & total_nms_integrals,radius_integrals, &
        radius_squared_integrals,inverse_radius_integrals,inverse_radius_squared_integrals,leading_nms_integrals, &
        orbital_density_profiles,level_density_profiles,nuclear_density_point_count)
        ! Module dependencies
        use calculation_file_units, only: one_body_cache_unit
        !                                                                      *
        !   IF angular coefficients must be calculated                         *
        !   This routine controls combines the radial and angular parts for the*
        !   calculation of the NMS parameter, the electron density at the      *
        !   origin and radial expectation values.
        !                                                                      *
        !   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
        !                        ITJPO, ONESCALAR                              *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !                                         Last revision: 10 Nov 1995   *
        !                                                                      *
        !   Modified by C. Naz\'e  Feb. 2012                                   *
        !   Modified by J. Ekman   Nov. 2013                                   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Gediminas Gaigalas 11/18/19
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: orbital_index_radix, max_orbitals
        use debug_output_settings, only: cutoff
        use eigenstate_expansions, only: csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_kappa
        use selected_level_indices, only: active_level_count
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        use integer_format, only: format_integer
        use yes_no_input, only: read_yes_no
        use csf_quantum_queries, only: csf_total_j_code
        use one_electron_scalar_coefficients, only: calculate_one_electron_scalar_coefficients
        implicit none (type, external)

        ! Local declarations
        integer :: l, level_mixing_offset, coefficient_count, packed_orbital_label

        ! Arguments
        integer, intent(in) :: nuclear_density_point_count
        real(real64), intent(in) :: contact_density_integrals(max_orbitals, max_orbitals), &
            & total_nms_integrals(max_orbitals, max_orbitals), &
            radius_integrals(max_orbitals, max_orbitals), radius_squared_integrals(max_orbitals, max_orbitals), &
            inverse_radius_integrals(max_orbitals, max_orbitals), inverse_radius_squared_integrals(max_orbitals, max_orbitals), &
            leading_nms_integrals(max_orbitals, max_orbitals)
        real(real64), intent(out)     :: level_density_profiles(active_level_count, nuclear_density_point_count) !  JE ADD
        real(real64), intent(in) :: orbital_density_profiles(max_orbitals, max_orbitals, nuclear_density_point_count) !  JE ADD
        integer, intent(in) :: angular_cache_write_mode

        ! Local constants

        ! Local declarations
        real(real64) :: shell_coefficients(max_orbitals), stored_shell_coefficients(max_orbitals)
        real(real64) :: contri1vec(nuclear_density_point_count), elemnt1vec(nuclear_density_point_count)
        real(real64) :: contact_density_element, total_nms_element, radius_element, radius_squared_element, &
            & inverse_radius_element,     &
            inverse_radius_squared_element, leading_nms_element, contact_density_contribution, total_nms_contribution, &
                & radius_contribution,     &
            radius_squared_contribution, inverse_radius_contribution, inverse_radius_squared_contribution, leading_nms_contribution
        character :: number_text*11
        integer :: contributing_orbitals(max_orbitals)
        integer :: operator_rank, operator_parity, include_closed_shell_terms, csf_index, text_length, &
            & column_two_j_plus_one, row_csf, orbital_index, &
            ket_orbital_index, i, level_index
        !-----------------------------------------------
        !
        ! DOIT: IF DOIT=1 angular coefficients will be stored after creation
        !
        level_density_profiles(:,:) = 0.0_real64                                     ! JE ADD
        !
        !   Set the rank (zero) and parity (even) for the one-particle
        !   coefficients
        !
        operator_rank = 0
        operator_parity = 1
        include_closed_shell_terms = 1
        !
        !   Allocate storage for the arrays in BUFFER
        !
        call manage_angular_coefficient_buffer (1)
        !
        !   Sweep through the Hamiltonian matrix to determine the
        !   sms parameter
        !
        do csf_index = 1,csf_count
            !
            !   Output IC on the screen to show how far the calculation has preceede
            !
            call format_integer (csf_index,number_text,text_length)
            if (mod(csf_index,100)==0) then
                print *, "Column "//number_text(1:text_length)//" complete;"
            end if
            !
            column_two_j_plus_one = csf_total_j_code (csf_index)

            do row_csf = csf_index,csf_count
                !
                !   Matrix elements are diagonal in J
                !
                if (csf_total_j_code(row_csf) == column_two_j_plus_one) then
                    !
                    !   Initialise the accumulator
                    !
                    contact_density_element = 0.0_real64
                    total_nms_element = 0.0_real64
                    radius_element = 0.0_real64
                    radius_squared_element = 0.0_real64
                    inverse_radius_element = 0.0_real64
                    inverse_radius_squared_element = 0.0_real64
                    leading_nms_element = 0.0_real64
                    elemnt1vec(:) = 0.0_real64                                      ! JE ADD
                    !
                    !   Call the MCT package to compute T coefficients
                    !
                    call calculate_one_electron_scalar_coefficients(csf_index,row_csf,orbital_index,ket_orbital_index, &
                        & shell_coefficients)
                    !GG            CALL TNSRJJ (KA,IOPAR,IC,IR,IA,IB,TSHELL)
                    if (orbital_index /= 0) then
                        if (orbital_index == ket_orbital_index) then
                            coefficient_count = 0
                            do orbital_index = 1,orbital_count
                                if (abs (shell_coefficients(orbital_index)) > cutoff) then
                                    coefficient_count = coefficient_count + 1
                                    stored_shell_coefficients(coefficient_count) = shell_coefficients(orbital_index)
                                    contributing_orbitals(coefficient_count) = orbital_index
                                    contact_density_element = contact_density_element + &
                                        & contact_density_integrals(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                                    do l = 2,nuclear_density_point_count                                    ! JE ADD
                                        elemnt1vec(l) = elemnt1vec(l)+orbital_density_profiles(orbital_index,orbital_index,l)*& ! JE ADD
                                            shell_coefficients(orbital_index)                        ! JE ADD
                                    end do  ! JE ADD
                                    total_nms_element = total_nms_element + total_nms_integrals(orbital_index, &
                                        & orbital_index)*shell_coefficients(orbital_index)
                                    radius_element = radius_element + radius_integrals(orbital_index,orbital_index)* &
                                        & shell_coefficients(orbital_index)
                                    radius_squared_element = radius_squared_element + &
                                        & radius_squared_integrals(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                                    inverse_radius_element = inverse_radius_element + &
                                        & inverse_radius_integrals(orbital_index,orbital_index)*shell_coefficients(orbital_index)
                                    inverse_radius_squared_element = inverse_radius_squared_element + &
                                        & inverse_radius_squared_integrals(orbital_index,orbital_index)* &
                                        & shell_coefficients(orbital_index)
                                    leading_nms_element = leading_nms_element + leading_nms_integrals(orbital_index, &
                                        & orbital_index)*shell_coefficients(orbital_index)
                                end if
                            end do
                            if (angular_cache_write_mode==1) write(one_body_cache_unit) csf_index,row_csf,coefficient_count
                            do i = 1,coefficient_count
                                packed_orbital_label = contributing_orbitals(i)*(orbital_index_radix + 1)
                                if (angular_cache_write_mode==1) write(one_body_cache_unit) stored_shell_coefficients(i), &
                                    & packed_orbital_label
                            end do
                        else
                            if (abs (shell_coefficients(1)) > cutoff) then
                                if (orbital_kappa(orbital_index)==orbital_kappa(ket_orbital_index)) then
                                    if (angular_cache_write_mode==1) write(one_body_cache_unit) csf_index,row_csf,1
                                    packed_orbital_label = orbital_index*orbital_index_radix + ket_orbital_index
                                    if (angular_cache_write_mode==1) write(one_body_cache_unit) shell_coefficients(1), &
                                        & packed_orbital_label
                                    contact_density_element = contact_density_element + &
                                        & contact_density_integrals(orbital_index,ket_orbital_index)*shell_coefficients(1)
                                    do l = 2,nuclear_density_point_count                                   ! JE ADD
                                        elemnt1vec(l) = elemnt1vec(l)+orbital_density_profiles(orbital_index,ket_orbital_index,l)*& ! JE ADD
                                            shell_coefficients(1)                         ! JE ADD
                                    end do  ! JE ADD
                                    total_nms_element = total_nms_element + total_nms_integrals(orbital_index, &
                                        & ket_orbital_index)*shell_coefficients(1)
                                    radius_element = radius_element + radius_integrals(orbital_index,ket_orbital_index)* &
                                        & shell_coefficients(1)
                                    radius_squared_element = radius_squared_element + &
                                        & radius_squared_integrals(orbital_index,ket_orbital_index)*shell_coefficients(1)
                                    inverse_radius_element = inverse_radius_element + &
                                        & inverse_radius_integrals(orbital_index,ket_orbital_index)*shell_coefficients(1)
                                    inverse_radius_squared_element = inverse_radius_squared_element + &
                                        & inverse_radius_squared_integrals(orbital_index,ket_orbital_index)*shell_coefficients(1)
                                    leading_nms_element = leading_nms_element + leading_nms_integrals(orbital_index, &
                                        & ket_orbital_index)*shell_coefficients(1)
                                end if
                            end if
                        end if
                    end if
                    do level_index = 1,active_level_count
                        level_mixing_offset = (level_index-1)*csf_count
                        contact_density_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+ &
                            & level_mixing_offset)*contact_density_element
                        contri1vec(:) = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+level_mixing_offset)*elemnt1vec(:)   ! JE ADD
                        total_nms_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+level_mixing_offset)* &
                            & total_nms_element
                        radius_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+level_mixing_offset)*radius_element
                        radius_squared_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+ &
                            & level_mixing_offset)*radius_squared_element
                        inverse_radius_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+ &
                            & level_mixing_offset)*inverse_radius_element
                        inverse_radius_squared_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+level_mixing_offset)*inverse_radius_squared_element
                        leading_nms_contribution = csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                            & csf_mixing_coefficients(row_csf+level_mixing_offset)* &
                            & leading_nms_element
                        if (row_csf/=csf_index) then
                            contact_density_contribution = 2.0_real64 * contact_density_contribution
                            contri1vec(:) = 2.0_real64 * contri1vec(:)                  ! JE ADD
                            total_nms_contribution = 2.0_real64 * total_nms_contribution
                            radius_contribution = 2.0_real64 * radius_contribution
                            radius_squared_contribution = 2.0_real64 * radius_squared_contribution
                            inverse_radius_contribution = 2.0_real64 * inverse_radius_contribution
                            inverse_radius_squared_contribution = 2.0_real64 * inverse_radius_squared_contribution
                            leading_nms_contribution = 2.0_real64 * leading_nms_contribution
                        end if
                        contact_density_by_level(level_index) = contact_density_by_level(level_index) + contact_density_contribution
                        do l = 2,nuclear_density_point_count                                         ! JE ADD
                            level_density_profiles(level_index,l) = level_density_profiles(level_index,l) + contri1vec(l)           ! JE ADD
                        end do  ! JE ADD
                        total_nms_by_level(level_index) = total_nms_by_level(level_index) + total_nms_contribution
                        radius_by_level(level_index) = radius_by_level(level_index) + radius_contribution
                        radius_squared_by_level(level_index) = radius_squared_by_level(level_index) + radius_squared_contribution
                        inverse_radius_by_level(level_index) = inverse_radius_by_level(level_index) + inverse_radius_contribution
                        inverse_radius_squared_by_level(level_index) = inverse_radius_squared_by_level(level_index) + &
                            & inverse_radius_squared_contribution
                        leading_nms_by_level(level_index) = leading_nms_by_level(level_index) + leading_nms_contribution
                    end do
                end if
            end do
        end do
        if (angular_cache_write_mode==1) write(one_body_cache_unit) -1
        !
        ! Empty the buffer and close file
        if (angular_cache_write_mode==1) close(one_body_cache_unit)
        !
        !   Deallocate storage for the arrays in BUFFER
        call manage_angular_coefficient_buffer (3)
        return
    end subroutine calculate_seltzer_field_shift_factors
end module calculate_seltzer_field_shift_factors_interface
