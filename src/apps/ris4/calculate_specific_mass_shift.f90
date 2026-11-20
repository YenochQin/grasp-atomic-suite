module calculate_specific_mass_shift_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_specific_mass_shift
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_specific_mass_shift(angular_cache_write_mode,vinti_integrals,recoil_corrected_vinti_integrals)
        ! Module dependencies
        use calculation_file_units, only: two_body_cache_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   This routine controls the main sequence of routine calls for the   *
        !   calculation  of the  sms parameter, the electron density at the    *
        !   origin.
        !                                                                      *
        !   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
        !                        ITJPO,                                        *
        !               [SMS92]: RINTISO, RINTDENS, VINTI                      *
        !               [LIBRANG]: RKCO_GG                                     *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !                                         Last revision: 10 Nov 1995   *
        !                                                                      *
        !   Modified by C. Naz\'e  Mai. 2012                                   *
        !                                                                      *
        !***********************************************************************
        !...Translated by Gediminas Gaigalas 11/18/19
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: orbital_index_radix, max_orbitals
        use debug_output_settings, only: cutoff, ldbpa
        use selected_level_indices, only: active_level_count
        use orbital_csf_basis, only: csf_count, max_orbitals, orbital_principal_n, orbital_symmetry_labels
        use angular_coefficient_buffer_state, only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
        use eigenstate_expansions, only: csf_mixing_coefficients
        use isotope_shift_options, only: contact_density_by_level, total_nms_by_level, radius_by_level, &
            & radius_squared_by_level, inverse_radius_by_level, inverse_radius_squared_by_level, leading_nms_by_level, &
            & leading_sms_by_level, total_sms_by_level
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer
        use csf_quantum_queries, only: csf_total_j_code
        use integer_format, only: format_integer
        use two_electron_angular_coefficients, only: calculate_two_electron_angular_coefficients
        use closed_shell_mcp, only: calculate_closed_shell_mcp
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: vinti_integrals(max_orbitals, max_orbitals), &
            & recoil_corrected_vinti_integrals(max_orbitals, max_orbitals)
        integer, intent(in) :: angular_cache_write_mode
        !-----------------------------------------------
        !  E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------

        ! Local constants

        ! Local declarations
        character(len=11) :: number_text
        real(real64) :: total_sms_contribution, leading_sms_contribution, angular_coefficient
        integer      :: i,level_index,k,packed_orbital_label,level_mixing_offset,csf_index,row_csf,column_two_j_plus_one, &
            & orbital_a,orbital_b,orbital_c,orbital_d,include_closed_shell_terms
        integer      :: text_length
        !-----------------------------------------------
        !
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
                !   Call the MCP package to generate V coefficients; ac and bd
                !   are the density pairs
                !
                !   Initialize
                !
                !
                !   Matrix elements are diagonal in J
                !
                if (csf_total_j_code(row_csf) == column_two_j_plus_one) then
                    angular_coefficient_count = 0
                    !GG            CALL RKCO (IC,IR,COR,CORD,INCOR)
                    call calculate_two_electron_angular_coefficients (csf_index, row_csf, calculate_closed_shell_mcp, &
                        & include_closed_shell_terms, 1)
                    !
                    do i = 1,angular_coefficient_count
                        angular_coefficient = angular_coefficients(i)
                        if (abs (angular_coefficient) > cutoff) then
                            orbital_a = angular_coefficient_labels(1,i)
                            orbital_b = angular_coefficient_labels(2,i)
                            orbital_c = angular_coefficient_labels(3,i)
                            orbital_d = angular_coefficient_labels(4,i)
                            k   = angular_coefficient_labels(5,i)
                            !
                            !   Only K = 1   LABEL(5,I) .EQ. 1
                            !
                            if (angular_coefficient_labels(5,i) == 1) then
                                if (ldbpa(2)) then
                                    call ensure_debug_output()
                                    write (angular_debug_unit, &
                                        & '(" V^[(",1i2,")]_[",1i3,",",1i3,"]"," (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,&
                                        &1a2,") = ",1pd19.12)') k,csf_index,row_csf,                            &
                                        orbital_principal_n(orbital_a),orbital_symmetry_labels(orbital_a), &
                                            & orbital_principal_n(orbital_b), &
                                        & orbital_symmetry_labels(orbital_b),                   &
                                        orbital_principal_n(orbital_c),orbital_symmetry_labels(orbital_c), &
                                            & orbital_principal_n(orbital_d), &
                                        & orbital_symmetry_labels(orbital_d),angular_coefficient
                                end if
                                if(angular_cache_write_mode==1) write(two_body_cache_unit) csf_index,row_csf
                                !** Storage sequence
                                packed_orbital_label  = ((orbital_a*orbital_index_radix + orbital_c)*orbital_index_radix+ &
                                    & orbital_b)*orbital_index_radix+orbital_d
                                if(angular_cache_write_mode==1) then
                                    write(two_body_cache_unit) angular_coefficient,packed_orbital_label
                                end if
                                !***
                                do level_index = 1,active_level_count
                                    level_mixing_offset = (level_index-1)*csf_count
                                    leading_sms_contribution = - csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                                        & csf_mixing_coefficients(row_csf+level_mixing_offset) &
                                        * angular_coefficient                                    &
                                        * vinti_integrals (angular_coefficient_labels(1,i),angular_coefficient_labels(3, &
                                            & i))              &
                                        * vinti_integrals (angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))
                                    total_sms_contribution = - csf_mixing_coefficients(csf_index+level_mixing_offset)* &
                                        & csf_mixing_coefficients(row_csf+level_mixing_offset) &
                                        * angular_coefficient                                    &
                                        * ( recoil_corrected_vinti_integrals(angular_coefficient_labels(1,i), &
                                            & angular_coefficient_labels(3,i))            &
                                        * vinti_integrals (angular_coefficient_labels(2,i),angular_coefficient_labels(4, &
                                            & i))              &
                                        + recoil_corrected_vinti_integrals(angular_coefficient_labels(2,i), &
                                            & angular_coefficient_labels(4,i))              &
                                        * vinti_integrals (angular_coefficient_labels(1,i),angular_coefficient_labels(3, &
                                            & i)) )/2.0_real64
                                    if (row_csf/=csf_index) then
                                        total_sms_contribution = 2.0_real64 * total_sms_contribution
                                        leading_sms_contribution = 2.0_real64 * leading_sms_contribution
                                    end if
                                    leading_sms_by_level(level_index) = leading_sms_by_level(level_index) + leading_sms_contribution
                                    total_sms_by_level(level_index) = total_sms_by_level(level_index) + total_sms_contribution
                                end do
                            end if
                        end if
                    end do
                end if
            end do
        end do
        if(angular_cache_write_mode==1) write(two_body_cache_unit) -1
        !
        ! Empty the buffer and close file
        if(angular_cache_write_mode==1) close(two_body_cache_unit)
        !
        !   Deallocate storage for the arrays in BUFFER
        !
        call manage_angular_coefficient_buffer (3)
        return
        !
    end subroutine calculate_specific_mass_shift
end module calculate_specific_mass_shift_interface
