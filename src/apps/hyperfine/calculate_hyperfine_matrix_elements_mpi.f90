! Distribute hyperfine and magnetic matrix-element evaluation across MPI ranks.
module calculate_hyperfine_matrix_elements_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_hyperfine_matrix_elements_mpi
contains

    ! MPI version of HFSGG. Each rank processes a strided subset of the
    ! outer IC loop and the accumulated hyperfine constants are reduced
    ! over all MPI ranks before the root rank writes the output files.
    subroutine calculate_hyperfine_matrix_elements_mpi
        ! Module dependencies
        use calculation_file_units, only: hyperfine_output_unit, summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals
        use array_memory_management, only: alloc, dalloc
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation
        use atomic_calculation_state, only: vacuum_speed_of_light_au
        use hyperfine_nuclear_factors, only: calculate_nuclear_factors
        use eigenstate_expansions, only: csf_mixing_coefficients
        use first_order_csf_cutoffs, only: first_order_csf_cutoff
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use mpi_state,           only: mpi_rank, mpi_process_count, mpi_comm_world
        use isotope_properties, only: nuclear_spin
        use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one

        ! Procedure interfaces
        use hyperfine_radial_integral_interface, only: hyperfine_radial_integral

        ! Module dependencies
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap

        ! Procedure interfaces
        use hyperfine_one_electron_element_interface, only: hyperfine_one_electron_element

        ! Module dependencies
        use csf_quantum_queries, only: csf_parity_code, csf_total_j_code
        use angular_momentum_triangle, only: satisfies_angular_triangle
        use one_electron_tensor_coefficients, only: calculate_one_electron_tensor_coefficients
        use racah_coefficients_scott_fischer, only: scott_fischer_racah_coefficient
        use mpi_f08, only: mpi_count_kind
        use mpi_element_counts, only: checked_mpi_count
        use mpi_utilities, only: reduce_mpi_real_sum_to_root
        implicit none (type, external)

        ! Local constants
        real(real64), parameter :: contribution_magnitude_cutoff = 1.0e-10_real64

        ! Local declarations
        integer :: level_pair_count, hyperfine_value_count
        integer :: lowest_two_f, highest_two_f, two_f, first_orbital_index, second_orbital_index, tensor_rank, operator_parity, &
                & column_csf_index, row_csf_index, column_csf_j_code, &
            & row_csf_j_code&
            , csf_j_code_difference, bra_orbital_index, ket_orbital_index, first_level_index, second_level_index, &
                & first_mixing_offset, second_mixing_offset, second_output_level, first_output_level, two_nuclear_spin, &
                & first_level_two_j_plus_one, second_level_two_j_plus_one, first_recoupled_level, second_recoupled_level, &
                & first_level_two_j, second_level_two_j, wrote_matrix_heading, &
            column_csf_occupation, row_csf_occupation, differing_occupation_count, increased_occupation_orbital, &
                & decreased_occupation_orbital
        real(real64) :: shell_coefficients(max_orbitals)
        real(real64) :: hyperfine_radial_integrals(2, max_orbitals, max_orbitals), hyperfine_angular_factors(2, max_orbitals, &
                & max_orbitals), &
            hyperfine_orbital_elements(2, max_orbitals, max_orbitals)
        real(real64) :: lande_g_radial_integrals(max_orbitals, max_orbitals), radial_overlap_integrals(max_orbitals, &
            & max_orbitals), &
            lande_g_angular_factors(max_orbitals, max_orbitals), lande_g_correction_angular_factors(max_orbitals, max_orbitals), &
            lande_g_orbital_elements(max_orbitals, max_orbitals), lande_correction_orbital_elements(max_orbitals, max_orbitals)
        ! .. Local pointer arrays
        real(real64), pointer :: reduced_hyperfine_elements(:, :)
        real(real64), pointer :: reduced_lande_g_elements(:), reduced_lande_correction_elements(:)
        real(real64) ::  hyperfine_angular_factor, lande_g_angular_factor, lande_g_correction_angular_factor, &
                & hyperfine_csf_matrix_element,  &
            & lande_g_csf_matrix_element, lande_g_correction_csf_matrix_element,&
            hyperfine_contribution, lande_g_contribution, lande_g_correction_contribution, &
                & magnetic_dipole_mhz_scale, electric_quadrupole_mhz_scale, electronic_j, &
            inverse_j_norm, dipole_same_j_scale, dipole_delta_j_one_scale, quadrupole_same_j_scale, &
                & quadrupole_delta_j_one_scale, quadrupole_delta_j_two_scale, lande_g, lande_g_correction, reduced_nuclear_dipole, &
            reduced_nuclear_quadrupole, recoupling_phase, recoupling_prefactor, dipole_racah_coefficient, &
                & quadrupole_racah_coefficient, dipole_hyperfine_element, quadrupole_hyperfine_element, &
                & first_level_column_amplitude, first_level_row_amplitude,&
            second_level_column_amplitude, second_level_row_amplitude
        ! Matrix elements smaller than CUTOFF are not accumulated
        !
        level_pair_count = checked_mpi_count(int(active_level_count,mpi_count_kind)*active_level_count, "HYPERFINE level pairs")
        hyperfine_value_count = checked_mpi_count(5_mpi_count_kind*level_pair_count, "HYPERFINE matrix reduction")
        ! Allocate storage for local arrays
        call alloc (reduced_hyperfine_elements, 5, level_pair_count, "reduced_hyperfine_elements", &
            "calculate_hyperfine_matrix_elements_mpi")
        call alloc (reduced_lande_g_elements, level_pair_count, "reduced_lande_g_elements", &
            "calculate_hyperfine_matrix_elements_mpi")
        call alloc (reduced_lande_correction_elements, level_pair_count, "reduced_lande_correction_elements", &
            "calculate_hyperfine_matrix_elements_mpi")
        ! Initialise
        reduced_hyperfine_elements(:,:level_pair_count) = 0.0_real64
        reduced_lande_g_elements(:level_pair_count) = 0.0_real64
        reduced_lande_correction_elements(:level_pair_count) = 0.0_real64
        lande_g_orbital_elements(:,:) = 0.0_real64
        lande_correction_orbital_elements(:,:) = 0.0_real64

        ! Calculate and save the radial integrals and angular
        ! matrix elements for the two multipolarities
        do tensor_rank = 1, 2
            do first_orbital_index = 1, orbital_count
                do second_orbital_index = 1, orbital_count
                    if (tensor_rank == 1) then
                        hyperfine_radial_integrals(tensor_rank,first_orbital_index, &
                & second_orbital_index) = hyperfine_radial_integral(first_orbital_index,second_orbital_index,-2)
                        lande_g_radial_integrals(first_orbital_index, &
                & second_orbital_index) = hyperfine_radial_integral(first_orbital_index,second_orbital_index,1)
                    else
                        hyperfine_radial_integrals(tensor_rank,first_orbital_index, &
                & second_orbital_index) = integrate_weighted_radial_overlap(first_orbital_index,second_orbital_index,-3)
                        radial_overlap_integrals(first_orbital_index, &
                & second_orbital_index) = integrate_weighted_radial_overlap(first_orbital_index,second_orbital_index,0)
                    end if
                    call hyperfine_one_electron_element (first_orbital_index, tensor_rank, second_orbital_index, &
                & hyperfine_angular_factor, lande_g_angular_factor, &
                        & lande_g_correction_angular_factor)
                    hyperfine_angular_factors(tensor_rank,first_orbital_index,second_orbital_index) = hyperfine_angular_factor
                    hyperfine_orbital_elements(tensor_rank,first_orbital_index,second_orbital_index) = hyperfine_angular_factor* &
                & hyperfine_radial_integrals(tensor_rank,first_orbital_index,second_orbital_index)
                    if (tensor_rank /= 1) cycle
                    lande_g_angular_factors(first_orbital_index,second_orbital_index) = lande_g_angular_factor
                    lande_g_correction_angular_factors(first_orbital_index,second_orbital_index) = lande_g_correction_angular_factor
                end do
            end do
        end do
        do first_orbital_index = 1, orbital_count
            do second_orbital_index = 1, orbital_count
                lande_g_orbital_elements(first_orbital_index,second_orbital_index) = lande_g_angular_factors( &
                & first_orbital_index,second_orbital_index) &
                                    & *lande_g_radial_integrals(first_orbital_index,second_orbital_index)
                lande_correction_orbital_elements(first_orbital_index,second_orbital_index) = &
                    & lande_g_correction_angular_factors(first_orbital_index,second_orbital_index) &
                                    & *radial_overlap_integrals(first_orbital_index,second_orbital_index)
            end do
        end do
        ! Set the parity of the one-body operators
        operator_parity = 1
        ! Sweep through the Hamiltonian matrix to determine the
        ! diagonal and off-diagonal hyperfine constants
        do column_csf_index = mpi_rank + 1, csf_count, mpi_process_count
            column_csf_j_code = csf_total_j_code(column_csf_index)
            do row_csf_index = 1, csf_count
                ! If LFORDR is .TRUE., a `first order' calculation is indicated;
                ! only the CSFs with serial numbers exceeding IC are treated specially
                ! only diagonal elements are evaluated for the `first order' CSFs
                if (use_first_order_csf_perturbation .and. column_csf_index>first_order_csf_cutoff .and. column_csf_index/= &
                & row_csf_index) cycle
                row_csf_j_code = csf_total_j_code(row_csf_index)
                csf_j_code_difference = column_csf_j_code - row_csf_j_code
                ! Loop over the multipolarities
                do tensor_rank = 1, 2
                    ! Initialise the accumulator
                    hyperfine_csf_matrix_element = 0.0_real64
                    lande_g_csf_matrix_element = 0.0_real64
                    lande_g_correction_csf_matrix_element = 0.0_real64
                    ! Consider 3 cases
                    !             (k)
                    ! (1) < J || T   || J >    , k = 1,2  and IR >= IC
                    !
                    !             (k)
                    ! (2) < J || T   || J-1 >  , k = 1,2
                    !
                    !             (k)
                    ! (3) < J || T   || J-2 >  , k = 2
                    if (.not.(csf_j_code_difference==0 .and. row_csf_index>=column_csf_index .or. csf_j_code_difference==2 .or. &
                & csf_j_code_difference==4&
                        .and. tensor_rank==2)) cycle
                    if (satisfies_angular_triangle(column_csf_j_code, row_csf_j_code, 2*tensor_rank + 1) == 0) then
                        cycle
                    end if
                    if (operator_parity /= 0 .and. csf_parity_code(column_csf_index)*csf_parity_code(row_csf_index)* &
                & operator_parity /= 1) then
                        cycle
                    end if
                    differing_occupation_count = 0
                    increased_occupation_orbital = 0
                    decreased_occupation_orbital = 0
                    do first_orbital_index = 1, orbital_count
                        column_csf_occupation = csf_subshell_occupations(first_orbital_index,column_csf_index)
                        row_csf_occupation = csf_subshell_occupations(first_orbital_index,row_csf_index)
                        if (column_csf_occupation == row_csf_occupation) cycle
                        differing_occupation_count = differing_occupation_count + 1
                        if (differing_occupation_count > 2) exit
                        if (iabs(column_csf_occupation - row_csf_occupation) > 1) then
                            differing_occupation_count = 3
                            exit
                        end if
                        if (column_csf_occupation > row_csf_occupation) then
                            increased_occupation_orbital = first_orbital_index
                        else
                            decreased_occupation_orbital = first_orbital_index
                        end if
                    end do
                    if (differing_occupation_count /= 0 .and. differing_occupation_count /= 2) then
                        cycle
                    end if
                    if (differing_occupation_count == 2) then
                        if (abs(hyperfine_orbital_elements(tensor_rank,increased_occupation_orbital, &
                & decreased_occupation_orbital)) <= contribution_magnitude_cutoff) then
                            if (.not.(tensor_rank == 1 .and. csf_j_code_difference == 0 .and. &
                                (abs(lande_g_orbital_elements(increased_occupation_orbital, &
                & decreased_occupation_orbital)) > contribution_magnitude_cutoff .or. &
                                abs(lande_correction_orbital_elements(increased_occupation_orbital, &
                & decreased_occupation_orbital)) > contribution_magnitude_cutoff))) then
                                cycle
                            end if
                        end if
                    end if
                    call calculate_one_electron_tensor_coefficients(tensor_rank,operator_parity,column_csf_index,row_csf_index, &
                & bra_orbital_index, &
                        & ket_orbital_index,shell_coefficients)
                    !    Accumulate the contribution from the one-body operators;
                    if (bra_orbital_index /= 0) then
                        if (bra_orbital_index == ket_orbital_index) then
                            if (tensor_rank/=1 .or. csf_j_code_difference/=0) then
                                do bra_orbital_index = 1, orbital_count
                                    if (abs(shell_coefficients(bra_orbital_index)) <= contribution_magnitude_cutoff) cycle
                                    if (abs(hyperfine_orbital_elements(tensor_rank,bra_orbital_index, &
                & bra_orbital_index)) <= contribution_magnitude_cutoff) cycle
                                    hyperfine_csf_matrix_element = hyperfine_csf_matrix_element + &
                & hyperfine_orbital_elements(tensor_rank,bra_orbital_index,bra_orbital_index)*shell_coefficients(bra_orbital_index)
                                    cycle
                                end do
                            else
                                do bra_orbital_index = 1, orbital_count
                                    if (abs(shell_coefficients(bra_orbital_index)) <= contribution_magnitude_cutoff) cycle
                                    if (abs(hyperfine_orbital_elements(tensor_rank,bra_orbital_index, &
                & bra_orbital_index)) <= contribution_magnitude_cutoff .and. &
                                        abs(lande_g_orbital_elements(bra_orbital_index,bra_orbital_index)) <= &
                                            & contribution_magnitude_cutoff .and. &
                                        abs(lande_correction_orbital_elements(bra_orbital_index,bra_orbital_index)) <= &
                                            & contribution_magnitude_cutoff) cycle
                                    hyperfine_csf_matrix_element = hyperfine_csf_matrix_element + &
                & hyperfine_orbital_elements(tensor_rank,bra_orbital_index,bra_orbital_index)*shell_coefficients(bra_orbital_index)
                                    lande_g_csf_matrix_element = lande_g_csf_matrix_element + &
                & lande_g_orbital_elements(bra_orbital_index,bra_orbital_index)*shell_coefficients(bra_orbital_index)
                                    lande_g_correction_csf_matrix_element = lande_g_correction_csf_matrix_element + &
                                & lande_correction_orbital_elements(bra_orbital_index,bra_orbital_index)* &
                & shell_coefficients(bra_orbital_index)
                                end do
                            end if
                        else
                            if (abs(shell_coefficients(1)) > contribution_magnitude_cutoff) then
                                hyperfine_csf_matrix_element = hyperfine_csf_matrix_element + &
                & hyperfine_orbital_elements(tensor_rank,bra_orbital_index,ket_orbital_index)*shell_coefficients(1)
                                if (tensor_rank==1 .and. csf_j_code_difference==0) then
                                    lande_g_csf_matrix_element = lande_g_csf_matrix_element + &
                                & lande_g_orbital_elements(bra_orbital_index,ket_orbital_index)*shell_coefficients(1)
                                    lande_g_correction_csf_matrix_element = lande_g_correction_csf_matrix_element + &
                                & lande_correction_orbital_elements(bra_orbital_index,ket_orbital_index)*shell_coefficients(1)
                                end if
                            end if
                        end if
                    end if
                    ! Multiply with the configuration expansion coefficients and add the
                    ! contributions from the matrix elements to obtain total contributions
                    if (abs(hyperfine_csf_matrix_element) <= contribution_magnitude_cutoff .and. &
                        & abs(lande_g_csf_matrix_element) <= contribution_magnitude_cutoff .and. &
                        abs(lande_g_correction_csf_matrix_element) <= contribution_magnitude_cutoff) then
                        cycle
                    end if
                    do first_level_index = 1, active_level_count
                        first_mixing_offset = (first_level_index - 1)*csf_count
                        first_level_column_amplitude = csf_mixing_coefficients(column_csf_index + first_mixing_offset)
                        first_level_row_amplitude = csf_mixing_coefficients(row_csf_index + first_mixing_offset)
                        do second_level_index = 1, active_level_count
                            second_mixing_offset = (second_level_index - 1)*csf_count
                            second_level_column_amplitude = csf_mixing_coefficients(column_csf_index + second_mixing_offset)
                            second_level_row_amplitude = csf_mixing_coefficients(row_csf_index + second_mixing_offset)
                            if (csf_j_code_difference==0 .and. row_csf_index/=column_csf_index) then
                                hyperfine_contribution = hyperfine_csf_matrix_element*(first_level_column_amplitude* &
                & second_level_row_amplitude + first_level_row_amplitude*second_level_column_amplitude)
                                lande_g_contribution = lande_g_csf_matrix_element*(first_level_column_amplitude* &
                & second_level_row_amplitude + first_level_row_amplitude*second_level_column_amplitude)
                                lande_g_correction_contribution = lande_g_correction_csf_matrix_element* &
                & (first_level_column_amplitude*second_level_row_amplitude + &
                                & first_level_row_amplitude*second_level_column_amplitude)
                            else
                                hyperfine_contribution = hyperfine_csf_matrix_element*first_level_column_amplitude* &
                & second_level_row_amplitude
                                lande_g_contribution = lande_g_csf_matrix_element*first_level_column_amplitude* &
                & second_level_row_amplitude
                                lande_g_correction_contribution = lande_g_correction_csf_matrix_element* &
                & first_level_column_amplitude*second_level_row_amplitude
                            end if
                            ! Magnetic dipole and the two operators of the g_j factor
                            if (tensor_rank == 1) then
                                if (csf_j_code_difference == 0) then
                                    reduced_hyperfine_elements(1,active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_hyperfine_elements(1,active_level_count*(first_level_index-1)+ &
                & second_level_index) + hyperfine_contribution
                                    reduced_lande_g_elements(active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_lande_g_elements(active_level_count*(first_level_index-1)+second_level_index) + &
                & lande_g_contribution
                                    reduced_lande_correction_elements(active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_lande_correction_elements(active_level_count*(first_level_index-1)+ &
                & second_level_index) + &
                                        lande_g_correction_contribution
                                else if (column_csf_j_code - row_csf_j_code == 2) then
                                    reduced_hyperfine_elements(2,active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_hyperfine_elements(2,active_level_count*(first_level_index-1)+ &
                & second_level_index) + hyperfine_contribution
                                end if
                            ! Electric quadrupole
                            else if (tensor_rank == 2) then
                                select case (csf_j_code_difference)
                                case (0)
                                    reduced_hyperfine_elements(3,active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_hyperfine_elements(3,active_level_count*(first_level_index-1)+ &
                & second_level_index) + hyperfine_contribution
                                case (2)
                                    reduced_hyperfine_elements(4,active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_hyperfine_elements(4,active_level_count*(first_level_index-1)+ &
                & second_level_index) + hyperfine_contribution
                                case (4)
                                    reduced_hyperfine_elements(5,active_level_count*(first_level_index-1)+ &
                & second_level_index) = reduced_hyperfine_elements(5,active_level_count*(first_level_index-1)+ &
                & second_level_index) + hyperfine_contribution
                                end select
                            end if
                        end do
                    end do
                end do
            end do
        end do
        call reduce_mpi_real_sum_to_root(reduced_hyperfine_elements(:,:), hyperfine_value_count)
        call reduce_mpi_real_sum_to_root(reduced_lande_g_elements(:), level_pair_count)
        call reduce_mpi_real_sum_to_root(reduced_lande_correction_elements(:), level_pair_count)
        if (mpi_rank == 0) then
            call calculate_nuclear_factors(magnetic_dipole_mhz_scale, electric_quadrupole_mhz_scale, &
                reduced_nuclear_dipole, reduced_nuclear_quadrupole)
            ! Output the hyperfine interaction constants
            write (summary_output_unit, '(/,/," Interaction constants:"/,/," Level1  J Parity  Level2  J Parity",8x,"A (MHz)",13x,&
                &"B (MHz)"/)')
            write (hyperfine_output_unit, '(/,/," Interaction constants:"/,/," Level1  J Parity ",8x,"A (MHz)",13x,"B (MHz)",13x,&
                &"g_J",14x,&
                &"delta g_J",11x,"total g_J"/)')
            do first_output_level = 1, active_level_count
                do second_output_level = 1, active_level_count
                    first_level_two_j_plus_one = level_two_j_plus_one(first_output_level)
                    second_level_two_j_plus_one = level_two_j_plus_one(second_output_level)
                    if (.not.(first_level_two_j_plus_one==second_level_two_j_plus_one .and. second_level_two_j_plus_one>1 .or. &
                & first_level_two_j_plus_one>second_level_two_j_plus_one)) cycle
                    electronic_j = 0.5_real64*dble(first_level_two_j_plus_one - 1)
                    inverse_j_norm = sqrt(1.0_real64/(electronic_j*(electronic_j + 1.0_real64)))

                    dipole_same_j_scale = magnetic_dipole_mhz_scale*inverse_j_norm

                    if (first_level_two_j_plus_one == 2) then
                        dipole_delta_j_one_scale = 0.0_real64
                    else
                        dipole_delta_j_one_scale = magnetic_dipole_mhz_scale*sqrt(1.0_real64/(electronic_j*(2.0_real64* &
                & electronic_j - 1.0_real64)))
                    end if
                    quadrupole_same_j_scale = electric_quadrupole_mhz_scale*sqrt((electronic_j*(2.0_real64*electronic_j - &
                & 1.0_real64))/((electronic_j + 1.0_real64)*(2.0_real64*&
                        electronic_j + 3.0_real64)))
                    if (first_level_two_j_plus_one == 2) then
                        quadrupole_delta_j_one_scale = 0.0_real64
                    else
                        quadrupole_delta_j_one_scale = 0.25_real64*electric_quadrupole_mhz_scale*sqrt((electronic_j* &
                & (electronic_j - 1.0_real64))/((electronic_j + 1.0_real64)*(&
                            2.0_real64*electronic_j - 1.0_real64)))
                    end if
                    if (first_level_two_j_plus_one == 4) then
                        quadrupole_delta_j_two_scale = 0.0_real64
                    else
                        quadrupole_delta_j_two_scale = 0.125_real64*electric_quadrupole_mhz_scale*sqrt((electronic_j* &
                & (electronic_j - 1.0_real64)*(2.0_real64*electronic_j - 1.0_real64&
                            ))/(2.0_real64*electronic_j - 3.0_real64))
                    end if
                    ! Diagonal (J,J) A and B factors
                    select case (first_level_two_j_plus_one - second_level_two_j_plus_one)
                    case (0)
                        if (first_output_level <= second_output_level) then
                            ! GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II)
                            ! DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II)

                            write (summary_output_unit, '(1x,1i3,5x,2a,2x,1i3,5x,2a,1p,2d20.10)') active_level_indices( &
                & first_output_level), &
                                & format_j_label(first_level_two_j_plus_one), &
                                parity_labels((level_parity_signs(first_output_level)+3)/2), &
                                active_level_indices(second_output_level), format_j_label(second_level_two_j_plus_one), &
                                    parity_labels((level_parity_signs(second_output_level)+3)/2), dipole_same_j_scale* &
                & reduced_hyperfine_elements(1,&
                                active_level_count*(first_output_level-1)+second_output_level), quadrupole_same_j_scale* &
                & reduced_hyperfine_elements(3,active_level_count*(first_output_level-1)+second_output_level)
                            ! Output diagonal hfs and g_j factors to file <name>.hyperfine or <name>.cihyperfine
                            if (first_output_level == second_output_level) then
                                lande_g = vacuum_speed_of_light_au*inverse_j_norm*reduced_lande_g_elements(active_level_count* &
                & (first_output_level-1)+second_output_level)
                                lande_g_correction = 0.001160_real64*inverse_j_norm* &
                & reduced_lande_correction_elements(active_level_count*(first_output_level-1)+second_output_level)
                                write (hyperfine_output_unit, '(1x,1i3,5x,2a,1p,5d20.10)') active_level_indices( &
                & first_output_level), &
                                    format_j_label(first_level_two_j_plus_one), &
                                    & parity_labels((level_parity_signs(first_output_level)+3)/2), &
                                    dipole_same_j_scale*reduced_hyperfine_elements(1,active_level_count*(first_output_level-1)+ &
                & second_output_level), quadrupole_same_j_scale*reduced_hyperfine_elements(3,active_level_count* &
                & (first_output_level-1)+second_output_level), &
                                    lande_g, lande_g_correction, lande_g + lande_g_correction
                            end if
                        end if
                    ! Off diagonal (J,J-1) A and B factors
                    case (2)
                        write (summary_output_unit, '(1x,1i3,5x,2a,2x,1i3,5x,2a,1p,2d20.10)') &
                            active_level_indices(first_output_level), format_j_label(first_level_two_j_plus_one), &
                                parity_labels((level_parity_signs(first_output_level)+3)/2), active_level_indices(&
                            second_output_level), format_j_label(second_level_two_j_plus_one), &
                & parity_labels((level_parity_signs(second_output_level)+3)/2), dipole_delta_j_one_scale* &
                & reduced_hyperfine_elements(2, &
                                active_level_count*(first_output_level-1&
                            )+second_output_level), quadrupole_delta_j_one_scale*reduced_hyperfine_elements(4, &
                & active_level_count*(first_output_level-1)+second_output_level)
                    ! Off diagonal (J,J-2) B factor
                    case (4)
                        write (summary_output_unit, '(1x,1i3,5x,2a,2x,1i3,5x,2a,1p,2d20.10)') &
                            active_level_indices(first_output_level), format_j_label(first_level_two_j_plus_one), &
                                parity_labels((level_parity_signs(first_output_level)+3)/2), active_level_indices(&
                            second_output_level), format_j_label(second_level_two_j_plus_one), &
                & parity_labels((level_parity_signs(second_output_level)+3)/2), 0.0_real64, &
                                quadrupole_delta_j_two_scale*reduced_hyperfine_elements(5,&
                            active_level_count*(first_output_level-1)+second_output_level)
                    end select
                end do
            end do
            ! II = 2*I
            two_nuclear_spin = nint(2.0_real64*nuclear_spin)
            ! Calculate the F-dependent matrix elemnts
            !
            ! Loop over the states
            do first_recoupled_level = 1, active_level_count
                do second_recoupled_level = 1, active_level_count
                    first_level_two_j = level_two_j_plus_one(first_recoupled_level) - 1
                    second_level_two_j = level_two_j_plus_one(second_recoupled_level) - 1
                    if (.not.(second_level_two_j==first_level_two_j .and. second_level_two_j>0 .or. first_level_two_j> &
                & second_level_two_j)) cycle
                    ! Determine the possible F quantum numbers for the matrix element
                    lowest_two_f = max(abs(second_level_two_j - two_nuclear_spin),abs(first_level_two_j - two_nuclear_spin))
                    highest_two_f = min(second_level_two_j + two_nuclear_spin,first_level_two_j + two_nuclear_spin)
                    ! Loop over the possible F quantum numbers
                    wrote_matrix_heading = 0
                    do two_f = lowest_two_f, highest_two_f, 2
                        ! Phase factor
                        if (mod((two_nuclear_spin + second_level_two_j - two_f)/2,2) == 1) then
                            recoupling_phase = -1.0_real64
                        else
                            recoupling_phase = 1.0_real64
                        end if
                        recoupling_prefactor = recoupling_phase*sqrt((dble(first_level_two_j) + 1.0_real64)* &
                & (dble(two_nuclear_spin) + 1.0_real64))
                        !    Determine the Racah W coefficients.
                        call scott_fischer_racah_coefficient (two_nuclear_spin,second_level_two_j,two_nuclear_spin, &
                & first_level_two_j,two_f,2,dipole_racah_coefficient)
                        call scott_fischer_racah_coefficient (two_nuclear_spin,second_level_two_j,two_nuclear_spin, &
                & first_level_two_j,two_f,4,quadrupole_racah_coefficient)
                        ! Obtain and output matrix elements for J,J
                        if (second_level_two_j - first_level_two_j == 0) then
                            dipole_hyperfine_element = recoupling_prefactor*dipole_racah_coefficient* &
                & reduced_hyperfine_elements(1,active_level_count*(first_recoupled_level-1)+second_recoupled_level)* &
                & reduced_nuclear_dipole
                            quadrupole_hyperfine_element = recoupling_prefactor*quadrupole_racah_coefficient* &
                & reduced_hyperfine_elements(3,active_level_count*(first_recoupled_level-1)+second_recoupled_level)* &
                & reduced_nuclear_quadrupole
                            if (abs(dipole_hyperfine_element + quadrupole_hyperfine_element) > &
                                & contribution_magnitude_cutoff*10.0e-5_real64) then
                                if (wrote_matrix_heading == 0) write (summary_output_unit, '(/,/," Matrix elements:"/,/,&
                                    &" Level1  J Parity  Level2  J Parity    F ",4x,"Matrix element (a.u.)"/)')
                                wrote_matrix_heading = 1
                                write (summary_output_unit, '(1x,1i3,5x,2a,2x,1i3,5x,2a,1x,a,4x,1p,1d20.10)') &
                                    active_level_indices(first_recoupled_level), format_j_label(first_level_two_j+1), &
                                        parity_labels((level_parity_signs(first_recoupled_level)+3)&
                                    /2), active_level_indices(second_recoupled_level), &
                & format_j_label(second_level_two_j+1), &
                                        parity_labels((level_parity_signs(second_recoupled_level)+3)/2), &
                                    format_j_label(two_f+1), dipole_hyperfine_element + quadrupole_hyperfine_element
                            end if
                        ! Obtain and output matrix elements for J,J-1
                        else if (abs(second_level_two_j - first_level_two_j) == 2) then
                            dipole_hyperfine_element = recoupling_prefactor*dipole_racah_coefficient* &
                & reduced_hyperfine_elements(2,active_level_count*(first_recoupled_level-1)+second_recoupled_level)* &
                & reduced_nuclear_dipole
                            quadrupole_hyperfine_element = recoupling_prefactor*quadrupole_racah_coefficient* &
                & reduced_hyperfine_elements(4,active_level_count*(first_recoupled_level-1)+second_recoupled_level)* &
                & reduced_nuclear_quadrupole
                            if (abs(dipole_hyperfine_element + quadrupole_hyperfine_element) > &
                                & contribution_magnitude_cutoff*10.0e-5_real64) then
                                if (wrote_matrix_heading == 0) write (summary_output_unit, '(/,/," Matrix elements:"/,/,&
                                    &" Level1  J Parity  Level2  J Parity    F ",4x,"Matrix element (a.u.)"/)')
                                wrote_matrix_heading = 1
                                write (summary_output_unit, '(1x,1i3,5x,2a,2x,1i3,5x,2a,1x,a,4x,1p,1d20.10)') &
                                    active_level_indices(first_recoupled_level), format_j_label(first_level_two_j+1), &
                                        parity_labels((level_parity_signs(first_recoupled_level)+3)&
                                    /2), active_level_indices(second_recoupled_level), &
                & format_j_label(second_level_two_j+1), &
                                        parity_labels((level_parity_signs(second_recoupled_level)+3)/2), &
                                    format_j_label(two_f+1), dipole_hyperfine_element + quadrupole_hyperfine_element
                            end if
                        ! Obtain and output matrix elements for J,J-2
                        else if (abs(second_level_two_j - first_level_two_j) == 4) then
                            dipole_hyperfine_element = 0.0_real64
                            quadrupole_hyperfine_element = recoupling_prefactor*quadrupole_racah_coefficient* &
                & reduced_hyperfine_elements(5,active_level_count*(first_recoupled_level-1)+second_recoupled_level)* &
                & reduced_nuclear_quadrupole
                            if (abs(dipole_hyperfine_element + quadrupole_hyperfine_element) > &
                                & contribution_magnitude_cutoff*10.0e-5_real64) then
                                if (wrote_matrix_heading == 0) write (summary_output_unit, '(/,/," Matrix elements:"/,/,&
                                    &" Level1  J Parity  Level2  J Parity    F ",4x,"Matrix element (a.u.)"/)')
                                wrote_matrix_heading = 1
                                write (summary_output_unit, '(1x,1i3,5x,2a,2x,1i3,5x,2a,1x,a,4x,1p,1d20.10)') &
                                    active_level_indices(first_recoupled_level), format_j_label(first_level_two_j+1), &
                                        parity_labels((level_parity_signs(first_recoupled_level)+3)&
                                    /2), active_level_indices(second_recoupled_level), &
                & format_j_label(second_level_two_j+1), &
                                        parity_labels((level_parity_signs(second_recoupled_level)+3)/2), &
                                    format_j_label(two_f+1), dipole_hyperfine_element + quadrupole_hyperfine_element
                            end if
                        else
                            dipole_hyperfine_element = 0.0_real64
                            quadrupole_hyperfine_element = 0.0_real64
                        end if
                    end do
                end do
            end do
        end if
        call dalloc (reduced_hyperfine_elements, "reduced_hyperfine_elements", "calculate_hyperfine_matrix_elements_mpi")
        call dalloc (reduced_lande_g_elements, "reduced_lande_g_elements", "calculate_hyperfine_matrix_elements_mpi")
        call dalloc (reduced_lande_correction_elements, "reduced_lande_correction_elements", &
            "calculate_hyperfine_matrix_elements_mpi")
    end subroutine calculate_hyperfine_matrix_elements_mpi
end module calculate_hyperfine_matrix_elements_mpi_interface
