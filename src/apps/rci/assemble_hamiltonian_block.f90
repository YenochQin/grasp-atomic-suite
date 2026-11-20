module assemble_hamiltonian_block_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_hamiltonian_block
contains

    !***********************************************************************
    !                                                                      *
    subroutine assemble_hamiltonian_block (mpi_rank, mpi_process_count, symmetry_block_index, &
        & closed_shell_breit_shift,first_local_csf_column, initial_stored_element_count, &
        inverse_nuclear_mass,slf_en)
        ! Module dependencies
        !                                                                      *
        !   Sets up the Hamiltonian matrix and determines the average energy.  *
        !
        !   Serial I/O moved out; able to run on single/multiple processors
        !   For this purpose a new common /setham_to_genmat2/ is created
        !                                                                      *
        !   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
        !               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
        !                                                                      *
        !   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
        !   Block version by Xinghong He          Last revision: 15 Jun 1998   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use grasp_parameters, only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc
        use breit_core_occupations, only: breit_core_orbital_occupations
        use breit_integral_tables, only: breit_integral_counts
        use angular_coefficient_buffer_state, only: angular_coefficient_count, angular_coefficient_labels, angular_coefficients
        use coulomb_integral_table, only: coulomb_integral_count
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation, include_normal_mass_shift, &
            estimate_self_energy, include_specific_mass_shift, include_transverse_interaction, include_vacuum_polarization
        use debug_output_settings, only: cutoff, angular_coefficient_debug_flag
        use atomic_calculation_state, only: nuclear_mass_in_electron_masses
        use eigenstate_expansions, only: hamiltonian_energy_shift
        use csf_block_cutoffs, only: active_csf_cutoffs
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, matrix_partner_indices, &
            block_breit_integral_counts, block_stored_element_count, block_coulomb_one_body_integral_count, &
            block_coulomb_one_body_coefficient_count, block_slater_integral_count, block_coulomb_two_body_coefficient_count, &
            block_breit_coefficient_count, block_breit_core_coefficient_count, block_vacuum_polarization_integral_count, &
            block_kinetic_integral_count, block_vinti_integral_count, block_assembly_csf_count, block_matrix_element_cutoff
        use kinetic_energy_integral_table, only: kinetic_integral_count
        use orbital_csf_basis, only: csf_count, orbital_count
        use selected_level_indices, only: active_level_count
        use vinti_integral_table, only: vinti_integral_count
        use vacuum_polarization_integral_table, only: vacuum_polarization_integral_count
        use radial_wavefunctions, only: max_orbitals
        use distributed_slater_integrals, only: local_slater_integral_values, rk_local_count
        use rci_prerun_settings, only: prerun_one_body_weight_cutoff, prerun_two_body_weight_cutoff, prerun_stage, prerun_csf_count
        use ci_restart_units, only: restart_matrix_unit
        use scf_previous_eigenvectors, only: previous_csf_mixing_coefficients
        use angular_coefficient_buffer, only: manage_angular_coefficient_buffer

        ! Procedure interfaces
        use one_electron_integral_interface, only: one_electron_integral
        use breit_radial_integral_1_interface, only: breit_radial_integral_1
        use breit_radial_integral_2_interface, only: breit_radial_integral_2
        use breit_radial_integral_3_interface, only: breit_radial_integral_3
        use breit_radial_integral_4_interface, only: breit_radial_integral_4
        use breit_radial_integral_5_interface, only: breit_radial_integral_5
        use breit_radial_integral_6_interface, only: breit_radial_integral_6

        ! Module dependencies
        use integer_format, only: format_integer
        use subshell_occupation_class, only: classify_subshell_occupation

        ! Procedure interfaces
        use kinetic_energy_integral_interface, only: kinetic_energy_integral
        use lookup_slater_integral_interface, only: lookup_slater_integral
        use calculate_vinti_integral_interface, only: calculate_vinti_integral
        use calculate_vacuum_polarization_integral_interface, only: calculate_vacuum_polarization_integral

        ! Module dependencies
        use one_electron_scalar_coefficients, only: calculate_one_electron_scalar_coefficients
        use two_electron_angular_coefficients, only: calculate_two_electron_angular_coefficients
        use closed_shell_mcp, only: calculate_closed_shell_mcp
        implicit none (type, external)

        ! Arguments
        ! Arguments and local declarations
        integer(int64)       :: initial_stored_element_count
        integer             :: symmetry_block_index, first_local_csf_column
        integer, intent(in) :: mpi_rank, mpi_process_count
        real(real64)        :: closed_shell_breit_shift, inverse_nuclear_mass
        real(real64) :: slf_en(*)

        ! Local declarations
        real(real64) :: shell_coefficients(max_orbitals)
        real(real64) :: tgrl1, tgrl2, radial_integral_value

        ! Local constants

        ! Local declarations
        !
        !     Matrix elements smaller than CUTOFF are not accumulated
        !
        !cjb cutoff is use associated and cannot be redeclared
        !     REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
        !cjb CUTOFF = 1.0D-12 below
        !cjb
        integer :: ipi, csf_index, inc1, inc2, kt, ipt, include_closed_shell_terms, ncoec, nctec, &
            i, j, nmcbp, closed_shell_contribution_count, ic, nelc, irstart, ket_csf_index, orbital_index, ket_orbital_index,  &
            itype, iia
        real(real64) :: elemnt, precoeff, tcoeff, vcoeff, contr
        !-----------------------------------------------------------------------
        print *, "Calling setham ..."
        !cjb
        !cjb  CUTOFF = 1.0D-12
        cutoff = 1.0e-12_real64

        stored_matrix_element_count = initial_stored_element_count

        !Per Bug fix 30/4-2013
        inverse_nuclear_mass = 1.0_real64/nuclear_mass_in_electron_masses
        !Per Bug fix 30/4-2013

        if (prerun_stage == 2) then
            do ipi = 1,active_level_count
                do csf_index = 1,csf_count
                    write (*,*) ipi,csf_index,previous_csf_mixing_coefficients(csf_index+(ipi-1)*csf_count)
                end do
            end do
        end if
        !
        !   Allocate storage to arrays in COMMON/BUFFER/; these are
        !   used for the Coulomb and transverse two-electron integrals
        !
        call manage_angular_coefficient_buffer (1)

        !     ...Locals
        call alloc (hamiltonian_values, csf_count, "EMT","SETHAM" )
        call alloc (matrix_partner_indices, csf_count, "IROW", "SETHAM")
        !
        inc1 = 1
        inc2 = 1
        !
        !   Initialisations for contributions from the Dirac-Coulomb
        !   operator
        !
        kt  = 0
        ipt = 1
        !
        include_closed_shell_terms = 1

        ncoec = 0
        !
        nctec   = 0

        if (include_transverse_interaction) then

            !        ...Initialisations for transverse interaction correction
            do i = 1, orbital_count
                breit_core_orbital_occupations(i) = 0
                do j = 1, csf_count
                    if (classify_subshell_occupation(i,j) <= 0) exit
                end do
                if (j > csf_count) breit_core_orbital_occupations(i) = 1
            end do

            nmcbp = 0
            closed_shell_contribution_count = 0
        end if

        ! Loop over rows of the Hamiltonian matrix - distributed

        do ic = first_local_csf_column, csf_count, mpi_process_count

            nelc = 0    ! counter - Number of non-zeros of this row

            ! Loop over columns of the current row

            irstart = 1
            do ket_csf_index = irstart, ic

                ! PER
                if (use_first_order_csf_perturbation .and. (ket_csf_index > active_csf_cutoffs(1))) then
                    if (ket_csf_index/=ic) cycle
                end if
                ! PER

                elemnt = 0.0_real64     ! accumulates various contributions to H

                !
                !   Generate the integral list for the matrix element of the
                !   one-body operators
                !
                if (prerun_stage == 1) then
                    inc1 = 0
                    inc2 = 0
                    if (ic<=prerun_csf_count .or. ic==ket_csf_index) then
                        inc1 = 1
                    end if
                end if

                if (prerun_stage == 2) then
                    !
                    !   Diagonal elements are always included
                    !   Off diagonal elements are included only if the
                    !   products of the weights from the prerun are larger
                    !   than the cutoff.
                    !
                    if (ic == ket_csf_index) then
                        inc1 = 1
                        inc2 = 1
                    else
                        inc1 = 0
                        inc2 = 0
                    end if
                    do ipi = 1,active_level_count
                        precoeff =                                           &
                            dabs(previous_csf_mixing_coefficients(ic+(ipi-1)*csf_count)* &
                                previous_csf_mixing_coefficients(ket_csf_index+(ipi-1)*csf_count))
                        if (precoeff > prerun_one_body_weight_cutoff) inc1 = 1
                        if (precoeff > prerun_two_body_weight_cutoff) inc2 = 1
                    end do
                end if

                !            ...INC1.EQ.1 ------------>
                if (inc1 == 1) then   !inc1 is always 1 without PRE-RUN
                    call calculate_one_electron_scalar_coefficients(ic,ket_csf_index,orbital_index,ket_orbital_index, &
                        & shell_coefficients)
                    !
                    !   Accumulate the contribution from the one-body operators:
                    !   kinetic energy, electron-nucleus interaction; update the
                    !   angular integral counter
                    !
                    if (orbital_index /= 0) then
                        if (orbital_index == ket_orbital_index) then
                            do orbital_index = 1,orbital_count
                                iia = orbital_index
                                tcoeff = shell_coefficients(orbital_index)
                                if (abs (tcoeff) > cutoff) then
                                    ncoec = ncoec + 1
                                    call one_electron_integral (iia, iia, radial_integral_value)
                                    !------------------------
                                    elemnt = elemnt + radial_integral_value*tcoeff
                                    if (include_normal_mass_shift) then
                                        call kinetic_energy_integral (iia,iia,radial_integral_value)
                                        !------------------------
                                        elemnt = elemnt + radial_integral_value*inverse_nuclear_mass*tcoeff
                                    end if
                                    if (include_vacuum_polarization) then
                                        call calculate_vacuum_polarization_integral (iia, iia, radial_integral_value)
                                        !------------------------
                                        elemnt = elemnt + radial_integral_value*tcoeff
                                    end if
                                end if
                            end do
                        else
                            tcoeff = shell_coefficients(1)
                            if (abs (tcoeff) > cutoff) then
                                ncoec = ncoec + 1
                                call one_electron_integral (orbital_index, ket_orbital_index, radial_integral_value)
                                !------------------------
                                elemnt = elemnt + radial_integral_value*tcoeff
                                if (include_normal_mass_shift) then
                                    call kinetic_energy_integral (orbital_index, ket_orbital_index, radial_integral_value)
                                    !------------------------
                                    elemnt = elemnt + radial_integral_value*inverse_nuclear_mass*tcoeff
                                end if
                                if (include_vacuum_polarization) then
                                    call calculate_vacuum_polarization_integral (orbital_index, ket_orbital_index, &
                                & radial_integral_value)
                                    !------------------------
                                    elemnt = elemnt + radial_integral_value*tcoeff
                                end if
                            end if
                        end if
                    end if
                    !
                    angular_coefficient_debug_flag = 0
                    !
                    !   Accumulate the contributions from the two-electron
                    !   Coulomb operator and the mass polarisation; the latter
                    !   is computed first because the orbital indices may be
                    !   permuted by RKINTC
                    !
                    angular_coefficient_count = 0
                    !
                    call calculate_two_electron_angular_coefficients (ic, ket_csf_index, calculate_closed_shell_mcp, &
                        & include_closed_shell_terms, 1)
                    !
                    do i = 1, angular_coefficient_count
                        vcoeff = angular_coefficients(i)
                        if (abs (vcoeff) > cutoff) then
                            nctec = nctec + 1
                            if (include_specific_mass_shift) then
                                if (angular_coefficient_labels(5,i) == 1) then
                                    call calculate_vinti_integral (angular_coefficient_labels(1,i), angular_coefficient_labels(3, &
                                        & i), tgrl1)
                                    call calculate_vinti_integral (angular_coefficient_labels(2,i), angular_coefficient_labels(4, &
                                        & i), tgrl2)
                                    elemnt = elemnt - tgrl1*tgrl2*inverse_nuclear_mass*vcoeff
                                end if
                            end if
                            call lookup_slater_integral (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),                     &
                                angular_coefficient_labels(5,i), radial_integral_value)
                            elemnt = elemnt + radial_integral_value*vcoeff
                        end if
                    end do
                    !
                    angular_coefficient_debug_flag = 0

                end if  !inc1 is always 1 without PRE-RUN
                !            ...INC1.EQ.1 <------------
                !***********************************************************************
                !            ...LTRANS .AND. (INC2.EQ.1) ------------>
                if (include_transverse_interaction .and. (inc2==1)) then
                    !IF (INC2 .EQ. 1) THEN  !inc2 is always 1 without PRE-RUN
                    !
                    !   Accumulate the contribution from the two-electron
                    !   transverse interaction operator
                    !
                    angular_coefficient_count = 0
                    !
                    call calculate_two_electron_angular_coefficients (ic, ket_csf_index, closed_shell_breit_callback, 1, 2)
                    !
                    do i = 1, angular_coefficient_count
                        if (dabs (angular_coefficients(i)) > cutoff) then
                            nmcbp = nmcbp + 1
                            itype = abs (angular_coefficient_labels(6,i))
                            if (itype == 1) then
                                call breit_radial_integral_1 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                    angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                    angular_coefficient_labels(5,i), radial_integral_value)
                            else if (itype == 2) then
                                call breit_radial_integral_2 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                    angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                    angular_coefficient_labels(5,i), radial_integral_value)
                            else if (itype == 3) then
                                call breit_radial_integral_3 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                    angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                    angular_coefficient_labels(5,i), radial_integral_value)
                            else if (itype == 4) then
                                call breit_radial_integral_4 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                    angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                    angular_coefficient_labels(5,i), radial_integral_value)
                            else if (itype == 5) then
                                call breit_radial_integral_5 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                    angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                    angular_coefficient_labels(5,i), radial_integral_value)
                            else if (itype == 6) then
                                call breit_radial_integral_6 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i), &
                                    angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                    angular_coefficient_labels(5,i), radial_integral_value)
                            end if
                            contr = angular_coefficients(i)*radial_integral_value
                            if (angular_coefficient_labels(6,i) > 0) then
                                elemnt = elemnt + contr
                            else
                                !                        ...It comes here only when ic=ir=1
                                !                           A negative type label marks a closed-shell Breit contribution.
                                closed_shell_contribution_count = closed_shell_contribution_count + 1
                                closed_shell_breit_shift = closed_shell_breit_shift + contr
                            end if
                        end if
                    end do
                    !
                    angular_coefficient_debug_flag = 0
                    !
                    !               ...ELSTO is a constant over all diagonals, thus its
                    !                  contribution to the total energy can be added later
                    !            IF (IR .EQ. IC) ELEMNT = ELEMNT + ELSTO
                    !
                    !ENDIF   !inc2 is always 1 without PRE-RUN
                end if
                !            ...LTRANS .AND. (INC2.EQ.1) <------------
                !***********************************************************************
                !
                ! Store this element if it is diagonal or its value is greater than
                ! CUTOFF
                !
                if ( (ket_csf_index == ic) .or. (dabs (elemnt) > cutoff) ) then
                    nelc       = nelc + 1
                    hamiltonian_values(nelc)  = elemnt
                    matrix_partner_indices(nelc) = ket_csf_index
                end if
                !
            end do
            ! zou
            !        print *, ic,SLF_EN(IC)
            if(estimate_self_energy) hamiltonian_values(nelc) = hamiltonian_values(nelc) + slf_en(ic)
            ! zou
            !
            !   This column is done; write it to disk
            !
            write (restart_matrix_unit) nelc, closed_shell_breit_shift, (hamiltonian_values(ket_csf_index), ket_csf_index = &
                & 1, nelc),            &
                (matrix_partner_indices(ket_csf_index), ket_csf_index = 1, nelc)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This EAV (and the above EMT) does not have ELSTO.
            hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(nelc)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !-----------------------------------------------------------------------
            !         IF (MOD (IC, 20) .EQ. 0 .OR.                                   &
            if (mod (ic, 100) == 0 .or.                                   &
                ic < mpi_process_count*2 .or. ic > (csf_count-mpi_process_count*2)) then
                print *, "Row ", ic, ": ", nelc, " nonzero elements;"       &
                    , "  block = ", symmetry_block_index
            end if
            !
            !   Update the counter for the total number of elements
            !
            stored_matrix_element_count = stored_matrix_element_count + nelc
            !
        end do
        !***********************************************************************
        !
        !   Deallocate storage for the arrays in /BUFFER/
        !
        call manage_angular_coefficient_buffer (3)

        !     ...Locals
        call dalloc (hamiltonian_values, "EMT", "SETHAM")
        call dalloc (matrix_partner_indices, "IROW", "SETHAM")

        !  Fill the common block /setham_to_genmat2/ for use in genmat2

        block_matrix_element_cutoff = cutoff
        block_coulomb_one_body_integral_count = coulomb_integral_count
        block_coulomb_one_body_coefficient_count = ncoec
        block_slater_integral_count = rk_local_count(size(local_slater_integral_values), mpi_rank, mpi_process_count)
        block_coulomb_two_body_coefficient_count = nctec
        block_breit_integral_counts = breit_integral_counts
        block_breit_coefficient_count = nmcbp
        block_breit_core_coefficient_count = closed_shell_contribution_count
        block_vacuum_polarization_integral_count = vacuum_polarization_integral_count
        block_kinetic_integral_count = kinetic_integral_count
        block_vinti_integral_count = vinti_integral_count
        block_stored_element_count = stored_matrix_element_count
        block_assembly_csf_count = csf_count

        return
    end subroutine assemble_hamiltonian_block
    ! Match the legacy CORD callback characteristics at the angular-library boundary.
    subroutine closed_shell_breit_callback(bra_csf_index, ket_csf_index, bra_shell_position, orbital_list_selector, &
        ket_shell_position)
        ! Procedure interfaces
        use accumulate_closed_shell_breit_interface, only: accumulate_closed_shell_breit
        implicit none (type, external)

        ! Arguments
        integer :: bra_csf_index, ket_csf_index
        integer, intent(in) :: bra_shell_position, orbital_list_selector, ket_shell_position

        call accumulate_closed_shell_breit(bra_csf_index, ket_csf_index, bra_shell_position, orbital_list_selector, &
            ket_shell_position)
    end subroutine closed_shell_breit_callback
end module assemble_hamiltonian_block_interface
