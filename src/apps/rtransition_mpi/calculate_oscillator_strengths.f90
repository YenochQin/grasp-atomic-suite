module calculate_oscillator_strengths_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_oscillator_strengths
contains

    !***********************************************************************
    !                                                                      *
    subroutine calculate_oscillator_strengths(name,fullname,tmpdir,startup_directory,block_label_text)
        ! Module dependencies
        use transition_matrix_workspace, only: final_mixing_unit, initial_mixing_unit
        use calculation_file_units, only: lsj_transition_output_unit, summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real32
        !   This routine controls the main sequence of routine calls for the   *
        !   calculation  of  data for transitions between multiconfiguration   *
        !   Dirac-Fock energy levels.                                          *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc
        use transition_basis_blocks, only: initial_transition_block_count, final_transition_block_count
        use atomic_calculation_state, only: speed_of_light_cm_per_s, electron_count, hartree_to_ev, hartree_to_wavenumber_cm, &
            speed_of_light_au
        use calculation_input_defaults, only: input_mixing_file_kind
        use eigenstate_expansions, only: final_hamiltonian_energy_shift, initial_hamiltonian_energy_shift, &
            final_shifted_eigenvalues, initial_shifted_eigenvalues, final_csf_mixing_coefficients, initial_csf_mixing_coefficients
        use transition_matrix_workspace, only: transition_matrix_capacity, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, &
            transition_multipole_rank, transition_multipole_type, transition_coefficient_capacity, &
            transition_angular_coefficients, total_coulomb_probabilities, total_babushkin_probabilities, &
            transition_bra_csf_indices, transition_ket_csf_indices, transition_integral_coefficient_offsets, &
            packed_transition_integral_labels, requested_multipole_count, requested_multipole_codes, transition_option_flags
        use selected_level_indices, only: final_level_indices, initial_level_indices, final_level_count, initial_level_count
        use level_symmetry_metadata, only: final_level_parity_signs, initial_level_parity_signs, final_level_two_j_plus_one, &
            initial_level_two_j_plus_one
        use transition_lsj_labels, only: nvectoti, nvectotf, iopen_status1, iopen_status2
        use mpi_state, only: mpi_allreduce, mpi_barrier, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_rank, mpi_sum

        ! Procedure interfaces
        use prepare_transition_mixing_coefficients_interface, only: prepare_transition_mixing_coefficients
        use resize_transition_matrix_arrays_interface, only: resize_transition_matrix_arrays
        use resize_transition_integral_arrays_interface, only: resize_transition_integral_arrays
        use map_merged_orbitals_interface, only: map_merged_orbitals
        use write_transition_angular_coefficients_interface, only: write_transition_angular_coefficients
        use read_transition_mixing_coefficients_interface, only: read_transition_mixing_coefficients
        use read_transition_angular_coefficients_interface, only: read_transition_angular_coefficients

        ! Module dependencies
        use angular_momentum_triangle, only: satisfies_angular_triangle

        ! Procedure interfaces
        use build_transition_bessel_functions_interface, only: build_transition_bessel_functions
        use write_transition_csf_weights_interface, only: write_transition_csf_weights
        use write_transition_probability_interface, only: write_transition_probability
        use write_lsj_transition_probability_interface, only: write_lsj_transition_probability

        ! Module dependencies
        use posix_services, only: terminate_process

        ! Procedure interfaces
        use load_initial_lsj_labels_interface, only: load_initial_lsj_labels
        use load_final_lsj_labels_interface, only: load_final_lsj_labels
        use write_input_summary_interface, only: write_input_summary

        ! Module dependencies
        use mpi_change_directory, only: change_mpi_directory
        implicit none (type, external)

        ! Arguments
        character (len = 24), intent(inout) :: name(2)
        character(len=3)  ::  block_label_text
        character (len = 128), intent(inout) :: fullname(2)
        character(len=128) :: startup_directory,tmpdir

        ! Local constants
        character(len=4), parameter :: hartree_unit_label = "Hart"
        character(len=4), parameter :: electronvolt_unit_label = " eV "
        character(len=4), parameter :: wavenumber_unit_label = "Kays"
        character(len=4), parameter :: hertz_unit_label = " Hz "
        character(len=4), parameter :: angstrom_unit_label = " A  "
        integer, parameter :: nca = 65536

        ! Local declarations
        integer :: copied_mixing_unit
        integer :: kka, transition_type_index, transition_angular_unit, io_status, my_nelec, initial_csf_count, &
            initial_orbital_count, initial_mixing_storage_size, iopar, initial_block_index, final_csf_count, &
            final_orbital_count, final_mixing_storage_size, final_block_index, level_pair_count, m, ielec, initial_level_index, &
            final_level_index, itkpo, itest, nlp, lines, ilbl, inum_ii, inum_ff, icount1, icount2, iii
        real(real64) :: factor, transition_angular_frequency, bessel_argument, coulomb_or_magnetic_amplitude, &
            & babushkin_amplitude, asfatmp, asfbtmp
        real(real64), pointer :: et(:), et1(:), ipr(:), ipr1(:), next(:)
        logical :: same_initial_and_final_state
        character :: iunits*4, mixing_file_signature*6
        !-----------------------------------------------
        !
        ! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
        ! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
        !
        ! for the case that the inital file and final file are same
        !
        same_initial_and_final_state = trim(fullname(1)) == trim(fullname(2))
        if (same_initial_and_final_state.and.mpi_rank==0) call prepare_transition_mixing_coefficients (fullname, &
            input_mixing_file_kind)
        call mpi_barrier(mpi_comm_world,mpi_error_code)
        !
        ! write header for the result file
        !
        !cjb fort.24
        if (mpi_rank == 0) then
            if (same_initial_and_final_state) then
                write(summary_output_unit,*) "Transition in file:"
                write(summary_output_unit,*) "f = ",trim(name(1))
            else
                write(summary_output_unit,*) "Transition between files:"
                write(summary_output_unit,*) "f1 = ",trim(name(1))
                write(summary_output_unit,*) "f2 = ",trim(name(2))
            end if
        end if ! myid .EQ. 0
        !
        call resize_transition_matrix_arrays (jja, jjb, hb1, hb2, hc1, hc2, hm1, &
            hm2, packed_transition_integral_labels, transition_integral_coefficient_offsets, transition_matrix_capacity, 1)
        call resize_transition_integral_arrays (transition_bra_csf_indices, transition_ket_csf_indices, &
            transition_angular_coefficients, transition_coefficient_capacity, 1)
        !
        !   Make a connection between the orbitals of the
        !   merged list and the initial and final state lists
        !
        call map_merged_orbitals
        !
        !   Set up units for printing transition energy
        !
        if (transition_option_flags(1)) then
            !
            !   Print transition energies in Angstroms
            !
            factor = hartree_to_wavenumber_cm
            factor = 1.0e8_real64/factor
            iunits = angstrom_unit_label
        !
        else if (transition_option_flags(2)) then
            !
            !   Print energies in eV
            !
            factor = hartree_to_ev
            iunits = electronvolt_unit_label
        !
        else if (transition_option_flags(3)) then
            !
            !   Print transition energies in Hartree Atomic Units
            !
            factor = 1.0_real64
            iunits = hartree_unit_label
        !
        else if (transition_option_flags(4)) then
            !
            !   Print transition energies in Hz
            !
            factor = hartree_to_wavenumber_cm
            factor = factor*speed_of_light_cm_per_s
            iunits = hertz_unit_label
        !
        else if (transition_option_flags(5)) then
            !
            !   Print transition energies in Kaysers
            !
            factor = hartree_to_wavenumber_cm
            iunits = wavenumber_unit_label
            !
        end if
        !
        !   Select type of transition
        !
        !   KK  =  0  for electric multipole.
        !       =  1  for magnetic multipole.
        !
        !   CF:  IOPAR  =  (-1)**N      Electric N-pole
        !               =  (-1)**(N+1)  Magnetic N-pole.
        !                               N > 0
        !
        kka = 1
        do transition_type_index = 1, requested_multipole_count
            !
            ! read the head of the file of mixing coef. for initial
            !
            if (same_initial_and_final_state) fullname(1) = trim(fullname(2))//"_CP"
            if (input_mixing_file_kind == 0) then
                open(newunit=initial_mixing_unit,file=trim(fullname(1))//".cbm",              &
                    form="UNFORMATTED",status="OLD")
            else
                open(newunit=initial_mixing_unit, file=trim(fullname(1))//".bm",              &
                    form="UNFORMATTED",status="OLD")
            end if
            if (same_initial_and_final_state) fullname(1) = fullname(2)
            read (initial_mixing_unit, iostat=io_status) mixing_file_signature
            if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
                write (*, *) "Not a GRASP mixing file"
                stop
            end if
            read (initial_mixing_unit) electron_count, initial_csf_count, initial_orbital_count, nvectoti, &
                & initial_mixing_storage_size, initial_transition_block_count
            write (*, *) "   nelec  = ", my_nelec
            write (*, *) "   ncftoti = ", initial_csf_count
            write (*, *) "   nwi     = ", initial_orbital_count
            write (*, *) "   nblocki = ", initial_transition_block_count
            write (*, *)
            !   lbl
            icount1 = 0
            call load_initial_lsj_labels (fullname(1))
            ! If not available generate angular coefficients for all pares of blocks
            call write_transition_angular_coefficients (iopar, transition_type_index, name, transition_angular_unit)
            do initial_block_index = 1, initial_transition_block_count
                call read_transition_mixing_coefficients (name, input_mixing_file_kind, 1)
                !   lbl
                icount1 = initial_level_count + icount1
                !
                ! read the head of the file of mixing coef. for final
                !
                if (input_mixing_file_kind == 0) then
                    open(newunit=final_mixing_unit,file=trim(fullname(2))//".cbm",          &
                        form="UNFORMATTED",status="OLD")
                else
                    open(newunit=final_mixing_unit, file=trim(fullname(2))//".bm",          &
                        form="UNFORMATTED",status="OLD")
                end if
                read (final_mixing_unit, iostat=io_status) mixing_file_signature
                if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
                    write (*, *) "Not a GRASP mixing file"
                    stop
                end if
                read (final_mixing_unit) electron_count, final_csf_count, final_orbital_count, nvectotf, &
                    & final_mixing_storage_size, final_transition_block_count
                write (*, *) "   nelec  = ", my_nelec
                write (*, *) "   ncftotf = ", final_csf_count
                write (*, *) "   nwf     = ", final_orbital_count
                write (*, *) "   nblockf = ", final_transition_block_count
                write (*, *)
                !GG  lbl
                if(initial_block_index == 1) then
                    call load_final_lsj_labels (fullname(2))
                    if(mpi_rank==0.and.iopen_status1==0.and.iopen_status2==0) then
                        iii = len_trim(startup_directory)
                        call change_mpi_directory(trim(startup_directory),iii,mpi_error_code)
                        if (mpi_error_code/=0) then
                            print *, "error changing dir!"
                            call terminate_process(1)
                        end if
                        ilbl = 1
                        call write_input_summary(name,input_mixing_file_kind,ilbl)
                        iii = len_trim(tmpdir)
                        call change_mpi_directory(trim(tmpdir)//"/"//block_label_text,iii+4,mpi_error_code)
                        if (mpi_error_code/=0) then
                            print *, "Error changing dir!"
                            call terminate_process(1)
                        end if
                        write(lsj_transition_output_unit,*) "Transition between files:"
                        write(lsj_transition_output_unit,*) trim(name(1))
                        write(lsj_transition_output_unit,*) trim(name(2))
                    end if
                end if
                icount2 = 0
                !GG   lbl end
                do final_block_index = 1, final_transition_block_count
                    call read_transition_mixing_coefficients (name, input_mixing_file_kind, 2)
                    !   lbl
                    icount2 = final_level_count + icount2
                    !
                    !   Allocate storage
                    !
                    call alloc (total_babushkin_probabilities, final_level_count, "TOTB", "OSCL")
                    call alloc (total_coulomb_probabilities, final_level_count, "TOTC", "OSCL")
                    !
                    level_pair_count = initial_level_count*final_level_count
                    call alloc (et, level_pair_count, "ET", "OSCL" )
                    call alloc (et1, level_pair_count, "ET1","OSCL" )
                    call alloc (ipr, level_pair_count, "IPR", "OSCL" )
                    call alloc (ipr1, level_pair_count,"IPR1", "OSCL" )
                    call alloc (next, level_pair_count, "NEXT", "OSCL")
                    !
                    !   Initialization for total decay rate
                    !
                    total_coulomb_probabilities(:final_level_count) = 0.0_real64
                    total_babushkin_probabilities(:final_level_count) = 0.0_real64
                    !
                    call read_transition_angular_coefficients (iopar, transition_type_index, name, transition_angular_unit)
                    !
                    if (transition_multipole_rank > 0) then
                        ielec = (-1)**transition_multipole_rank
                        if (ielec == iopar) then
                            transition_multipole_type = 0
                            kka = 0
                        else
                            transition_multipole_type = 1
                            ielec = -ielec
                        end if
                        !
                        !   Set up list of levels for calculation of oscillator strengths
                        !   sort list into increasing order of energy if option 6 set
                        !
                        if (mpi_rank==0.and.initial_block_index==1.and.final_block_index==1) then
                            if (transition_multipole_type == 0) then
                                write (summary_output_unit, '(/,/," Electric 2**(",i2,")-pole transitions")') &
                                    transition_multipole_rank
                            else
                                write (summary_output_unit, '(/,/," Magnetic 2**(",i2,")-pole transitions")') &
                                    transition_multipole_rank
                            end if
                            write (summary_output_unit, '(1x,33("="))')
                            if (transition_option_flags(1)) then
                                write (summary_output_unit, &
                                    & '(/,"   Upper state        Lower state  ",8x,"Gauge",8x,"Wavelength",13x,&
                                    &"Einstein coefficients",13x,"Oscillator")')
                                if (.not.transition_option_flags(7)) then
                                    write (summary_output_unit, &
                                        & '(81x,"-1",15x,"3 -2 -1",/," Level  J Parity",4x,"Level  J ","Parity",21x,&
                                        &"(Angstroms)",10x,"A (s  )",9x,"gB (m s  J  )",7x,"strength gf"/)')
                                else
                                    write (summary_output_unit, '(" Level  J Parity",4x,"Level  J Parity",21x,"(Angstroms)",10x,&
                                        &"A (au)",13x,&
                                        &"gB (au)",10x,"strength gf"/)')
                                end if
                            else
                                write (summary_output_unit, '(/," Upper       Lower ")')
                                if (.not.transition_option_flags(7)) then
                                    write (summary_output_unit, &
                                        '(" Lev  J P",3x,"Lev  J P",7x,"E (",a4,")",9x,"A (s-1)",10x,"gf",12x,"S")') &
                                        iunits
                                else
                                    write (summary_output_unit, &
                                        & '(" Level  J Parity",4x,"Level  J Parity",23x,"(",a4,")",13x,"A (au)",13x,&
                                        &"gB (au)",10x,"strength gf"/)') iunits
                                end if
                            end if
                        end if !myid=0
                        !
                        do initial_level_index = 1, initial_level_count
                            !   lbl
                            inum_ii = icount1 - initial_level_count + initial_level_index
                            do final_level_index = 1, final_level_count
                                !   lbl
                                inum_ff = icount2 - final_level_count + final_level_index
                                !
                                !   Check for consistent parity and J
                                !
                                itkpo = transition_multipole_rank + transition_multipole_rank + 1
                                if (satisfies_angular_triangle(initial_level_two_j_plus_one(initial_level_index), &
                                    final_level_two_j_plus_one(final_level_index), itkpo) == 0) &
                                    cycle
                                itest = initial_level_parity_signs(initial_level_index)* &
                                    final_level_parity_signs(final_level_index)*ielec
                                if (itest < 0) cycle
                                !
                                !   Calculate and print transition probability data
                                !
                                nlp = 70 - 8
                                lines = nlp
                                !
                                if (lines >= nlp) lines = 0
                                !
                                m = final_level_index + final_level_count*(initial_level_index - 1)
                                !             M = LEVFF+NVECII*(LEVII-1)
                                et(m) = final_shifted_eigenvalues(final_level_index) + final_hamiltonian_energy_shift - &
                                    initial_shifted_eigenvalues(initial_level_index) - initial_hamiltonian_energy_shift
                                if (same_initial_and_final_state .and. et(m)<=0.0_real32) cycle
                                transition_angular_frequency = -et(m)
                                bessel_argument = transition_angular_frequency/speed_of_light_au
                                call build_transition_bessel_functions (bessel_argument)
                                !
                                !  Calculate oscillator strength between the ASFs
                                !
                                call write_transition_csf_weights (asfatmp,asfbtmp,initial_level_index,final_level_index)
                                call mpi_allreduce(asfatmp,coulomb_or_magnetic_amplitude,1,mpi_double_precision, &
                                    mpi_sum, mpi_comm_world, mpi_error_code)
                                call mpi_allreduce(asfbtmp,babushkin_amplitude,1,mpi_double_precision, &
                                    mpi_sum, mpi_comm_world, mpi_error_code)
                                if (mpi_rank == 0) then
                                    call write_transition_probability(coulomb_or_magnetic_amplitude, babushkin_amplitude, &
                                        initial_level_index, final_level_index, transition_angular_frequency, factor, &
                                        lines,same_initial_and_final_state)
                                    if(iopen_status1==0 .and. iopen_status2 ==0) then
                                        call write_lsj_transition_probability &
                                & (inum_ii, inum_ff, coulomb_or_magnetic_amplitude, babushkin_amplitude, &
                                    initial_level_index, final_level_index, &
                                            transition_angular_frequency,factor)
                                    end if
                                end if
                                !              WRITE (24,317)
                            end do
                        end do
                    end if
                    !
                    !   Deallocate storage; this is local to OSCL
                    !
                    call dalloc (total_babushkin_probabilities, "TOTB", "OSCL")
                    call dalloc (total_coulomb_probabilities, "TOTC", "OSCL")
                    !
                    call dalloc (et, "ET", "OSCL")
                    call dalloc (et1, "ET!", "OSCL")
                    call dalloc (ipr, "IPR", "OSCL")
                    call dalloc (ipr1, "IPR1", "OSCL")
                    call dalloc (next, "NEXT", "OSCL")
                    !
                    !   Deallocate storage; this is allocated in READMIX
                    !
                    call dalloc (final_shifted_eigenvalues, "EVALFF", "OSCL")
                    call dalloc (final_csf_mixing_coefficients, "EVECFF", "OSCL")
                    call dalloc (final_level_indices, "IVECFF", "OSCL")
                    call dalloc (final_level_two_j_plus_one, "IATJPOFF", "OSCL")
                    call dalloc (final_level_parity_signs, "IASPARFF", "OSCL")

                end do
                close(final_mixing_unit)
                !
                !   Deallocate storage; this is allocated in READMIX
                !
                call dalloc (initial_shifted_eigenvalues, "EVALII", "OSCL")
                call dalloc (initial_csf_mixing_coefficients, "EVECII", "OSCL")
                call dalloc (initial_level_indices, "IVECII", "OSCL")
                call dalloc (initial_level_two_j_plus_one, "IATJPOII", "OSCL")
                call dalloc (initial_level_parity_signs, "IASPARII", "OSCL")

            end do
            close(initial_mixing_unit)
            close(transition_angular_unit)
        end do
        call dalloc (requested_multipole_codes, "KP", "OSCL")
        !
        ! close and delete duplicated mixing file
        if (same_initial_and_final_state) then
            if (input_mixing_file_kind == 0) then
                open(newunit=copied_mixing_unit, file=trim(fullname(2))//"_CP.cbm")
            else
                open(newunit=copied_mixing_unit, file=trim(fullname(2))//"_CP.bm")
            end if
            !GG         CLOSE(68, STATUS='delete')
        end if
        !
        call resize_transition_matrix_arrays (jja, jjb, hb1, hb2, hc1, hc2, hm1, hm2, packed_transition_integral_labels, &
            transition_integral_coefficient_offsets, transition_matrix_capacity, 3)
        call resize_transition_integral_arrays (transition_bra_csf_indices, transition_ket_csf_indices, &
            transition_angular_coefficients, transition_coefficient_capacity, 3)
        !
        !   This was allocated in LOAD
        !
        !      CALL DALLOC (PF, 'QF', 'OSLC')
        !      CALL DALLOC (QF, 'QF', 'OSCL')
        !      CALL DALLOC (IVEC, 'IVEC', 'OSCL')
        !      CALL DALLOC (IATJPO, 'IATJPO', 'OSCL')
        !      CALL DALLOC (IASPAR, 'IASPAR', 'OSCL')
        !      CALL DALLOC (EVAL, 'EVAL', 'OSCL')
        !      CALL DALLOC (EVEC, 'EVEC', 'OSCL')
        !
        !   Close all files
        !
        if (mpi_rank == 0) close(summary_output_unit)
        !
        return
        !
        return
        !
    end subroutine calculate_oscillator_strengths
end module calculate_oscillator_strengths_interface
