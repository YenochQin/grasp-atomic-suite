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
    use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
        nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
    use atomic_calculation_state, only: speed_of_light_cm_per_s, electron_count, hartree_to_ev, &
        hartree_to_wavenumber_cm, speed_of_light_au
    use calculation_input_defaults, only: inpci, ndef, ndump
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use orbital_csf_basis, only: max_orbitals
    use transition_matrix_workspace, only: nsdim, hb1, hb2, hc1, hc2, hm1, hm2, jja, jjb, lk, kk, ntdim, &
        xsldr, totc, totb, isldr, isldr1, nint, ninteg, nptr, lab, nkp, kp, ltc
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
    use calculation_header, only: idate, ihed, itime
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use transition_lsj_labels, only: rlev_ener_1, rlev_ener_2, string_csf1, string_csf2, nvectoti, nvectotf, &
        iopen_status1, iopen_status2, lev_pos_1, lev_j_1, lev_par_1, lev_pos_2, lev_j_2, lev_par_2
    use mpi_state, only: mpi_allreduce, mpi_barrier, mpi_comm_world, mpi_double_precision, mpi_error_code, &
        & mpi_rank, mpi_sum

    ! Procedure interfaces
    use prepare_transition_mixing_coefficients_interface, only: prepare_transition_mixing_coefficients
    use resize_transition_matrix_arrays_interface, only: resize_transition_matrix_arrays
    use resize_transition_integral_arrays_interface, only: resize_transition_integral_arrays
    use map_merged_orbitals_interface, only: map_merged_orbitals
    use write_transition_angular_coefficients_interface, only: write_transition_angular_coefficients
    use read_transition_mixing_coefficients_interface, only: read_transition_mixing_coefficients
    use read_transition_angular_coefficients_interface, only: read_transition_angular_coefficients

    ! Module dependencies
    use angular_momentum_triangle, only: itrig

    ! Procedure interfaces
    use build_transition_bessel_functions_interface, only: build_transition_bessel_functions
    use write_transition_csf_weights_interface, only: write_transition_csf_weights
    use write_transition_probability_interface, only: write_transition_probability
    use write_lsj_transition_probability_interface, only: write_lsj_transition_probability

    ! Module dependencies
    use posix_services, only: exit => terminate_process

    ! Procedure interfaces
    use load_initial_lsj_labels_interface, only: load_initial_lsj_labels
    use load_final_lsj_labels_interface, only: load_final_lsj_labels
    use write_input_summary_interface, only: write_input_summary

    ! Module dependencies
    use mpi_change_directory, only: sys_chdir
    implicit none (type, external)

    ! Arguments
    character (len = 24), intent(inout) :: name(2)
    character(len=3)  ::  block_label_text
    character (len = 128), intent(inout) :: fullname(2)
    character(len=128) :: startup_directory,tmpdir

    ! Local constants
    character(len=4), parameter :: iau = "Hart"
    character(len=4), parameter :: iev = " eV "
    character(len=4), parameter :: icm = "Kays"
    character(len=4), parameter :: ihz = " Hz "
    character(len=4), parameter :: iang = " A  "
    integer, parameter :: nca = 65536

    ! Local declarations
    integer :: copied_mixing_unit
    integer :: kka, jkp, transition_angular_unit, io_status, my_nelec, ncftoti, nwi, &
        nvecsizi, iopar, iblki, ncftotf, nwf,&
        nvecsizf, iblkf, nvecpr, i, m,       &
        ielec, levii, levff, itkpo, itest, nlp, lines, &
        ilbl, inum_ii, inum_ff, icount1, icount2, iii
    real(real64) :: factor, omega, argu, asfa, asfb, asfatmp, asfbtmp
    real(real64), pointer :: et(:), et1(:), ipr(:), ipr1(:), next(:)
    logical :: avail, lsame
    character :: answ*1, iunits*4, mixing_file_signature*6
    !-----------------------------------------------
    !
    ! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
    ! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
    !
    ! for the case that the inital file and final file are same
    !
    lsame = trim(fullname(1)) == trim(fullname(2))
    if (lsame.and.mpi_rank==0) call prepare_transition_mixing_coefficients (fullname, inpci)
    call mpi_barrier(mpi_comm_world,mpi_error_code)
    !
    ! write header for the result file
    !
    !cjb fort.24
    if (mpi_rank == 0) then
        if (lsame) then
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
        hm2, lab, nptr, nsdim, 1)
    call resize_transition_integral_arrays (isldr, isldr1, xsldr, ntdim, 1)
    !
    !   Make a connection between the orbitals of the
    !   merged list and the initial and final state lists
    !
    call map_merged_orbitals
    !
    !   Set up units for printing transition energy
    !
    if (ltc(1)) then
        !
        !   Print transition energies in Angstroms
        !
        factor = hartree_to_wavenumber_cm
        factor = 1.0e8_real64/factor
        iunits = iang
    !
    else if (ltc(2)) then
        !
        !   Print energies in eV
        !
        factor = hartree_to_ev
        iunits = iev
    !
    else if (ltc(3)) then
        !
        !   Print transition energies in Hartree Atomic Units
        !
        factor = 1.0_real64
        iunits = iau
    !
    else if (ltc(4)) then
        !
        !   Print transition energies in Hz
        !
        factor = hartree_to_wavenumber_cm
        factor = factor*speed_of_light_cm_per_s
        iunits = ihz
    !
    else if (ltc(5)) then
        !
        !   Print transition energies in Kaysers
        !
        factor = hartree_to_wavenumber_cm
        iunits = icm
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
    do jkp = 1, nkp
        !
        ! read the head of the file of mixing coef. for initial
        !
        if (lsame) fullname(1) = trim(fullname(2))//"_CP"
        if (inpci == 0) then
            open(newunit=initial_mixing_unit,file=trim(fullname(1))//".cbm",              &
                form="UNFORMATTED",status="OLD")
        else
            open(newunit=initial_mixing_unit, file=trim(fullname(1))//".bm",              &
                form="UNFORMATTED",status="OLD")
        end if
        if (lsame) fullname(1) = fullname(2)
        read (initial_mixing_unit, iostat=io_status) mixing_file_signature
        if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
            write (*, *) "Not a GRASP mixing file"
            stop
        end if
        read (initial_mixing_unit) electron_count, ncftoti, nwi, nvectoti, nvecsizi, nblocki
        write (*, *) "   nelec  = ", my_nelec
        write (*, *) "   ncftoti = ", ncftoti
        write (*, *) "   nwi     = ", nwi
        write (*, *) "   nblocki = ", nblocki
        write (*, *)
        !   lbl
        icount1 = 0
        call load_initial_lsj_labels (fullname(1))
        ! If not available generate angular coefficients for all pares of blocks
        call write_transition_angular_coefficients (iopar, jkp, name, transition_angular_unit)
        do iblki = 1, nblocki
            call read_transition_mixing_coefficients (name, inpci, 1)
            !   lbl
            icount1 = nvecii + icount1
            !
            ! read the head of the file of mixing coef. for final
            !
            if (inpci == 0) then
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
            read (final_mixing_unit) electron_count, ncftotf, nwf, nvectotf, nvecsizf, nblockf
            write (*, *) "   nelec  = ", my_nelec
            write (*, *) "   ncftotf = ", ncftotf
            write (*, *) "   nwf     = ", nwf
            write (*, *) "   nblockf = ", nblockf
            write (*, *)
            !GG  lbl
            if(iblki == 1) then
                call load_final_lsj_labels (fullname(2))
                if(mpi_rank==0.and.iopen_status1==0.and.iopen_status2==0) then
                    iii = len_trim(startup_directory)
                    call sys_chdir(trim(startup_directory),iii,mpi_error_code)
                    if (mpi_error_code/=0) then
                        print *, "error changing dir!"
                        call exit(1)
                    end if
                    ilbl = 1
                    call write_input_summary(name,inpci,ilbl)
                    iii = len_trim(tmpdir)
                    call sys_chdir(trim(tmpdir)//"/"//block_label_text,iii+4,mpi_error_code)
                    if (mpi_error_code/=0) then
                        print *, "Error changing dir!"
                        call exit(1)
                    end if
                    write(lsj_transition_output_unit,*) "Transition between files:"
                    write(lsj_transition_output_unit,*) trim(name(1))
                    write(lsj_transition_output_unit,*) trim(name(2))
                end if
            end if
            icount2 = 0
            !GG   lbl end
            do iblkf = 1, nblockf
                call read_transition_mixing_coefficients (name, inpci, 2)
                !   lbl
                icount2 = nvecff + icount2
                !
                !   Allocate storage
                !
                call alloc (totb, nvecff, "TOTB", "OSCL")
                call alloc (totc, nvecff, "TOTC", "OSCL")
                !
                nvecpr = nvecii*nvecff
                call alloc (et, nvecpr, "ET", "OSCL" )
                call alloc (et1, nvecpr, "ET1","OSCL" )
                call alloc (ipr, nvecpr, "IPR", "OSCL" )
                call alloc (ipr1, nvecpr,"IPR1", "OSCL" )
                call alloc (next, nvecpr, "NEXT", "OSCL")
                !
                !   Initialization for total decay rate
                !
                totc(:nvecff) = 0.0_real64
                totb(:nvecff) = 0.0_real64
                !
                call read_transition_angular_coefficients (iopar, jkp, name, transition_angular_unit)
                !
                if (lk > 0) then
                    ielec = (-1)**lk
                    if (ielec == iopar) then
                        kk = 0
                        kka = 0
                    else
                        kk = 1
                        ielec = -ielec
                    end if
                    !
                    !   Set up list of levels for calculation of oscillator strengths
                    !   sort list into increasing order of energy if option 6 set
                    !
                    if (mpi_rank==0.and.iblki==1.and.iblkf==1) then
                        if (kk == 0) then
                            write (summary_output_unit, '(/,/," Electric 2**(",i2,")-pole transitions")') lk
                        else
                            write (summary_output_unit, '(/,/," Magnetic 2**(",i2,")-pole transitions")') lk
                        end if
                        write (summary_output_unit, '(1x,33("="))')
                        if (ltc(1)) then
                            write (summary_output_unit, '(/,"   Upper state        Lower state  ",8x,"Gauge",8x,"Wavelength",13x,&
                                &"Einstein coefficients",13x,"Oscillator")')
                            if (.not.ltc(7)) then
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
                            if (.not.ltc(7)) then
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
                    do levii = 1, nvecii
                        !   lbl
                        inum_ii = icount1 - nvecii + levii
                        do levff = 1, nvecff
                            !   lbl
                            inum_ff = icount2 - nvecff + levff
                            !
                            !   Check for consistent parity and J
                            !
                            itkpo = lk + lk + 1
                            if (itrig(iatjpoii(levii),iatjpoff(levff),itkpo) == 0) &
                                cycle
                            itest = iasparii(levii)*iasparff(levff)*ielec
                            if (itest < 0) cycle
                            !
                            !   Calculate and print transition probability data
                            !
                            nlp = 70 - 8
                            lines = nlp
                            !
                            if (lines >= nlp) lines = 0
                            !
                            m = levff + nvecff*(levii - 1)
                            !             M = LEVFF+NVECII*(LEVII-1)
                            et(m) = evalff(levff) + eavff - evalii(levii) - eavii
                            if (lsame .and. et(m)<=0.0_real32) cycle
                            omega = -et(m)
                            argu = omega/speed_of_light_au
                            call build_transition_bessel_functions (argu)
                            !
                            !  Calculate oscillator strength between the ASFs
                            !
                            call write_transition_csf_weights (asfatmp,asfbtmp,levii,levff)
                            call mpi_allreduce(asfatmp,asfa,1,mpi_double_precision, &
                                mpi_sum, mpi_comm_world, mpi_error_code)
                            call mpi_allreduce(asfbtmp,asfb,1,mpi_double_precision, &
                                mpi_sum, mpi_comm_world, mpi_error_code)
                            if (mpi_rank == 0) then
                                call write_transition_probability(asfa,asfb,levii,levff,omega,factor, &
                                    lines,lsame)
                                if(iopen_status1==0 .and. iopen_status2 ==0) then
                                    call write_lsj_transition_probability (inum_ii,inum_ff,asfa,asfb,levii,levff,&
                                        omega,factor)
                                end if
                            end if
                            !              WRITE (24,317)
                        end do
                    end do
                end if
                !
                !   Deallocate storage; this is local to OSCL
                !
                call dalloc (totb, "TOTB", "OSCL")
                call dalloc (totc, "TOTC", "OSCL")
                !
                call dalloc (et, "ET", "OSCL")
                call dalloc (et1, "ET!", "OSCL")
                call dalloc (ipr, "IPR", "OSCL")
                call dalloc (ipr1, "IPR1", "OSCL")
                call dalloc (next, "NEXT", "OSCL")
                !
                !   Deallocate storage; this is allocated in READMIX
                !
                call dalloc (evalff, "EVALFF", "OSCL")
                call dalloc (evecff, "EVECFF", "OSCL")
                call dalloc (ivecff, "IVECFF", "OSCL")
                call dalloc (iatjpoff, "IATJPOFF", "OSCL")
                call dalloc (iasparff, "IASPARFF", "OSCL")

            end do
            close(final_mixing_unit)
            !
            !   Deallocate storage; this is allocated in READMIX
            !
            call dalloc (evalii, "EVALII", "OSCL")
            call dalloc (evecii, "EVECII", "OSCL")
            call dalloc (ivecii, "IVECII", "OSCL")
            call dalloc (iatjpoii, "IATJPOII", "OSCL")
            call dalloc (iasparii, "IASPARII", "OSCL")

        end do
        close(initial_mixing_unit)
        close(transition_angular_unit)
    end do
    call dalloc (kp, "KP", "OSCL")
    !
    ! close and delete duplicated mixing file
    if (lsame) then
        if (inpci == 0) then
            open(newunit=copied_mixing_unit, file=trim(fullname(2))//"_CP.cbm")
        else
            open(newunit=copied_mixing_unit, file=trim(fullname(2))//"_CP.bm")
        end if
        !GG         CLOSE(68, STATUS='delete')
    end if
    !
    call resize_transition_matrix_arrays (jja, jjb, hb1, hb2, hc1, hc2, hm1, hm2, lab, nptr, nsdim, 3)
    call resize_transition_integral_arrays (isldr, isldr1, xsldr, ntdim, 3)
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
