module read_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_initial_radial_wavefunctions
contains

    ! This subroutine loads  radial wavefunctions from the  .rwf  file
    ! and performs some related setup.
    !
    ! Written by Per Jonsson                              June 1996
    subroutine read_initial_radial_wavefunctions(state_name, global_trace_level)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_state_files, only: orbital_input_path
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals, max_radial_points
        use array_memory_management, only: alloc, dalloc
        use biorthogonal_orbital_sets, only: max_orbitals, initial_orbital_symmetry_labels, &
            initial_orbital_energy_parameters, initial_orbital_origin_exponents, initial_orbital_count, &
            initial_orbital_principal_n, initial_orbital_kappa
        use atomic_calculation_state,           only: nuclear_charge, speed_of_light_au
        use debug_output_settings, only: rmcdhf_debug_flags
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use biorthogonal_subshell_layout,         only: kappa_group_count, initial_orbital_counts_by_group, &
            initial_orbitals_by_group
        use radial_wavefunctions,          only: initial_radial_large_components, initial_radial_small_components, &
            initial_orbital_origin_parameters, initial_orbital_last_grid_index

        ! Procedure interfaces
        use interpolate_initial_orbital_interface, only: interpolate_initial_orbital
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        integer, intent(in) :: global_trace_level
        character, intent(in) :: state_name*128
        integer  :: nak(max_orbitals), np(max_orbitals)
        integer :: local_trace_level, trace_level, j, k, i, nwin, io_status, npy, naky, my, jj, kk, jjj&
            , kkk, ierr
        real(real64) :: gama(max_orbitals)
        real(real64) :: con, fkk, ey, pzy, dnorm
        real(real64), pointer :: pa(:), qa(:), ra(:)
        character  :: nh(max_orbitals)*2
        character :: radial_file_signature*6
        ! Common relevant for the initial state
        local_trace_level = 0
        trace_level = max0(local_trace_level,global_trace_level)
        trace_level = 0

        write (output_unit, *) "Loading Radial WaveFunction File for initial state..."
        ! Open the radial wave function file
        open(newunit=wavefunction_input_unit, file=orbital_input_path(state_name), form="UNFORMATTED", status="OLD", &
            position="asis")
        ! Save NAK, NP and NH
        nak(:initial_orbital_count) = initial_orbital_kappa(:initial_orbital_count)
        np(:initial_orbital_count) = initial_orbital_principal_n(:initial_orbital_count)
        nh(:initial_orbital_count) = initial_orbital_symmetry_labels(:initial_orbital_count)
        ! Allocate storage to orbital arrays
        call alloc (initial_radial_large_components, max_radial_points,initial_orbital_count, "initial_radial_large_components", &
            "read_initial_radial_wavefunctions")
        call alloc (initial_radial_small_components, max_radial_points,initial_orbital_count, "initial_radial_small_components", &
            "read_initial_radial_wavefunctions")

        con = nuclear_charge/speed_of_light_au
        con = con*con
        do j = 1, initial_orbital_count
            initial_radial_large_components(:max_radial_points,j) = 0.0_real64
            initial_radial_small_components(:max_radial_points,j) = 0.0_real64
            k = abs(nak(j))
            if (nuclear_model_parameter_count /= 0) cycle
            fkk = dble(k*k)
            if (fkk >= con) then
                gama(j) = sqrt(fkk - con)
            else
                write (output_unit, *) "read_initial_radial_wavefunctions: Imaginary gamma parameter"
                write (output_unit, *) " for ", np(j), nh(j), " orbital; the"
                write (output_unit, *) " point model for the nucleus"
                write (output_unit, *) " is inappropriate for Z > ", speed_of_light_au, "."
                stop
            end if
        end do
        ! Read orbital information from Read Orbitals File;
        nwin = 0
        read (wavefunction_input_unit, iostat=io_status) radial_file_signature
        if (io_status/=0 .or. radial_file_signature/="G92RWF") then
            write (output_unit, *) "This is not a Radial WaveFunction File;"
            close(wavefunction_input_unit)
        end if

        if (trace_level >= 100) then
            write (*, *) "******************"
            write (*, *) " Entering lodrwfi"
            write (*, *) "******************"
        end if

        do
            read (wavefunction_input_unit, iostat=io_status) npy, naky, ey, my
            if (io_status /= 0) exit
            call alloc (pa, my, "PA", "read_initial_radial_wavefunctions")
            call alloc (qa, my, "QA", "read_initial_radial_wavefunctions")
            call alloc (ra, my, "RA", "read_initial_radial_wavefunctions")
            read (wavefunction_input_unit) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
            read (wavefunction_input_unit) (ra(i),i=1,my)
            ! Orbital order as defined in kapdata
            jj = 0
            do k = 1, kappa_group_count
                if (k > 1) jj = initial_orbital_counts_by_group(k-1) + jj
                do j = 1, initial_orbital_counts_by_group(k)
                    kk = initial_orbitals_by_group(k,j)
                    if (npy/=np(kk) .or. naky/=nak(kk)) cycle
                    jjj = jj + j
                    initial_orbital_origin_parameters(jjj) = pzy
                    initial_orbital_energy_parameters(jjj) = ey
                    initial_orbital_kappa(jjj) = nak(kk)
                    initial_orbital_principal_n(jjj) = np(kk)
                    initial_orbital_symmetry_labels(jjj) = nh(kk)
                    initial_orbital_origin_exponents(jjj) = gama(kk)
                    call interpolate_initial_orbital (pa, qa, my, ra, jjj, dnorm)
                    if (trace_level >= 100) write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') &
                        initial_orbital_principal_n(jjj), initial_orbital_symmetry_labels(jjj), &
                        initial_orbital_energy_parameters(jjj)&
                        , dnorm
                    if (trace_level > 1000) then
                        write (*, *) "PF              QF             RA"
                        do kkk = 1, initial_orbital_last_grid_index(jjj)
                            write (*, *) initial_radial_large_components(kkk,jjj), initial_radial_small_components(kkk,jjj), ra(kkk)
                        end do
                    end if
                    nwin = nwin + 1
                end do
            end do
            call dalloc (pa, "PA", "read_initial_radial_wavefunctions")
            call dalloc (qa, "QA", "read_initial_radial_wavefunctions")
            call dalloc (ra, "RA", "read_initial_radial_wavefunctions")
        end do
        if (rmcdhf_debug_flags(3)) then
            call ensure_debug_output()
            write (angular_debug_unit, *) " orbitals renormalised;"
        end if
        ! Stop with an error message if all orbitals are not known
        if (nwin < initial_orbital_count) then
            write (output_unit, *) "read_initial_radial_wavefunctions: All required orbitals not"
            write (output_unit, *) " found."
            ierr = 1
        else
            write (output_unit, *) " ... load complete;"
        end if
        close(wavefunction_input_unit)
        if (trace_level >= 100) then
            write (*, *) "Sorted order should be the same as from kapdat"
            do j = 1, initial_orbital_count
                write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') initial_orbital_principal_n(j), &
                    initial_orbital_symmetry_labels(j), initial_orbital_energy_parameters(j), dnorm
            end do
            write (*, *)
            write (*, *) "*****************"
            write (*, *) " Leaving lodrwfi"
            write (*, *) "*****************"
        end if

        return
        return
    end subroutine read_initial_radial_wavefunctions
end module read_initial_radial_wavefunctions_interface
