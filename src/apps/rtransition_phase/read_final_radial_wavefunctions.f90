module read_final_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_final_radial_wavefunctions
contains

    ! This subroutine loads  radial wavefunctions from the  .rwf  file
    ! and performs some related setup.
    !
    ! Written by Per Jonsson                              June 1996
    subroutine read_final_radial_wavefunctions(state_name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use array_memory_management, only: alloc, dalloc
        use atomic_calculation_state, only: max_radial_points, nuclear_charge, speed_of_light_au
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use biorthogonal_orbital_sets, only: final_orbital_symmetry_labels, final_orbital_energy_parameters, &
            final_orbital_origin_exponents, final_orbital_count, final_orbital_principal_n, final_orbital_kappa
        use radial_wavefunctions, only: final_orbital_last_grid_index, final_orbital_origin_parameters, &
            final_radial_large_components, final_radial_small_components
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        character , intent(in) :: state_name*24
        integer :: j, i, k, nwin, io_status, npyff, nakyff, myff, jj
        real(real64) :: con, fkk, eyff, pzy
        real(real64), pointer :: pa(:), qa(:), ra(:)
        character :: radial_file_signature*6
        ! Common relevant for the final state
        write (output_unit, *) "Loading Radial WaveFunction File for final state..."
        ! Open the radial wave function file
        j = index(state_name," ")
        open(newunit=wavefunction_input_unit, file=state_name(1:j-1)//".bw", form="UNFORMATTED", status="OLD", &
            position="asis")
        ! Allocate storage to orbital arrays
        call alloc (final_radial_large_components, max_radial_points, final_orbital_count, "final_radial_large_components", &
            "read_final_radial_wavefunctions")
        call alloc (final_radial_small_components, max_radial_points, final_orbital_count, "final_radial_small_components", &
            "read_final_radial_wavefunctions")
        ! Setup: (1) Orbital arrays to zero
        !        (2) Array E to -1 (no orbitals estimated)
        !        (3) Parameters GAMMA for each orbital
        con = nuclear_charge/speed_of_light_au
        con = con*con
        do j = 1, final_orbital_count
            final_radial_large_components(:max_radial_points,j) = 0.0_real64
            final_radial_small_components(:max_radial_points,j) = 0.0_real64
            final_orbital_energy_parameters(j) = -1.0_real64
            k = abs(final_orbital_kappa(j))
            if (nuclear_model_parameter_count > 0) then
                final_orbital_origin_exponents(j) = dble(k)
            else if (nuclear_model_parameter_count == 0) then
                fkk = dble(k*k)
                if (fkk >= con) then
                    final_orbital_origin_exponents(j) = sqrt(fkk - con)
                else
                    write (output_unit, *) "read_final_radial_wavefunctions: Imaginary gamma parameter"
                    write (output_unit, *) " for ", final_orbital_principal_n(j), final_orbital_symmetry_labels(j), " orbital; the"
                    write (output_unit, *) " point model for the nucleus"
                    write (output_unit, *) " is inappropriate for Z > ", speed_of_light_au, "."
                    stop
                end if
            end if
        end do
        ! Read orbital information from Read Orbitals File;
        nwin = 0
        read (wavefunction_input_unit, iostat=io_status) radial_file_signature
        if (io_status/=0 .or. radial_file_signature/="G92RWF") then
            write (output_unit, *) "This is not a Radial WaveFunction File;"
            close(wavefunction_input_unit)
        end if

        do
            read (wavefunction_input_unit, iostat=io_status) npyff, nakyff, eyff, myff
            if (io_status /= 0) exit
            call alloc (pa, myff, "PA", "read_final_radial_wavefunctions" )
            call alloc (qa, myff, "QA", "read_final_radial_wavefunctions" )
            call alloc (ra, myff, "RA", "read_final_radial_wavefunctions" )
            read (wavefunction_input_unit) pzy, (pa(i),i=1,myff), (qa(i),i=1,myff)
            read (wavefunction_input_unit) (ra(i),i=1,myff)

            do j = 1, final_orbital_count
                if (.not.(final_orbital_energy_parameters(j)<0.0_real64 .and. npyff==final_orbital_principal_n(j) .and. &
                    nakyff==final_orbital_kappa(j)&
                    )) cycle
                final_orbital_origin_parameters(j) = pzy
                final_orbital_energy_parameters(j) = eyff
                final_orbital_last_grid_index(j) = myff
                do jj = 1, final_orbital_last_grid_index(j)
                    final_radial_large_components(jj,j) = pa(jj)
                    final_radial_small_components(jj,j) = qa(jj)
                end do
                nwin = nwin + 1
            end do
            call dalloc (pa, "PA", "read_final_radial_wavefunctions")
            call dalloc (qa, "QA", "read_final_radial_wavefunctions")
            call dalloc (ra, "RA", "read_final_radial_wavefunctions")
        end do
        ! Stop with an error message if all orbitals are not known
        if (nwin < final_orbital_count) then
            write (output_unit, *) "read_final_radial_wavefunctions: All required orbitals not"
            write (output_unit, *) " found."
            stop
        end if
        write (output_unit, *) " ... load complete;"
        close(wavefunction_input_unit)

        return
    end subroutine read_final_radial_wavefunctions
end module read_final_radial_wavefunctions_interface
