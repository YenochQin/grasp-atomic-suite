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
    subroutine read_initial_radial_wavefunctions(state_name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use array_memory_management, only: alloc, dalloc
        use atomic_calculation_state, only: max_radial_points, nuclear_charge, speed_of_light_au
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use biorthogonal_orbital_sets, only: initial_orbital_symmetry_labels, initial_orbital_energy_parameters, &
            initial_orbital_origin_exponents, initial_orbital_count, initial_orbital_principal_n, initial_orbital_kappa
        use radial_wavefunctions, only: initial_orbital_last_grid_index, initial_orbital_origin_parameters, &
            initial_radial_large_components, initial_radial_small_components
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: wavefunction_input_unit
        character , intent(in) :: state_name*128
        integer :: j, i, k, nwin, io_status, npyii, nakyii, myii, jj
        real(real64) :: con, fkk, eyii, pzy
        real(real64), pointer :: pa(:), qa(:), ra(:)
        character :: radial_file_signature*6
        ! Common relevant for the initial state
        !
        ! Open the radial wave function file
        j = index(state_name," ")
        open(newunit=wavefunction_input_unit, file=state_name(1:j-1)//".bw", form="UNFORMATTED", status="OLD", &
            position="asis")
        ! Allocate storage to orbital arrays
        call alloc (initial_radial_large_components, max_radial_points, initial_orbital_count, "initial_radial_large_components", &
            "read_initial_radial_wavefunctions")
        call alloc (initial_radial_small_components, max_radial_points, initial_orbital_count, "initial_radial_small_components", &
            "read_initial_radial_wavefunctions")
        con = nuclear_charge/speed_of_light_au
        con = con*con
        write (*, *) "NWII", initial_orbital_count
        do j = 1, initial_orbital_count
            write (*, *) initial_orbital_kappa(j), initial_orbital_principal_n(j), initial_orbital_symmetry_labels(j)
            initial_radial_large_components(:max_radial_points,j) = 0.0_real64
            initial_radial_small_components(:max_radial_points,j) = 0.0_real64
            initial_orbital_energy_parameters(j) = -1.0_real64
            k = abs(initial_orbital_kappa(j))
            if (nuclear_model_parameter_count > 0) then
                initial_orbital_origin_exponents(j) = dble(k)
            else if (nuclear_model_parameter_count == 0) then
                fkk = dble(k*k)
                if (fkk >= con) then
                    initial_orbital_origin_exponents(j) = sqrt(fkk - con)
                else
                    write (output_unit, *) "read_initial_radial_wavefunctions: Imaginary gamma parameter"
                    write (output_unit, *) " for ", initial_orbital_principal_n(j), initial_orbital_symmetry_labels(j), &
                        " orbital; the"
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
            read (wavefunction_input_unit, iostat=io_status) npyii, nakyii, eyii, myii
            if (io_status /= 0) exit
            call alloc (pa, myii, "PA", "read_initial_radial_wavefunctions")
            call alloc (qa, myii, "QA", "read_initial_radial_wavefunctions")
            call alloc (ra, myii, "RA", "read_initial_radial_wavefunctions")
            read (wavefunction_input_unit) pzy, (pa(i),i=1,myii), (qa(i),i=1,myii)
            read (wavefunction_input_unit) (ra(i),i=1,myii)

            do j = 1, initial_orbital_count
                if (.not.(initial_orbital_energy_parameters(j)<0.0_real64 .and. npyii==initial_orbital_principal_n(j) .and. &
                    nakyii==initial_orbital_kappa(j)&
                    )) cycle
                initial_orbital_origin_parameters(j) = pzy
                initial_orbital_energy_parameters(j) = eyii
                initial_orbital_last_grid_index(j) = myii
                do jj = 1, initial_orbital_last_grid_index(j)
                    initial_radial_large_components(jj,j) = pa(jj)
                    initial_radial_small_components(jj,j) = qa(jj)
                end do
                nwin = nwin + 1
            end do
            call dalloc (pa, "PA", "read_initial_radial_wavefunctions")
            call dalloc (qa, "QA", "read_initial_radial_wavefunctions")
            call dalloc (ra, "RA", "read_initial_radial_wavefunctions")
        end do
        ! Stop with an error message if all orbitals are not known
        if (nwin < initial_orbital_count) then
            write (output_unit, *) "read_initial_radial_wavefunctions: All required orbitals not"
            write (output_unit, *) " found."
            stop
        end if

        close(wavefunction_input_unit)

        return
    end subroutine read_initial_radial_wavefunctions
end module read_initial_radial_wavefunctions_interface
