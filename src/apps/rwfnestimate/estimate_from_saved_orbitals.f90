

!***********************************************************************
!                                                                      *
subroutine estimate_from_saved_orbitals(matched_orbital_indices, matched_subshell_count, file_path)
    ! Module dependencies
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This subroutine loads  radial wavefunctions from the  .rwf  file   *
    !   and performs some related setup.                                   *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, OPENFL.                *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use array_memory_management, only: alloc, dalloc
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use initial_orbital_selection, only: max_orbitals, set
    use orbital_csf_basis, only: orbital_energy_parameters, orbital_principal_n, orbital_kappa, orbital_symmetry_labels
    use radial_wavefunctions, only: orbital_origin_parameters
    use initial_orbital_source, only: source
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use file_open, only: openfl
    use radial_orbital_interpolation, only: intrpq
    implicit none (type, external)

    ! Arguments and local declarations
    integer :: wavefunction_input_unit
    integer , intent(in) :: matched_subshell_count
    character  :: file_path*(*)
    integer , intent(in) :: matched_orbital_indices(max_orbitals)
    integer :: open_status, io_status, npy, naky, my, j, loc, i, lenth
    real(real64) :: ey, dnorm
    real(real64), pointer :: pa(:), qa(:), ra(:)
    logical :: found
    character :: radial_file_signature*6
    !-----------------------------------------------
    !
    !

    !FORM = 'UNFORMATTED'
    !STATUS = 'OLD'

    call openfl (wavefunction_input_unit, file_path, "UNFORMATTED", "OLD", open_status)
    if (open_status == 1) then
        write (stderr_unit, *) 'Error openning file "', file_path(1:len_trim(file_path)), &
            '"'
        close(wavefunction_input_unit)
        stop
    end if
    !
    !   Check the file; if not as expected, try again
    !
    read (wavefunction_input_unit, iostat=io_status) radial_file_signature
    if (io_status/=0 .or. radial_file_signature/="G92RWF") then
        write (stderr_unit, *) "This is not a Radial WaveFunction File;"
        close(wavefunction_input_unit)
        stop
    end if
    !
    !   Read orbital information from Read Orbitals File; write summary
    !   to  .dbg  file if option set
    !
    if (ldbpr(3)) then
        call ensure_debug_output()
        write (angular_debug_unit, '(/,"From SUBROUTINE FRMRWF:"/," Orbital",8x,"Eigenvalue",19x,"Norm")')
    end if
    do
        found = .false.
        read (wavefunction_input_unit, iostat=io_status) npy, naky, ey, my
        if (io_status /= 0) exit
        do j = 1, matched_subshell_count
            loc = matched_orbital_indices(j)
            if (.not.(.not.set(loc) .and. orbital_principal_n(loc)==npy .and. orbital_kappa(loc)==naky)) &
                cycle
            found = .true.
            orbital_energy_parameters(loc) = ey
            call alloc (pa, my, "PA", "FRMFRW")
            call alloc (qa, my, "QA", "FRMFRW")
            call alloc (ra, my, "RA", "FRMFRW")
            read (wavefunction_input_unit) orbital_origin_parameters(loc), (pa(i),i=1,my), (qa(i),i=1,my)
            read (wavefunction_input_unit) (ra(i),i=1,my)
            call intrpq (pa, qa, my, ra, loc, dnorm)
            if (ldbpr(3)) then
                call ensure_debug_output()
                write (angular_debug_unit, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') &
                    orbital_principal_n(loc), orbital_symmetry_labels(loc), &
                    & orbital_energy_parameters(loc), dnorm
            end if
            call dalloc (pa, "PA", "FRMFRW")
            call dalloc (qa, "QA", "FRMFRW")
            call dalloc (ra, "RA", "FRMFRW")
            lenth = len_trim(file_path)
            set(loc) = .true.
            source(loc) = file_path(1:3)
            exit
        end do
        if (.not.found) then
            read (wavefunction_input_unit)
            read (wavefunction_input_unit)
        end if
    end do
    if (ldbpr(3)) then
        call ensure_debug_output()
        write (angular_debug_unit, *) " orbitals renormalised;"
    end if
    !
    close(wavefunction_input_unit)
    !
    return
    !
    return
    !
end subroutine estimate_from_saved_orbitals
