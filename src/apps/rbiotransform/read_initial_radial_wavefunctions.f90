!***********************************************************************
!                                                                      *
subroutine read_initial_radial_wavefunctions(state_name, ntestg)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use calculation_state_files, only: orbital_input_path
    use debug_output_settings, only: angular_debug_unit, ensure_debug_output
    !                                                                      *
    !   This subroutine loads  radial wavefunctions from the  .rwf  file   *
    !   and performs some related setup.                                   *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
    !                                                                      *
    !   Written by Per Jonsson                              June 1996      *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals, max_radial_points
    use array_memory_management, only: alloc, dalloc
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
    use atomic_calculation_state,           only: nuclear_charge, speed_of_light_au
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use biorthogonal_subshell_layout,         only: kamax, nshlii, nshlpii
    use radial_wavefunctions,          only: pfii, qfii, pzii, mfii

    ! Procedure interfaces
    use interpolate_initial_orbital_interface, only: interpolate_initial_orbital
    implicit none (type, external)

    ! Arguments and local declarations
    integer :: wavefunction_input_unit
    integer , intent(in) :: ntestg
    character , intent(in) :: state_name*24
    integer  :: nak(max_orbitals), np(max_orbitals)
    integer :: ntestl, ntest, j, k, i, nwin, io_status, npy, naky, my, jj, kk, jjj&
        , kkk, ierr
    real(real64) :: e(max_orbitals), gama(max_orbitals)
    real(real64) :: con, fkk, ey, pzy, dnorm
    real(real64), pointer :: pa(:), qa(:), ra(:)
    character  :: nh(max_orbitals)*2
    character :: radial_file_signature*6
    !-----------------------------------------------
    !
    !  Common relevant for the initial state
    !
    !
    ntestl = 0
    ntest = max0(ntestl,ntestg)
    ntest = 0

    !
    !   Write entry message
    !
    write (output_unit, *) "Loading Radial WaveFunction File for initial state..."
    !
    !   Open the radial wave function file
    !
    open(newunit=wavefunction_input_unit, file=orbital_input_path(state_name), form="UNFORMATTED", status="OLD", &
        position="asis")
    !
    !   Save NAK, NP and NH
    !
    nak(:nwii) = nakii(:nwii)
    np(:nwii) = npii(:nwii)
    nh(:nwii) = nhii(:nwii)
    !
    !   Allocate storage to orbital arrays
    !
    call alloc (pfii, max_radial_points,nwii, "PFII", "LODRWII")
    call alloc (qfii, max_radial_points,nwii, "QFII", "LODRWII")

    !
    con = nuclear_charge/speed_of_light_au
    con = con*con
    !
    do j = 1, nwii
        pfii(:max_radial_points,j) = 0.0_real64
        qfii(:max_radial_points,j) = 0.0_real64
        !
        k = abs(nak(j))
        if (nuclear_model_parameter_count /= 0) cycle
        fkk = dble(k*k)
        if (fkk >= con) then
            gama(j) = sqrt(fkk - con)
        else
            write (output_unit, *) "LODRWF: Imaginary gamma parameter"
            write (output_unit, *) " for ", np(j), nh(j), " orbital; the"
            write (output_unit, *) " point model for the nucleus"
            write (output_unit, *) " is inappropriate for Z > ", speed_of_light_au, "."
            stop
        end if
        !
    end do
    !
    !   Read orbital information from Read Orbitals File;
    !
    nwin = 0
    read (wavefunction_input_unit, iostat=io_status) radial_file_signature
    if (io_status/=0 .or. radial_file_signature/="G92RWF") then
        write (output_unit, *) "This is not a Radial WaveFunction File;"
        close(wavefunction_input_unit)
    end if

    if (ntest >= 100) then
        write (*, *) "******************"
        write (*, *) " Entering lodrwfi"
        write (*, *) "******************"
    end if

    do
        read (wavefunction_input_unit, iostat=io_status) npy, naky, ey, my
        if (io_status /= 0) exit
        call alloc (pa, my, "PA", "LODRWFII")
        call alloc (qa, my, "QA", "LODRWFII")
        call alloc (ra, my, "RA", "LODRWFII")
        read (wavefunction_input_unit) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
        read (wavefunction_input_unit) (ra(i),i=1,my)
        !
        !    Orbital order as defined in kapdata
        !
        jj = 0
        do k = 1, kamax
            if (k > 1) jj = nshlii(k-1) + jj
            do j = 1, nshlii(k)
                kk = nshlpii(k,j)
                if (npy/=np(kk) .or. naky/=nak(kk)) cycle
                jjj = jj + j
                pzii(jjj) = pzy
                eii(jjj) = ey
                nakii(jjj) = nak(kk)
                npii(jjj) = np(kk)
                nhii(jjj) = nh(kk)
                gamaii(jjj) = gama(kk)
                call interpolate_initial_orbital (pa, qa, my, ra, jjj, dnorm)
                if (ntest >= 100) write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') npii(jjj), nhii(jjj), eii(jjj)&
                    , dnorm
                if (ntest > 1000) then
                    write (*, *) "PF              QF             RA"
                    do kkk = 1, mfii(jjj)
                        write (*, *) pfii(kkk,jjj), qfii(kkk,jjj), ra(kkk)
                    end do
                end if
                nwin = nwin + 1
            end do
        end do
        call dalloc (pa, "PA", "LODRWII")
        call dalloc (qa, "QA", "LODRWII")
        call dalloc (ra, "RA", "LODRWII")
    end do
    if (ldbpr(3)) then
        call ensure_debug_output()
        write (angular_debug_unit, *) " orbitals renormalised;"
    end if
    !
    !   Stop with an error message if all orbitals are not known
    !
    if (nwin < nwii) then
        write (output_unit, *) "LODRWF: All required orbitals not"
        write (output_unit, *) " found."
        ierr = 1
    else
        write (output_unit, *) " ... load complete;"
    end if
    close(wavefunction_input_unit)
    if (ntest >= 100) then
        write (*, *) "Sorted order should be the same as from kapdat"
        do j = 1, nwii
            write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') npii(j), nhii(j), eii(j), dnorm
        end do
        write (*, *)
        write (*, *) "*****************"
        write (*, *) " Leaving lodrwfi"
        write (*, *) "*****************"
    end if

    return
    !
    return
    !
end subroutine read_initial_radial_wavefunctions
