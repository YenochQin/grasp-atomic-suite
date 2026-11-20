!***********************************************************************
!                                                                      *
subroutine prepare_isotope_calculation(name)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
    !                                                                      *
    !   Interactively determines the data governing the SMS problem.       *
    !                                                                      *
    !   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
    !               [RCI92]: SETISO, SETRWF.                               *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_state_files, only: orbital_input_path
    use grasp_parameters,    only: max_radial_points
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use atomic_calculation_state, only: accy, max_radial_points, nuclear_charge, speed_of_light_au, vacuum_speed_of_light_au
    use calculation_input_defaults, only: inpci, ndef, ndump
    use first_order_csf_cutoffs, only: iccut
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
    use nuclear_potential_state, only: first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
    use orbital_csf_basis, only: max_orbitals
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use breit_frequency_scale, only: wfact
    use yes_no_input, only: getyn
    use isotope_data_input, only: setiso
    use radial_stencil_weights, only: setqic
    use radial_grid_setup, only: radgrd
    use nuclear_coulomb_potential, only: nucpot
    use radial_wavefunction_input, only: setrwfa
    implicit none (type, external)

    ! Arguments
    character, intent(in) :: name*24

    ! Local declarations
    logical :: yes
    !-----------------------------------------------
    !
    !
    !   Open, check, load data from, and close the  .iso  file
    !
    call setiso ("isodata")
    !
    !   Determine the physical effects specifications
    !

    if (ndef /= 0) then
        write (output_unit, *) "The physical speed of light in"
        write (output_unit, *) " atomic units is", vacuum_speed_of_light_au, ";"
        write (output_unit, *) " revise this value?"
        yes = getyn()
        if (yes) then
            write (output_unit, *) "Enter the revised value:"
            read (input_unit, *) speed_of_light_au
        else
            speed_of_light_au = vacuum_speed_of_light_au
        end if
    else
        speed_of_light_au = vacuum_speed_of_light_au
    end if
    !
    if (ndef /= 0) then
        !
        write (output_unit, *) "Treat contributions of some CSFs"
        write (output_unit, *) " as first-order perturbations?"
        yes = getyn()
        if (yes) then
            lfordr = .true.
            write (output_unit, *) "The contribution of CSFs"
            write (output_unit, *) " 1 -- ICCUT will be treated"
            write (output_unit, *) " variationally; the remainder"
            write (output_unit, *) " perturbatively; enter ICCUT:"
            read (input_unit, *) iccut
        else
            lfordr = .false.
            iccut = 0
        end if
    else
        lfordr = .false.
        iccut = 0
    end if
    !
    !   Determine the parameters controlling the radial grid
    !
    if (nuclear_model_parameter_count == 0) then
        radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
        transformed_grid_step = 0.5_real64**4
        radial_point_count = min(220,max_radial_points)
    else
        !CFF     .. should be Z-dependent
        radial_grid_scale = 2.0e-6_real64/nuclear_charge
        transformed_grid_step = 5.0e-2_real64
        radial_point_count = max_radial_points
    end if
    asymptotic_grid_step = 0.0_real64
    if (ndef /= 0) then
        write (output_unit, *) "The default radial grid parameters"
        write (output_unit, *) " for this case are:"
        write (output_unit, *) " RNT = ", radial_grid_scale, ";"
        write (output_unit, *) " H = ", transformed_grid_step, ";"
        write (output_unit, *) " HP = ", asymptotic_grid_step, ";"
        write (output_unit, *) " N = ", radial_point_count, ";"
        write (output_unit, *) " revise these values?"
        yes = getyn()
        if (yes) then
            write (output_unit, *) "Enter RNT:"
            read (input_unit, *) radial_grid_scale
            write (output_unit, *) "Enter H:"
            read (input_unit, *) transformed_grid_step
            write (output_unit, *) "Enter HP:"
            read (input_unit, *) asymptotic_grid_step
            write (output_unit, *) "Enter N:"
            read (input_unit, *) radial_point_count
        end if
    end if
    !
    !   ACCY is an estimate of the accuracy of the numerical procedures
    !
    accy = transformed_grid_step**6
    !
    !   Set up the coefficients for the numerical procedures
    !
    call setqic
    !
    !   Generate the radial grid and all associated arrays
    !
    call radgrd
    !
    !   Generate $- r \times V_nuc (r)$
    !
    call nucpot
    !
    !   Load the radial wavefunctions
    !
    !      CALL SETRWFA(NAME)
    call setrwfa (orbital_input_path(name))
    !
    return
end subroutine prepare_isotope_calculation
