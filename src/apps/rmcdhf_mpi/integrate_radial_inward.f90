module integrate_radial_inward_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: integrate_radial_inward
contains

    !***********************************************************************
    !                                                                      *
    subroutine integrate_radial_inward(orbital_index, matching_grid_index, large_component_values, small_component_values, &
        & last_radial_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   This program computes the solution of an  inhomogeneous  pair of   *
        !   radial Dirac equations in the tail region. A simple extension of   *
        !   the method of  C Froese Fischer [C Froese, Can J Phys, 41 (1963)   *
        !   1895]  is  used. The  equations  are treated as a boundary value   *
        !   problem, with  the value of P(r) given at the inner boundary and   *
        !   required to  be  sufficiently small  at the outer boundary.  The   *
        !   location  of this second boundary is determined in the course of   *
        !   the calculation. The same  finite  difference equations are used   *
        !   as in the outward integration.                                     *
        !                                                                      *
        !   Arguments:                                                         *
        !                                                                      *
        !      IORB: (Input) Index of orbital                                  *
        !      JP:   (Input) Tabulation point at which the outward             *
        !            integration terminated                                    *
        !      P,Q:  (Input and Output) arrays containing, respectively, the   *
        !            large and small components of the solution                *
        !      MTP:  (Output) Maximum tabulation point of functions P, Q       *
        !                                                                      *
        !   When  written  in  matrix form the system of linear equations iS   *
        !   M*W = V, where  W  is a vector consisting of the elements Q(JP),   *
        !   P(JP+1), Q(JP+1), ... . The matrix  M  is  of band type, with  5   *
        !   elements  in each row.  M is  expressed as a product of two tri-   *
        !   angular matrices L and  U (the L/U decomposition, carried out by   *
        !   the Crout procedure), variable elements of L being stored in the   *
        !   arrays TC, TD  and TE, and  variable elements of U in TH, TI and   *
        !   TJ. The  solution vector W  is then obtained by solving two tri-   *
        !   angular systems, L*Z = V AND U*W = Z.                              *
        !                                                                      *
        !   Call(s) to: [LIB92]: CONVRT.                                       *
        !                                           Last update: 10 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use atomic_calculation_state,            only: working_tolerance
        use radial_grid_state, only: transformed_grid_step, radial_point_count, radial_derivative_over_radius
        use radial_dirac_workspace,            only: large_equation_coefficients, small_equation_coefficients, &
            large_radial_source, small_radial_source
        use orbital_csf_basis,            only: orbital_kappa
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index
        integer , intent(in) :: matching_grid_index
        integer , intent(out) :: last_radial_index
        real(real64) , intent(inout) :: large_component_values(max_radial_points)
        real(real64) , intent(inout) :: small_component_values(max_radial_points)

        ! Local declarations
        integer :: elimination_index, radial_index, text_length, first_zero_grid_index
        real(real64) :: small_component_pivots(max_radial_points), large_component_pivots(max_radial_points), &
            & large_equation_coupling(max_radial_points), small_component_rhs(max_radial_points), &
            large_component_rhs(max_radial_points)
        real(real64) :: tail_tolerance, half_step_times_kappa, radial_kappa_factor, one_plus_grid_kappa, &
            & one_minus_grid_kappa, next_one_plus_grid_kappa, next_one_minus_grid_kappa, &
            & large_from_current_small, current_tail_amplitude, small_from_previous_large, &
            large_from_previous_large, previous_tail_amplitude
        character :: number_text*5
        !-----------------------------------------------
        !
        !   Global initializations
        !
        !ww      EPS = 0.1D 00*ACCY
        tail_tolerance = 0.01_real64*working_tolerance
        half_step_times_kappa = 0.5_real64*transformed_grid_step*dble(orbital_kappa(orbital_index))
        !
        !   Initialize counters
        !
        elimination_index = 1
        radial_index = matching_grid_index
        !
        !   Other initializations
        !
        radial_kappa_factor = half_step_times_kappa*radial_derivative_over_radius(radial_index)
        one_plus_grid_kappa = 1.0_real64 + radial_kappa_factor
        one_minus_grid_kappa = 1.0_real64 - radial_kappa_factor
        radial_kappa_factor = half_step_times_kappa*radial_derivative_over_radius(radial_index+1)
        next_one_plus_grid_kappa = 1.0_real64 + radial_kappa_factor
        next_one_minus_grid_kappa = 1.0_real64 - radial_kappa_factor
        !
        !   Compute required elements of first two rows of  L  and  U
        !
        small_component_pivots(elimination_index) = one_plus_grid_kappa
        large_from_current_small = -large_equation_coefficients(radial_index)/small_component_pivots(elimination_index)
        large_component_pivots(elimination_index) = (-next_one_plus_grid_kappa) + &
            & large_from_current_small*small_equation_coefficients(radial_index+1)
        large_equation_coupling(elimination_index) = next_one_minus_grid_kappa*large_from_current_small - &
            & large_equation_coefficients(radial_index+1)
        !
        !   First elements of solution vector Z
        !
        small_component_rhs(elimination_index) = (-small_radial_source(radial_index)) + small_equation_coefficients(radial_index)* &
            large_component_values(radial_index)
        large_component_rhs(elimination_index) = (-large_radial_source(radial_index)) - &
            & one_minus_grid_kappa*large_component_values(radial_index) - large_from_current_small* &
            small_component_rhs(elimination_index)
        !
        current_tail_amplitude = abs(large_component_rhs(elimination_index)/large_component_pivots(elimination_index))
        !
        do
            elimination_index = elimination_index + 1
            radial_index = radial_index + 1
            !
            !   Failure if tables not long enough
            !
            if (radial_index >= radial_point_count) then
                call format_integer (radial_point_count, number_text, text_length)
                write (output_unit, *) "IN: maximum tabulation point exceeds"
                write (output_unit, *) " dimensional limit (currently "//number_text(1:text_length)//");"
                write (output_unit, *) " radial wavefunction may indicate a"
                write (output_unit, *) " continuum state."
                error stop
            end if
            !
            !   Compute required elements of remaining rows of  L  and  U
            !
            one_plus_grid_kappa = next_one_plus_grid_kappa
            one_minus_grid_kappa = next_one_minus_grid_kappa
            radial_kappa_factor = half_step_times_kappa*radial_derivative_over_radius(radial_index+1)
            next_one_plus_grid_kappa = 1.0_real64 + radial_kappa_factor
            next_one_minus_grid_kappa = 1.0_real64 - radial_kappa_factor
            !
            small_from_previous_large = -small_equation_coefficients(radial_index)/large_component_pivots(elimination_index-1)
            small_component_pivots(elimination_index) = one_plus_grid_kappa - &
                & small_from_previous_large*large_equation_coupling(elimination_index-1)
            large_from_previous_large = one_minus_grid_kappa/large_component_pivots(elimination_index-1)
            large_from_current_small = ((-large_equation_coefficients(radial_index)) - &
                & large_from_previous_large*large_equation_coupling(elimination_index-1)) &
                & /small_component_pivots(elimination_index)
            large_component_pivots(elimination_index) = (-next_one_plus_grid_kappa) + &
                & large_from_current_small*small_equation_coefficients(radial_index+1)
            large_equation_coupling(elimination_index) = (-large_equation_coefficients(radial_index+1)) + &
                & next_one_minus_grid_kappa*large_from_current_small
            !
            !   Solution of  L*Z = V
            !
            small_component_rhs(elimination_index) = (-small_radial_source(radial_index)) - &
                & small_from_previous_large*large_component_rhs(elimination_index-1)
            large_component_rhs(elimination_index) = (-large_radial_source(radial_index)) - &
                & large_from_previous_large*large_component_rhs(elimination_index-1) - &
                    & large_from_current_small*small_component_rhs(elimination_index)
            !
            !   Test for outer boundary
            !
            previous_tail_amplitude = current_tail_amplitude
            current_tail_amplitude = abs(large_component_rhs(elimination_index)/large_component_pivots(elimination_index))
            !cff  Stop if J=NNNP
            if (current_tail_amplitude + previous_tail_amplitude <= tail_tolerance .or. radial_index == max_radial_points-1 ) then
                last_radial_index = radial_index
                exit
            end if
        end do
        !
        !   Reset counter
        !
        elimination_index = elimination_index - 1
        !
        !   Last two rows of solution of  U*W = Z ; evaluation of Q(J)
        !
        small_component_values(radial_index) = 0.0_real64
        large_component_values(radial_index) = large_component_rhs(elimination_index)/large_component_pivots(elimination_index)
        small_component_values(radial_index-1) = &
            & (small_component_rhs(elimination_index) + &
                & small_equation_coefficients(radial_index)*large_component_values(radial_index)) &
            & /small_component_pivots(elimination_index)
        !
        !   Solution of  U*W = Z
        !
        do
            radial_index = radial_index - 1
            elimination_index = elimination_index - 1
            !
            if (elimination_index <= 0) exit
            large_component_values(radial_index) = &
                & (large_component_rhs(elimination_index) &
                & - large_equation_coupling(elimination_index)*small_component_values(radial_index)) &
                & /large_component_pivots(elimination_index)
            small_component_values(radial_index-1) = &
                & (small_component_rhs(elimination_index) + &
                    & (1.0_real64-half_step_times_kappa*radial_derivative_over_radius(radial_index)) &
                & *small_component_values(radial_index) + &
                small_equation_coefficients(radial_index) &
                & *large_component_values(radial_index))/small_component_pivots(elimination_index)
        end do
        !
        !   Complete tables with zeroes
        !
        first_zero_grid_index = last_radial_index + 1
        large_component_values(first_zero_grid_index:radial_point_count) = 0.0_real64
        small_component_values(first_zero_grid_index:radial_point_count) = 0.0_real64
        !
        return
    end subroutine integrate_radial_inward
end module integrate_radial_inward_interface
