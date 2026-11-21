module normal_mass_shift_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: normal_mass_shift_radial_integral
contains

    !***********************************************************************
    !                                                                      *
    subroutine normal_mass_shift_radial_integral (first_orbital,second_orbital,relativistic_nms_integral,leading_nms_integral)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real32
        !                                                                      *
        !   The value of this  function is the one-electron integral I (J,K)   *
        !   for  orbitals  J, K. The analytical expression for this quantity   *
        !   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
        !   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
        !                                                                      *
        !   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
        !                                                                      *
        !   Modified by C. Naz\'e  Oct. 2011                                   *
        !***********************************************************************
        !...Translated by Gediminas Gaigalas 11/18/19
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: radial_workspace_capacity
        use atomic_calculation_state, only: nuclear_charge, speed_of_light_au
        use radial_grid_state, only: radial_workspace_capacity, radial_grid, radial_grid_derivative, transformed_grid_step
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use orbital_csf_basis, only: orbital_kappa, orbital_principal_n, orbital_symmetry_labels
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_capacity, radial_workspace_a, &
            radial_workspace_b
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_wavefunction_derivatives, only: differentiate_radial_wavefunctions
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        real(real64), intent(out) :: relativistic_nms_integral, leading_nms_integral
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital

        ! Local declarations
        real(real64) :: pa(radial_workspace_capacity), pb(radial_workspace_capacity), qa(radial_workspace_capacity), &
            qb(radial_workspace_capacity)
        real(real64) :: a1, a2, piece1, piece2, piece3, piece4
        integer :: kap1, kap2, i, l1, l2, j1, j2, l2_tilda
        !-----------------------------------------------
        !
        !
        !   Set KAP1 and KAP2
        !
        !C*      RINTI_NMS = 0.0D 00
        leading_nms_integral   = 0.0_real64
        relativistic_nms_integral = 0.0_real64
        kap1 = orbital_kappa(first_orbital)
        kap2 = orbital_kappa(second_orbital)
        !
        !   Stop if orbitals J and K have different kappa values
        !
        if (kap1 /= kap2) then
            write (*,'("normal_mass_shift_radial_integral: Attempt to calculate I(",1i2,1a2,",",1i2,1a2,")")') &
                orbital_principal_n(first_orbital),orbital_symmetry_labels(first_orbital),orbital_principal_n(second_orbital), &
                & orbital_symmetry_labels(second_orbital)
            stop
        end if
        !
        !   Determine the l quantum numbers
        !
        if (kap1 > 0) then
            l1 =  kap1
        else
            l1 = -kap1-1
        end if
        !
        if (kap2 > 0) then
            l2 =  kap2
        else
            l2 = -kap2-1
        end if
        if (l1 /= l2) return
        !
        !   Determine the j quantum numbers and l_2 tilda
        !
        j1 = 2*iabs(orbital_kappa(first_orbital))-1
        j2 = 2*iabs(orbital_kappa(second_orbital))-1
        l2_tilda = j2 - l2
        !
        integration_last_index = max (orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        !
        !   Kinetic energy contribution
        !
        !   Piece involving derivatives
        !
        call differentiate_radial_wavefunctions (first_orbital)
        pa(1) = 0.0_real64
        qa(1) = 0.0_real64
        do i = 2,integration_last_index
            pa(i) = radial_workspace_a(i)
            qa(i) = radial_workspace_b(i)
        end do
        call differentiate_radial_wavefunctions (second_orbital)
        pb(1) = 0.0_real64
        qb(1) = 0.0_real64
        do i = 2,integration_last_index
            pb(i) = radial_workspace_a(i)
            qb(i) = radial_workspace_b(i)
        end do
        radial_workspace_a(1) = 0.0_real64
        do i = 2,integration_last_index
            radial_workspace_a(i) = (pa(i)*pb(i)+qa(i)*qb(i))/radial_grid_derivative(i)
        end do
        call integrate_radial_grid (piece1)
        piece1 = piece1/(transformed_grid_step*transformed_grid_step)
        !
        a1 = dble(l2*(l2+1))
        a2 = dble(l2_tilda*(l2_tilda+1))
        radial_workspace_a(1) = 0.0_real64
        do i = 2,integration_last_index
            radial_workspace_a(i) = radial_grid_derivative(i)                                                  &
                *(a1*radial_large_components(i,first_orbital)*radial_large_components(i, &
                & second_orbital)+a2*radial_small_components(i,first_orbital)* &
                & radial_small_components(i,second_orbital))                 &
                /(radial_grid(i)*radial_grid(i))
        end do
        call integrate_radial_grid (piece2)
        !
        !   Pieces not involving derivatives
        !
        radial_workspace_a(1) = 0.0_real64
        do i = 2,integration_last_index
            radial_workspace_a(i) = (radial_small_components(i,first_orbital)*pb(i)+radial_small_components(i, &
                second_orbital)*pa(i))/radial_grid(i)
        end do
        call integrate_radial_grid (piece3)
        piece3 = -2.0_real32*nuclear_charge*piece3/(speed_of_light_au*transformed_grid_step)
        !
        if (kap2 /= 1) then
            radial_workspace_a(1) = 0.0_real64
            do i = 2,integration_last_index
                radial_workspace_a(i) = radial_grid_derivative(i)                                               &
                    * (radial_small_components(i,first_orbital)*radial_large_components(i, &
                    & second_orbital)+radial_large_components(i,first_orbital)* &
                    & radial_small_components(i,second_orbital))                   &
                    /(radial_grid(i)*radial_grid(i))
            end do
            call integrate_radial_grid (piece4)
            piece4 = -piece4*nuclear_charge*(dble (kap2)-1.0_real32)/speed_of_light_au
        else
            piece4 = 0.0_real64
        end if
        leading_nms_integral = 0.5_real32*(piece1+piece2)
        relativistic_nms_integral = 0.5_real32*(piece1+piece2+piece3+piece4)
        !
        !   Debug printout
        !
        !C*      LDBPR(5)=.TRUE.
        if (rmcdhf_debug_flags(5)) then
            call ensure_debug_output()
            write (angular_debug_unit,'(/" K (",1i2,1a2,",",1i2,1a2,") = ",1pd19.12)') &
                orbital_principal_n(first_orbital),orbital_symmetry_labels(first_orbital),orbital_principal_n(second_orbital), &
                & orbital_symmetry_labels(second_orbital),relativistic_nms_integral
        end if
        !
        !C*      RETURN
        !
        !
        return
    end subroutine normal_mass_shift_radial_integral
end module normal_mass_shift_radial_integral_interface
