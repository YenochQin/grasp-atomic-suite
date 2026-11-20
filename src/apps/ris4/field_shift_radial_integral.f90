!***********************************************************************
!                                                                      *
real(real64) function field_shift_radial_integral (first_orbital, second_orbital)
    ! Module dependencies
    !                                                                      *
    !   The value of RINTDENS is an approximation to:                      *
    !                                                                      *
    !                                                                      *
    !      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
    !                     I     J       I     J                            *
    !                                                                      *
    !   Call(s) to: [SMS92]:  POLINT                                       *
    !                                                                      *
    !   Written by Per Jonsson             Last revision: 24 Dec 1992      *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use atomic_calculation_state,           only: vacuum_speed_of_light_au, pi
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii

    ! Procedure interfaces
    use interpolate_electron_density_interface, only: interpolate_electron_density
    implicit none (type, external)

    ! Arguments
    integer  :: first_orbital, second_orbital

    ! Local declarations
    integer :: l, ii, nmin
    real(real64) :: xa(3), ya(3)
    real(real64) :: dens
    !-----------------------------------------------
    !
    nmin = 10

    ii = 1
    do l = nmin,nmin+4,2
        xa(ii) = radial_grid(l)
        ya(ii) = (radial_large_components(l,first_orbital)*radial_large_components(l, &
            & second_orbital) + radial_small_components(l,first_orbital)* &
            & radial_small_components(l,second_orbital))/(4.0_real64*pi*radial_grid(l)*radial_grid(l))
        ii = ii + 1
    end do
    call interpolate_electron_density (xa,ya,dens)
    field_shift_radial_integral = dens
    !
    return
end function field_shift_radial_integral
