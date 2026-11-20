!***********************************************************************
!                                                                      *
real(real64) function breit_radial_integral (itype, orbital_a, orbital_c, orbital_b, orbital_d, k)
    ! Module dependencies
    !                                                                      *
    !   This routine evaluates the transverse interaction integrals:       *
    !                                                                      *
    !      ITYPE = 1: General R (k; a c | b d )                            *
    !            = 2: General S (k; a c | b d )                            *
    !            = 3: R (k; a , b d )                                      *
    !            = 4: F (k; a , b )                                        *
    !            = 5: G (k; a , b ) type integral                          *
    !            = 6: H (k; a , b ) type integral                          *
    !                                                                      *
    !   Call(s) to: [RCI92]: RKINT, SKINT.                                 *
    !                                                                      *
    !                                           Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_radial_points
    use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii

    ! Procedure interfaces
    use coulomb_radial_integral_interface, only: coulomb_radial_integral
    use transverse_radial_integral_interface, only: transverse_radial_integral
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: itype
    integer, intent(in) :: orbital_a
    integer, intent(in) :: orbital_c
    integer, intent(in) :: orbital_b
    integer, intent(in) :: orbital_d
    integer  :: k

    ! Local declarations
    integer :: mxrac, i, mxrbd
    logical :: paired_orbitals, exchanged_orbitals
    real(real64) :: rac(max_radial_points), rbd(max_radial_points)
    !-----------------------------------------------
    !
    mxrac = min(orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
    do i = 1, mxrac
        rac(i) = radial_large_components(i,orbital_a)*radial_small_components(i,orbital_c)
    end do
    !
    mxrbd = min(orbital_last_grid_index(orbital_b),orbital_last_grid_index(orbital_d))
    do i = 1, mxrbd
        rbd(i) = radial_large_components(i,orbital_b)*radial_small_components(i,orbital_d)
    end do
    !
    paired_orbitals = orbital_a == orbital_b .and. orbital_c == orbital_d
    exchanged_orbitals = orbital_a == orbital_d .and. orbital_c == orbital_b
    select case (itype)
    case (2)
        if (paired_orbitals .or. exchanged_orbitals) then
            breit_radial_integral = transverse_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1)
        else
            breit_radial_integral = (transverse_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k, &
                & 1) + transverse_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b, &
                & orbital_d,k,2))*&
                half
        end if
    case (3)
        if (orbital_a == orbital_c) then
            do i = 1, mxrbd
                rbd(i) = rbd(i) + radial_large_components(i,orbital_d)*radial_small_components(i,orbital_b)
            end do
            breit_radial_integral = (coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k, &
                & 0) + coulomb_radial_integral(rbd,orbital_b,orbital_d,rac,orbital_a, &
                & orbital_c,k,0&
                ) + coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,2) + coulomb_radial_integral(rbd, &
                & orbital_b,orbital_d,rac,orbital_a, &
                & orbital_c,k,2)&
                )*half
            return
        end if
        do i = 1, mxrac
            rac(i) = rac(i) + radial_large_components(i,orbital_c)*radial_small_components(i,orbital_a)
        end do
        breit_radial_integral = (coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k, &
            & 1) + coulomb_radial_integral(rbd,orbital_b,orbital_d,rac,orbital_a, &
            & orbital_c,k,1)&
            + coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,0) + coulomb_radial_integral(rbd, &
            & orbital_b,orbital_d,rac,orbital_a,orbital_c, &
            & k,0))*&
            half

    case (4)
        breit_radial_integral = coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k, &
            & 0) + coulomb_radial_integral(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c,k,0)
    case (6)
        breit_radial_integral = transverse_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1)
    case default
        !        Type 1, type 5 and the legacy fallthrough for other values.
        if (itype == 5) then
            paired_orbitals = .not. exchanged_orbitals
        end if
        if (paired_orbitals) then
            breit_radial_integral = two*coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k,1)
        else if (exchanged_orbitals) then
            breit_radial_integral = coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k, &
                & 1) + coulomb_radial_integral(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c,k,1)
        else
            breit_radial_integral = (coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b,orbital_d,k, &
                & 1) + coulomb_radial_integral(rac,orbital_a,orbital_c,rbd,orbital_b, &
                & orbital_d,k,2)&
                + coulomb_radial_integral(rbd,orbital_b,orbital_d,rac,orbital_a,orbital_c,k,1) + coulomb_radial_integral(rbd, &
                & orbital_b,orbital_d,rac,orbital_a,orbital_c, &
                & k,2))*&
                half
        end if
    end select
end function breit_radial_integral
