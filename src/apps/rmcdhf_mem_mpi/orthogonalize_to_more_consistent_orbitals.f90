! Project an updated orbital against the more self-consistent orbitals of the same symmetry.
module orthogonalize_to_more_consistent_orbitals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: orthogonalize_to_more_consistent_orbitals
contains

    !    This routine Schmidt orthogonalizes orbital  J  to  all orbitals
    !    which  have  better  self-consistency.  Note that fixed orbitals
    !    have the best self-consistency.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 22 Dec 1992
    !
    !  INV not initialized here but used as
    !  INV = INV + 1
    !  A bug regarding IORDER fixed and new schemes introduced. See comments
    !  Anyway this routine was not used anywhere.
    ! XHH 1997.02.21
    subroutine orthogonalize_to_more_consistent_orbitals(orbital_index, orbital_sign_reversal_count)
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use atomic_calculation_state, only: working_tolerance, max_radial_points
        use orbital_csf_basis, only: orbital_count, orbital_kappa
        use scf_orbital_order, only: orbital_update_order
        use orbital_inversion_flags, only: suppress_radial_inversion
        use fixed_orbital_flags, only: is_orbital_fixed
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap
        use radial_node_count, only: count_radial_nodes
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        integer , intent(inout) :: orbital_sign_reversal_count

        ! Local declarations
        integer :: kraw, k, mtp, nncff
        real(real64) :: eps, ovrlap, original_norm_squared, factor, sgn
        logical :: check, changed
        ! XHH Added common /orba/ and /iounit/

        ! EPS = ACCY*0.1D 00
        eps = working_tolerance*0.01_real64
        check = .not.suppress_radial_inversion(orbital_index)
        ! XHH
        !  Bug fixed. Here J is actually IORDER(J_raw), thus K should be
        !  treated the same way.
        changed = .false.
        do kraw = 1, orbital_count
            k = orbital_update_order(kraw)
            ! XHH orbitals with higher self-consistency are considered
            !      :        (NAK(K) .EQ. NAK(J)) .AND.
            !      :        (SCNSTY(K) .LT. SCNSTY(J)) ) THEN
            ! XHH All orbitals are considered
            !      :        (NAK(K) .EQ. NAK(J)) ) THEN
            ! XHH orbitals before the current and unchanged ones are considered
            if (.not.(orbital_kappa(k)==orbital_kappa(orbital_index) .and. (k<orbital_index .or. is_orbital_fixed(k)))) cycle

            changed = .true.
            ! Compute overlap
            ovrlap = integrate_weighted_radial_overlap(orbital_index,k,0)
            ! Schmidt orthogonalise
            orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index) - &
                & ovrlap*orbital_origin_parameters(k)
            mtp = max(orbital_last_grid_index(orbital_index),orbital_last_grid_index(k))
            radial_large_components(:mtp,orbital_index) = radial_large_components(:mtp, &
                & orbital_index) - ovrlap*radial_large_components(:mtp,k)
            radial_small_components(:mtp,orbital_index) = radial_small_components(:mtp, &
                & orbital_index) - ovrlap*radial_small_components(:mtp,k)
        end do
        ! XHH
        if (changed) then
            ! Normalise
            orbital_last_grid_index(orbital_index) = mtp
            original_norm_squared = integrate_weighted_radial_overlap(orbital_index,orbital_index,0)
            factor = 1.0_real64/sqrt(original_norm_squared)
            ! Determine if inversion is necessary
            if (check) then
                call count_radial_nodes (radial_large_components(:max_radial_points,orbital_index), mtp, nncff, sgn)
                if (sgn < 0.0_real64) then
                    orbital_sign_reversal_count = orbital_sign_reversal_count + 1
                    factor = -factor
                end if
            end if
            ! Perform normalization and/or inversion
            orbital_origin_parameters(orbital_index) = factor*orbital_origin_parameters(orbital_index)
            radial_large_components(2:mtp,orbital_index) = factor*radial_large_components(2:mtp,orbital_index)
            radial_small_components(2:mtp,orbital_index) = factor*radial_small_components(2:mtp,orbital_index)
            ! Find new MF(J)
            mtp = mtp + 1
            do
                mtp = mtp - 1
                if (abs(radial_large_components(mtp,orbital_index)) < eps) then
                    radial_large_components(mtp,orbital_index) = 0.0_real64
                    radial_small_components(mtp,orbital_index) = 0.0_real64
                else
                    orbital_last_grid_index(orbital_index) = mtp
                    exit
                end if
            end do
            ! XHH two lines moved forward and added an ENDIF
        end if
        return
    end subroutine orthogonalize_to_more_consistent_orbitals
end module orthogonalize_to_more_consistent_orbitals_interface
