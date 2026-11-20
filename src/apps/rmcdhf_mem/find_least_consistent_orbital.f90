module find_least_consistent_orbital_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: find_least_consistent_orbital
contains

    !***********************************************************************
    !                                                                      *
    subroutine find_least_consistent_orbital(orbital_index)
        ! Module dependencies
        !                                                                      *
        !   This subroutine finds the least self-consistent orbital            *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
        !                                                                      *
        ! J initialized to zero
        ! XHH 1997.02.14
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use fixed_orbital_flags, only: is_orbital_fixed
        use orbital_csf_basis, only: orbital_count
        use scf_iteration_state, only: orbital_self_consistency_errors
        implicit none (type, external)

        ! Arguments
        integer , intent(out) :: orbital_index

        ! Local declarations
        integer :: candidate_orbital_index
        real(real64) :: largest_self_consistency_error
        !-----------------------------------------------
        !
        !
        orbital_index = 0
        largest_self_consistency_error = 0.0_real64
        do candidate_orbital_index = 1, orbital_count
            if (is_orbital_fixed(candidate_orbital_index)) cycle
            if (orbital_self_consistency_errors(candidate_orbital_index) <= largest_self_consistency_error) cycle
            largest_self_consistency_error = orbital_self_consistency_errors(candidate_orbital_index)
            orbital_index = candidate_orbital_index
        end do
        !
        return
    end subroutine find_least_consistent_orbital
end module find_least_consistent_orbital_interface
