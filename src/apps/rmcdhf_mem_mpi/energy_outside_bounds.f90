module energy_outside_bounds_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: energy_outside_bounds
contains

    !***********************************************************************
    !                                                                      *
    logical function energy_outside_bounds (trial_energy)
        ! Module dependencies
        !                                                                      *
        !   This  subprogram determines whether the trial eigenvalue etry is   *
        !   within the bounds (EPSMIN,EPSMAX)                                  *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use scf_iteration_state, only: lower_trial_energy_bound, upper_trial_energy_bound
        implicit none (type, external)

        ! Arguments
        real(real64) , intent(in) :: trial_energy
        !
        !
        if (trial_energy>lower_trial_energy_bound .and. trial_energy<upper_trial_energy_bound) then
            energy_outside_bounds = .false.
        else
            energy_outside_bounds = .true.
        end if
        !
        return
    end function energy_outside_bounds
end module energy_outside_bounds_interface
