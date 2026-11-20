module transition_output_cutoff
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: transition_probability_cutoff
end module transition_output_cutoff
