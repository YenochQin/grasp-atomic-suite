module rci_prerun_settings
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! Stage 0 uses the full matrix, stage 1 computes prerun weights, and stage 2 applies the cutoffs.
    integer :: prerun_stage, prerun_csf_count
    real(real64) :: prerun_one_body_weight_cutoff, prerun_two_body_weight_cutoff
end module rci_prerun_settings
