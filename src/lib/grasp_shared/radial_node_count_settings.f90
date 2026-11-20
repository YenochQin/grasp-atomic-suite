! Stores the radial node-counting threshold and diagnostic output settings.
module radial_node_count_settings
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: relative_extremum_cutoff
    integer :: count_context = 0
    ! The current rank is zero in serial. Only rank zero writes the node trace.
    integer :: count_trace_rank = 0
end module radial_node_count_settings
