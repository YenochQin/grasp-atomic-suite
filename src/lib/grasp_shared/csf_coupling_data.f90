! Stores subshell quantum numbers and cumulative coupling codes for the active and paired CSF bases.
module csf_coupling_data
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int8
    use grasp_parameters, only: csf_total_j_slot, csf_parity_slot, csf_angular_code_count
    implicit none (type, external)

    ! Module interface
    private :: int8

    ! Module state
    ! Subshell codes are indexed by orbital, then by (seniority, w, 2J+1), then by CSF.
    ! The active arrays use signed bytes. The paired left arrays use default integers.
    integer, allocatable :: left_subshell_quantum_codes(:, :, :)
    integer(int8), allocatable :: csf_subshell_quantum_codes(:, :, :)
    ! The first orbital_count entries are cumulative 2K+1 codes in orbital order.
    ! The total 2J+1 and parity occupy separate slots beyond orbital capacity.
    integer, allocatable :: left_angular_codes(:, :)
    integer(int8), allocatable :: csf_angular_codes(:, :)
end module csf_coupling_data
