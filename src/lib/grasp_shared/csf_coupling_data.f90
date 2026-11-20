! Stores subshell quantum numbers and cumulative coupling codes for the active and paired CSF bases.
module csf_coupling_data
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: int8
    implicit none (type, external)

    ! Module interface
    private :: int8

    ! Module state
    ! Subshell codes are indexed by orbital, then by (seniority, w, 2J+1), then by CSF.
    ! The active arrays use signed bytes. The paired left arrays use default integers.
    integer, pointer :: jqsar(:, :, :) => null()
    integer(int8), pointer :: csf_subshell_quantum_codes(:, :, :) => null()
    ! Each coupling column packs nontrivial intermediate 2J+1 codes in sequence. The fixed
    ! last slot (max_orbitals) stores parity times the total 2J+1.
    integer, pointer :: jcupar(:, :) => null()
    integer(int8), pointer :: csf_coupling_quantum_codes(:, :) => null()
end module csf_coupling_data
