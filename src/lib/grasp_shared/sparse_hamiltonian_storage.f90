module sparse_hamiltonian_storage
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none (type, external)

    ! Module interface
    private :: real64, int64

    ! Module state
    ! Each stored column segment ends with its diagonal. End offsets are cumulative within
    ! the local arrays. MPI ranks store only their assigned columns.
    ! Matrix positions use int64. Partner indices identify CSFs, using the basis index kind.
    ! The storage owns all three arrays. Column offsets are allocated with a zero sentinel.
    integer(int64) :: stored_matrix_element_count
    real(real64), allocatable :: hamiltonian_values(:)
    integer(int64), allocatable :: stored_segment_end_offsets(:)
    integer, allocatable :: matrix_partner_indices(:)
    ! RCI stores the latest assembly statistics for reporting and MPI reduction.
    integer :: block_breit_integral_counts(6)
    integer(int64) :: block_stored_element_count
    integer :: block_coulomb_one_body_integral_count, block_coulomb_one_body_coefficient_count, block_slater_integral_count, &
        block_coulomb_two_body_coefficient_count, block_breit_coefficient_count,  &
        block_breit_core_coefficient_count, block_vacuum_polarization_integral_count, block_kinetic_integral_count, &
            block_vinti_integral_count, block_assembly_csf_count
    real(real64) :: block_matrix_element_cutoff
end module sparse_hamiltonian_storage
