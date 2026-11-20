module mcp_coefficient_storage
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: mcp_coefficient_count, mcp_integral_count
    ! Consecutive end offsets delimit the coefficients for each packed orbital pair.
    ! Each coefficient has one bra CSF index and one ket CSF index.
    integer, pointer :: mcp_bra_csf_indices(:), mcp_ket_csf_indices(:), mcp_packed_orbital_pairs(:), mcp_coefficient_end_offsets(:)
    real(real64), pointer :: mcp_coefficient_values(:)
end module mcp_coefficient_storage
