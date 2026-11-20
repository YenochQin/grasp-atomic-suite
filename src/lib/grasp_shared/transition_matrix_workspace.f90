module transition_matrix_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: transition_matrix_capacity
    real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)
    integer, pointer :: jja(:), jjb(:)
    integer :: transition_multipole_rank, transition_multipole_type
    integer :: transition_coefficient_capacity
    real(real64), pointer :: transition_angular_coefficients(:), total_coulomb_probabilities(:), total_babushkin_probabilities(:)
    integer, pointer :: transition_bra_csf_indices(:), transition_ket_csf_indices(:)
    integer :: stored_transition_integral_count, transition_integral_count
    integer, pointer :: transition_integral_coefficient_offsets(:), packed_transition_integral_labels(:)
    integer :: requested_multipole_count
    integer, pointer :: requested_multipole_codes(:)
    logical :: transition_option_flags(10)
    integer :: initial_mixing_unit, final_mixing_unit
end module transition_matrix_workspace
