module biorthogonal_subshell_layout
    ! Module dependencies
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module constants
    integer, parameter :: max_biorthogonal_kappa_groups = 20

    ! Module state
    ! Group membership maps each orbital to its kappa group.
    ! Each group stores its ordered orbitals and the inverse position map.
    integer :: initial_kappa_group_by_orbital(max_orbitals)
    integer :: initial_orbital_counts_by_group(max_biorthogonal_kappa_groups)
    integer :: initial_orbitals_by_group(max_biorthogonal_kappa_groups, max_biorthogonal_kappa_groups)
    integer :: final_kappa_group_by_orbital(max_orbitals)
    integer :: final_orbital_counts_by_group(max_biorthogonal_kappa_groups)
    integer :: final_orbitals_by_group(max_biorthogonal_kappa_groups, max_biorthogonal_kappa_groups)
    integer :: initial_orbital_positions_by_group(max_biorthogonal_kappa_groups, max_orbitals), &
        final_orbital_positions_by_group(max_biorthogonal_kappa_groups, max_orbitals)
    ! Record each group’s kappa and the core-orbital counts in both bases.
    integer :: initial_core_counts_by_group(max_biorthogonal_kappa_groups), &
        final_core_counts_by_group(max_biorthogonal_kappa_groups), kappa_by_group(max_biorthogonal_kappa_groups)
    integer :: kappa_group_count
    integer :: biorthogonal_angular_unit
    integer :: kappa_coefficient_units(max_biorthogonal_kappa_groups)
contains
    subroutine set_unpaired_rotations(order, orbital, ci)
        use, intrinsic :: iso_fortran_env, only: real64
        integer, intent(in) :: order
        real(real64), intent(out) :: orbital(order, order), ci(order, order)
        integer :: i

        orbital = 0.0_real64
        do i = 1, order
            orbital(i,i) = 1.0_real64
        end do
        ci = orbital
    end subroutine set_unpaired_rotations
end module biorthogonal_subshell_layout
