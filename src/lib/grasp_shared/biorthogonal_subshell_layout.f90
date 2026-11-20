module biorthogonal_subshell_layout
    ! Module dependencies
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module constants
    integer, parameter :: nlmax = 20

    ! Module state
    ! For each basis, nakinv maps orbital to kappa group and nshl counts orbitals per group.
    ! nshlp maps (group, position) to orbital, and nshlpp provides the inverse within each group.
    integer :: nakinvii(max_orbitals)
    integer :: nshlii(nlmax)
    integer :: nshlpii(nlmax, nlmax)
    integer :: nakinvff(max_orbitals)
    integer :: nshlff(nlmax)
    integer :: nshlpff(nlmax, nlmax)
    integer :: nshlppii(nlmax, max_orbitals), nshlppff(nlmax, max_orbitals)
    ! ninii and ninff count core orbitals per group. ikappa stores the corresponding kappa values.
    integer :: ninii(nlmax), ninff(nlmax), ikappa(nlmax)
    integer :: kamax
    integer :: biorthogonal_angular_unit
    integer :: kappa_coefficient_units(nlmax)
end module biorthogonal_subshell_layout
