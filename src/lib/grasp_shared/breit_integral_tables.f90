module breit_integral_tables
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: breit_integral_capacities(6), breit_integral_counts(6)
    integer, pointer :: breit_integral_indices_1(:), breit_integral_indices_2(:), breit_integral_indices_3(:), &
        breit_integral_indices_4(:), &
        breit_integral_indices_5(:), breit_integral_indices_6(:)
    real(real64), pointer :: breit_integral_values_1(:), breit_integral_values_2(:), breit_integral_values_3(:), &
        breit_integral_values_4(:), &
        breit_integral_values_5(:), breit_integral_values_6(:)
    logical :: first_breit_integral_lookup(6)
end module breit_integral_tables
