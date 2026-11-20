! Stores the bra and ket subshell quantum numbers used by the spin-angular tensor routines.
module jj_tensor_coupling_tables
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! The shell indices are positions in the current tensor coupling.
    integer :: bra_shell_codes_1(7), bra_shell_codes_2(7), bra_shell_codes_3(7), bra_shell_codes_4(7)
    integer :: ket_shell_codes_1(7), ket_shell_codes_2(7), ket_shell_codes_3(7), ket_shell_codes_4(7)
    real(real64) :: bra_shell_values_1(3), bra_shell_values_2(3), bra_shell_values_3(3), bra_shell_values_4(3)
    real(real64) :: ket_shell_values_1(3), ket_shell_values_2(3), ket_shell_values_3(3), ket_shell_values_4(3)
end module jj_tensor_coupling_tables
