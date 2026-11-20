module write_hamiltonian_records_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_hamiltonian_records
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_hamiltonian_records(restart_matrix_unit)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        !                                                                      *
        !   Routine for printing the Hamiltonian matrix.  File restart_matrix_unit must be   *
        !   positioned correctly before a call is made to this module.         *
        !                                                                      *
        !   Call(s) to: [LIB92]: CONVRT.                                       *
        !                                                                      *
        !   Written by Farid A Parpia             Last revision: 14 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use sparse_hamiltonian_storage, only: hamiltonian_values, matrix_partner_indices
        use orbital_csf_basis, only: csf_count
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: restart_matrix_unit

        ! Local declarations
        integer      :: csf_index, nelc, ir, lenthr, lenthc
        real(real64) :: elsto
        character    :: cir*8, cic*8
        !-----------------------------------------------
        !
        do csf_index = 1, csf_count
            read (restart_matrix_unit) nelc, elsto, (hamiltonian_values(ir),ir=1,nelc), (matrix_partner_indices(ir),ir=1,nelc)
            do ir = 1, nelc
                call format_integer (matrix_partner_indices(ir), cir, lenthr)
                call format_integer (csf_index, cic, lenthc)
                call ensure_debug_output()
                write (angular_debug_unit, '(" H (",a,",",a,") = ",1p,1d19.12)') cir(1:lenthr), cic(1:lenthc), &
                    & hamiltonian_values(ir)
            end do
        end do
        !
        return
        !
    end subroutine write_hamiltonian_records
end module write_hamiltonian_records_interface
