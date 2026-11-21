module load_biorthogonal_csf_basis_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_biorthogonal_csf_basis_mpi
contains

    subroutine load_biorthogonal_csf_basis_mpi(name,core_subshell_count,igg)
        ! Module dependencies
        use mpi_csf_basis_input, only: load_csf_basis_mpi
        use orbital_csf_basis, only: csf_subshell_occupations
        use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(in) :: name
        integer :: core_subshell_count, igg

        ! Local declarations
        integer :: status
        if (igg /= 1) then
            if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
            if (associated(csf_subshell_quantum_codes)) deallocate(csf_subshell_quantum_codes)
            if (associated(csf_angular_codes)) deallocate(csf_angular_codes)
        end if
        call load_csf_basis_mpi(trim(name),core_subshell_count)
        status=0
        if (status /= 0) error stop "Biotransform descriptor load failed"
    end subroutine load_biorthogonal_csf_basis_mpi
end module load_biorthogonal_csf_basis_mpi_interface
