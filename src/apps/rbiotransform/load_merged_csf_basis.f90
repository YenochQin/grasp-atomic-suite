module load_merged_csf_basis_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: load_merged_csf_basis
contains

    subroutine load_merged_csf_basis(name,core_subshell_count,igg)
        ! Module dependencies
        use csf_full_basis_reader, only: load_csf_basis
        use orbital_csf_basis, only: csf_subshell_occupations
        use csf_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(in) :: name
        integer :: core_subshell_count, igg

        ! Local declarations
        integer :: status
        if (igg /= 1) then
            if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
            if (associated(csf_subshell_quantum_codes)) deallocate(csf_subshell_quantum_codes)
            if (associated(csf_coupling_quantum_codes)) deallocate(csf_coupling_quantum_codes)
        end if
        call load_csf_basis(trim(name),core_subshell_count,status)
        if (status /= 0) error stop "Biotransform descriptor load failed"
    end subroutine load_merged_csf_basis
end module load_merged_csf_basis_interface
