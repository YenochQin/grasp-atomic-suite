! Load a descriptor basis for property programs and stop if the load fails.
module csf_basis_loader
    implicit none (type, external)

    ! Module interface
    private
    public :: load_angular_csf_basis
contains
    subroutine load_angular_csf_basis(state_name,core_subshell_count)
        ! Module dependencies
        use csf_full_basis_reader, only: load_csf_basis
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(in) :: state_name
        integer, intent(out) :: core_subshell_count

        ! Local declarations
        integer :: status
        call load_csf_basis(trim(state_name),core_subshell_count,status)
        if (status /= 0) error stop "Descriptor basis load failed"
    end subroutine load_angular_csf_basis
end module csf_basis_loader
