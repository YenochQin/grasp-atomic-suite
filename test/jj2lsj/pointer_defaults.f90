program pointer_defaults
    use jj2lsj_conversion, only: cs_function_ls, csf_basis_ls, as_function_ls, asf_basis_ls, lsj_list
    implicit none (type, external)
    type(csf_basis_ls) :: basis
    type(as_function_ls) :: state
    type(asf_basis_ls) :: states
    type(lsj_list) :: list
    type(cs_function_ls), allocatable :: rows(:)
    integer :: i
    if (associated(basis%shell).or.associated(basis%csf).or.associated(basis%parent)) error stop "basis defaults"
    if (associated(state%eigenvector).or.associated(states%asf).or.associated(list%items)) error stop "ASF defaults"
    do i=1,2
        allocate(rows(2))
        call check_row(rows(1))
        call check_row(rows(2))
        allocate(rows(1)%occupation(3))
        rows(1)%occupation=2
        if (any(rows(1)%occupation/=2)) error stop "pointer allocation"
        deallocate(rows(1)%occupation)
        if (associated(rows(1)%occupation)) error stop "pointer deallocation"
        deallocate(rows)
    end do
contains
    subroutine check_row(row)
        type(cs_function_ls), intent(in) :: row
        if (associated(row%occupation).or.associated(row%seniority).or.associated(row%w).or. &
            associated(row%shelll).or.associated(row%shells).or.associated(row%shelllx).or. &
            associated(row%shellsx)) error stop "CSF component defaults"
    end subroutine check_row
end program pointer_defaults
