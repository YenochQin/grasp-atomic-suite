subroutine setcsla(name,ncore,igg)
   use csf_basis_mpi, only: load_csf_basis_mpi
   use orb_c, only: csf_subshell_occupations
   use stat_c, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   implicit none (type, external)
   character(len=24), intent(in) :: name
   integer :: ncore, igg
   integer :: status
   if (igg /= 1) then
      if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
      if (associated(csf_subshell_quantum_codes)) deallocate(csf_subshell_quantum_codes)
      if (associated(csf_coupling_quantum_codes)) deallocate(csf_coupling_quantum_codes)
   end if
   call load_csf_basis_mpi(trim(name),ncore)
      status=0
   if (status /= 0) error stop "Biotransform descriptor load failed"
end subroutine setcsla
