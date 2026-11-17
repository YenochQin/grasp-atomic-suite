! Whole-basis adapter used by ordinary property programs.
module csf_basis_input
   use iso_fortran_env, only: int8, int64
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, read_descriptor_block
   use rmcdhf_descriptor_input, only: descriptor_path
   use orb_c, only: csf_count, orbital_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa, &
                    orbital_l, orbital_two_j, orbital_symmetry_labels
   use stat_c, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use def_c, only: electron_count
   use blk_c, only: nblock, ncfblk
   use parameter_def, only: max_orbitals
   implicit none (type, external)
   private
   public :: load_csf_basis
contains
   subroutine load_csf_basis(name,ncore,status)
      character(len=*), intent(in) :: name
      integer, intent(out) :: ncore, status
      type(descriptor_source) :: source
      integer :: b, first, last, i, kept
      integer(int64) :: offset
      integer(int8), pointer :: new_occupations(:,:), new_quantum(:,:,:), new_couplings(:,:)
      call open_descriptor_source(source,descriptor_path(name),size(ncfblk),status=status)
      if (status /= 0) return
      allocate(new_occupations(max_orbitals,source%total_csf_count),stat=status)
      if (status /= 0) goto 900
      allocate(new_quantum(max_orbitals,3,source%total_csf_count),stat=status)
      if (status /= 0) then
         deallocate(new_occupations)
         goto 900
      end if
      allocate(new_couplings(max_orbitals,source%total_csf_count),stat=status)
      if (status /= 0) then
         deallocate(new_occupations,new_quantum)
         goto 900
      end if
      offset=0_int64
      do b=1,source%block_count
         first=int(offset+1_int64)
         offset=offset+int(source%block_lengths(b),int64)
         last=int(offset)
         call read_descriptor_block(source,b,new_occupations(:,first:last), &
              new_quantum(:,:,first:last),new_couplings(:,first:last),status)
         if (status /= 0) then
            deallocate(new_occupations,new_quantum,new_couplings)
            goto 900
         end if
      end do
      csf_count=source%total_csf_count
      ncore=source%core_count
      electron_count=source%electron_count
      nblock=source%block_count
      ncfblk=0
      last=0
      do b=1,nblock
         last=last+source%block_lengths(b)
         ncfblk(b)=last
      end do
      kept=0
      do i=1,source%orbital_count
         if (i > ncore) then
            if (all(new_occupations(i,:) == 0)) cycle
         end if
         kept=kept+1
         orbital_principal_n(kept)=source%principal_numbers(i)
         orbital_kappa(kept)=source%relativistic_kappas(i)
         orbital_l(kept)=source%orbital_momenta(i)
         orbital_two_j(kept)=source%single_electron_two_j(i)
         orbital_symmetry_labels(kept)=source%orbital_labels(i)
         if (kept /= i) then
            new_occupations(kept,:)=new_occupations(i,:)
            new_quantum(kept,:,:)=new_quantum(i,:,:)
         end if
      end do
      if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
      if (associated(csf_subshell_quantum_codes)) deallocate(csf_subshell_quantum_codes)
      if (associated(csf_coupling_quantum_codes)) deallocate(csf_coupling_quantum_codes)
      csf_subshell_occupations=>new_occupations
      csf_subshell_quantum_codes=>new_quantum
      csf_coupling_quantum_codes=>new_couplings
      orbital_count=kept
      new_occupations(kept+1:,:)=0
      new_quantum(kept+1:,:,:)=0
900   continue
      call close_descriptor_source(source)
   end subroutine load_csf_basis
end module csf_basis_input
