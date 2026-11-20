! Load a full CSF basis and compact away unused orbitals outside the closed core.
module csf_full_basis_reader
   use iso_fortran_env, only: int8, int64
   use csf_descriptor_reader, only: descriptor_source, open_descriptor_source, close_descriptor_source, read_descriptor_block
   use rmcdhf_csf_input, only: descriptor_path
   use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations, orbital_principal_n, orbital_kappa, &
                    orbital_l, orbital_two_j, orbital_symmetry_labels
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use atomic_calculation_state, only: electron_count
   use transition_basis_blocks, only: symmetry_block_count => nblock, cumulative_csf_counts => ncfblk
   use grasp_parameters, only: max_orbitals
   implicit none (type, external)
   private
   public :: load_csf_basis
contains
   subroutine load_csf_basis(state_name,core_subshell_count,status)
      character(len=*), intent(in) :: state_name
      integer, intent(out) :: core_subshell_count, status
      type(descriptor_source) :: basis_reader
      integer :: block_index, first_csf, last_csf, source_orbital_index, retained_orbital_count
      integer(int64) :: csf_offset
      integer(int8), pointer :: new_occupations(:,:), new_subshell_quantum_codes(:,:,:), new_couplings(:,:)
      call open_descriptor_source(basis_reader,descriptor_path(state_name),size(cumulative_csf_counts),status=status)
      if (status /= 0) return
      allocate(new_occupations(max_orbitals,basis_reader%total_csf_count),stat=status)
      if (status /= 0) goto 900
      allocate(new_subshell_quantum_codes(max_orbitals,3,basis_reader%total_csf_count),stat=status)
      if (status /= 0) then
         deallocate(new_occupations)
         goto 900
      end if
      allocate(new_couplings(max_orbitals,basis_reader%total_csf_count),stat=status)
      if (status /= 0) then
         deallocate(new_occupations,new_subshell_quantum_codes)
         goto 900
      end if
      csf_offset=0_int64
      do block_index=1,basis_reader%block_count
         first_csf=int(csf_offset+1_int64)
         csf_offset=csf_offset+int(basis_reader%block_lengths(block_index),int64)
         last_csf=int(csf_offset)
         call read_descriptor_block(basis_reader,block_index,new_occupations(:,first_csf:last_csf), &
              new_subshell_quantum_codes(:,:,first_csf:last_csf),new_couplings(:,first_csf:last_csf),status)
         if (status /= 0) then
            deallocate(new_occupations,new_subshell_quantum_codes,new_couplings)
            goto 900
         end if
      end do
      csf_count=basis_reader%total_csf_count
      core_subshell_count=basis_reader%core_count
      electron_count=basis_reader%electron_count
      symmetry_block_count=basis_reader%block_count
      cumulative_csf_counts=0
      last_csf=0
      do block_index=1,symmetry_block_count
         last_csf=last_csf+basis_reader%block_lengths(block_index)
         cumulative_csf_counts(block_index)=last_csf
      end do
      retained_orbital_count=0
      do source_orbital_index=1,basis_reader%orbital_count
         if (source_orbital_index > core_subshell_count) then
            if (all(new_occupations(source_orbital_index,:) == 0)) cycle
         end if
         retained_orbital_count=retained_orbital_count+1
         orbital_principal_n(retained_orbital_count)=basis_reader%principal_numbers(source_orbital_index)
         orbital_kappa(retained_orbital_count)=basis_reader%relativistic_kappas(source_orbital_index)
         orbital_l(retained_orbital_count)=basis_reader%orbital_momenta(source_orbital_index)
         orbital_two_j(retained_orbital_count)=basis_reader%single_electron_two_j(source_orbital_index)
         orbital_symmetry_labels(retained_orbital_count)=basis_reader%orbital_labels(source_orbital_index)
         if (retained_orbital_count /= source_orbital_index) then
            new_occupations(retained_orbital_count,:)=new_occupations(source_orbital_index,:)
            new_subshell_quantum_codes(retained_orbital_count,:,:)=new_subshell_quantum_codes(source_orbital_index,:,:)
         end if
      end do
      if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
      if (associated(csf_subshell_quantum_codes)) deallocate(csf_subshell_quantum_codes)
      if (associated(csf_coupling_quantum_codes)) deallocate(csf_coupling_quantum_codes)
      csf_subshell_occupations=>new_occupations
      csf_subshell_quantum_codes=>new_subshell_quantum_codes
      csf_coupling_quantum_codes=>new_couplings
      orbital_count=retained_orbital_count
      new_occupations(retained_orbital_count+1:,:)=0
      new_subshell_quantum_codes(retained_orbital_count+1:,:,:)=0
900   continue
      call close_descriptor_source(basis_reader)
   end subroutine load_csf_basis
end module csf_full_basis_reader
