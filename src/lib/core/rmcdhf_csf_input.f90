! Load RMCDHF occupations after validating the descriptor header against the active MCP basis.
module rmcdhf_csf_input
   use iso_fortran_env, only: int8, int64
   use csf_descriptor_reader, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                  descriptor_block_label, read_descriptor_block
   use grasp_parameters, only: max_orbitals
   use orbital_csf_basis, only: csf_count, orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
                    orbital_symmetry_labels, csf_subshell_occupations
   use atomic_calculation_state, only: electron_count
   use closed_core_subshells, only: closed_core_subshell_count
   use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block
   use level_symmetry_metadata, only: packed_block_symmetries => jpgg
   implicit none (type, external)
   private
   public :: load_rmcdhf_descriptors, descriptor_path, rmcdhf_csf_symmetry
contains
   ! Resolve a global CSF index through cumulative block lengths.
   integer function rmcdhf_csf_symmetry(csf_index) result(packed_symmetry)
      integer, intent(in) :: csf_index
      integer :: block_index
      integer(int64) :: cumulative_csf_count
      if (csf_index < 1) error stop "Invalid CSF index"
      cumulative_csf_count=0_int64
      do block_index=1,symmetry_block_count
         cumulative_csf_count=cumulative_csf_count+int(csf_counts_by_block(block_index),int64)
         if (int(csf_index,int64) <= cumulative_csf_count) then
            packed_symmetry=packed_block_symmetries(block_index)
            if (packed_symmetry == 0) error stop "Missing CSF block symmetry"
            return
         end if
      end do
      error stop "CSF index exceeds basis"
   end function rmcdhf_csf_symmetry

   ! Accept an explicit Parquet path or derive one from the state name.
   function descriptor_path(state_name) result(path)
      character(len=*), intent(in) :: state_name
      character(len=:), allocatable :: path, trimmed_name
      integer :: name_length
      trimmed_name = trim(state_name)
      name_length = len(trimmed_name)
      if (name_length >= 8) then
         if (trimmed_name(name_length-7:) == ".parquet") then
            path=trimmed_name
            return
         end if
      end if
      path = trimmed_name//"_descriptors.parquet"
   end function descriptor_path

   subroutine load_rmcdhf_descriptors(state_name, core_subshell_count, block_labels, status)
      character(len=*), intent(in) :: state_name
      integer, intent(out) :: core_subshell_count, status
      character(len=8), intent(in) :: block_labels(:)
      type(descriptor_source) :: basis_reader
      integer :: block_index, first_csf, last_csf
      integer(int64) :: csf_offset
      integer(int8), pointer :: new_occupations(:,:)
      call open_descriptor_source(basis_reader, descriptor_path(state_name), size(packed_block_symmetries), status=status)
      if (status /= 0) return
      status = 1
      if (basis_reader%core_count /= closed_core_subshell_count .or. basis_reader%orbital_count /= orbital_count .or. &
          basis_reader%electron_count /= electron_count .or. basis_reader%total_csf_count /= csf_count .or. &
          basis_reader%block_count /= symmetry_block_count) goto 900
      if (lbound(csf_counts_by_block,1) /= 0) goto 900
      if (csf_counts_by_block(0) /= 0) goto 900
      if (size(block_labels) < basis_reader%block_count) goto 900
      if (any(basis_reader%block_lengths /= csf_counts_by_block(1:symmetry_block_count))) goto 900
      do block_index = 1, basis_reader%block_count
         if (block_labels(block_index) /= descriptor_block_label(basis_reader%block_two_j(block_index), &
            & basis_reader%block_parities(block_index))) goto 900
      end do
      allocate(new_occupations(max_orbitals,csf_count), stat=status)
      if (status /= 0) goto 900
      csf_offset = 0_int64
      do block_index = 1, basis_reader%block_count
         first_csf = int(csf_offset+1_int64)
         csf_offset = csf_offset+int(basis_reader%block_lengths(block_index),int64)
         last_csf = int(csf_offset)
         call read_descriptor_block(basis_reader, block_index, new_occupations(:,first_csf:last_csf), status=status)
         if (status /= 0) then
            deallocate(new_occupations)
            goto 900
         end if
      end do
      if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
      csf_subshell_occupations=>new_occupations
      core_subshell_count = basis_reader%core_count
      orbital_principal_n(:orbital_count) = basis_reader%principal_numbers
      orbital_kappa(:orbital_count) = basis_reader%relativistic_kappas
      orbital_l(:orbital_count) = basis_reader%orbital_momenta
      orbital_two_j(:orbital_count) = basis_reader%single_electron_two_j
      orbital_symmetry_labels(:orbital_count) = basis_reader%orbital_labels
      packed_block_symmetries = 0
      packed_block_symmetries(:symmetry_block_count) = (basis_reader%block_two_j+1)*basis_reader%block_parities
900   continue
      call close_descriptor_source(basis_reader)
      if (status /= 0) write(*,*) "RMCDHF descriptor loading or MCP contract failed"
   end subroutine load_rmcdhf_descriptors
end module rmcdhf_csf_input
