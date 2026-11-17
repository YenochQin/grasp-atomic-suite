module rmcdhf_descriptor_input
   use iso_fortran_env, only: int8, int64
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                  descriptor_block_label, read_descriptor_block
   use parameter_def, only: max_orbitals
   use orb_c, only: csf_count, orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
                    orbital_symmetry_labels, csf_subshell_occupations
   use def_c, only: electron_count
   use core_c, only: closed_core_subshell_count
   use hblock_c, only: symmetry_block_count, csf_counts_by_block
   use syma_c, only: jpgg
   implicit none (type, external)
   private
   public :: load_rmcdhf_descriptors, descriptor_path, rmcdhf_csf_symmetry
contains
   integer function rmcdhf_csf_symmetry(csf_index) result(code)
      integer, intent(in) :: csf_index
      integer :: b
      integer(int64) :: offset
      if (csf_index < 1) error stop "Invalid CSF index"
      offset=0_int64
      do b=1,symmetry_block_count
         offset=offset+int(csf_counts_by_block(b),int64)
         if (int(csf_index,int64) <= offset) then
            code=jpgg(b)
            if (code == 0) error stop "Missing CSF block symmetry"
            return
         end if
      end do
      error stop "CSF index exceeds basis"
   end function rmcdhf_csf_symmetry

   function descriptor_path(name) result(path)
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: path, stem
      integer :: n
      stem = trim(name)
      n = len(stem)
      if (n >= 8) then
         if (stem(n-7:) == ".parquet") then
            path=stem
            return
         end if
      end if
      path = stem//"_descriptors.parquet"
   end function descriptor_path

   subroutine load_rmcdhf_descriptors(name, ncore, labels, status)
      character(len=*), intent(in) :: name
      integer, intent(out) :: ncore, status
      character(len=8), intent(in) :: labels(:)
      type(descriptor_source) :: source
      integer :: b, first, last
      integer(int64) :: offset
      integer(int8), pointer :: candidate(:,:)
      call open_descriptor_source(source, descriptor_path(name), size(jpgg), status=status)
      if (status /= 0) return
      status = 1
      if (source%core_count /= closed_core_subshell_count .or. source%orbital_count /= orbital_count .or. &
          source%electron_count /= electron_count .or. source%total_csf_count /= csf_count .or. &
          source%block_count /= symmetry_block_count) goto 900
      if (lbound(csf_counts_by_block,1) /= 0) goto 900
      if (csf_counts_by_block(0) /= 0) goto 900
      if (size(labels) < source%block_count) goto 900
      if (any(source%block_lengths /= csf_counts_by_block(1:symmetry_block_count))) goto 900
      do b = 1, source%block_count
         if (labels(b) /= descriptor_block_label(source%block_two_j(b),source%block_parities(b))) goto 900
      end do
      allocate(candidate(max_orbitals,csf_count), stat=status)
      if (status /= 0) goto 900
      offset = 0_int64
      do b = 1, source%block_count
         first = int(offset+1_int64)
         offset = offset+int(source%block_lengths(b),int64)
         last = int(offset)
         call read_descriptor_block(source, b, candidate(:,first:last), status=status)
         if (status /= 0) then
            deallocate(candidate)
            goto 900
         end if
      end do
      if (associated(csf_subshell_occupations)) deallocate(csf_subshell_occupations)
      csf_subshell_occupations=>candidate
      ncore = source%core_count
      orbital_principal_n(:orbital_count) = source%principal_numbers
      orbital_kappa(:orbital_count) = source%relativistic_kappas
      orbital_l(:orbital_count) = source%orbital_momenta
      orbital_two_j(:orbital_count) = source%single_electron_two_j
      orbital_symmetry_labels(:orbital_count) = source%orbital_labels
      jpgg = 0
      jpgg(:symmetry_block_count) = (source%block_two_j+1)*source%block_parities
900   continue
      call close_descriptor_source(source)
      if (status /= 0) write(*,*) "RMCDHF descriptor loading or MCP contract failed"
   end subroutine load_rmcdhf_descriptors
end module rmcdhf_descriptor_input
