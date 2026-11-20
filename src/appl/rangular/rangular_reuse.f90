! Exact, order-preserving memoization for the ordinary scalar/Coulomb pair.
! This module owns no scientific cutoffs and never changes the MCP writer.
module rangular_reuse
   use iso_fortran_env, only: int64
   use, intrinsic :: iso_fortran_env, only: real64
   use grasp_parameters, only: max_orbitals
   use orbital_csf_basis, only: orbital_count, csf_count, csf_subshell_occupations, orbital_kappa
   use csf_coupling_data, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   use angular_coefficient_buffer_state, only: &
       angular_coefficient_count, angular_buffer_capacity, angular_coefficient_labels, angular_coefficients
   use debug_output_settings, only: ibug1, ibug2, ibug3, ldbpa, ldbpg
   use one_electron_scalar_coefficients, only: onescalar
   use two_electron_angular_coefficients, only: rkco_gg
   use closed_shell_mcp, only: cord
   use angular_coefficient_buffer, only: alcbuf
   implicit none (type, external)
   private
   public :: configure_angular_reuse, begin_reuse_block, end_reuse_block, load_reused_pair, reuse_requested

   type :: pair_entry
      integer(int64) :: hash
      integer :: key_start, key_size, term_start, term_count, first, second, bra, ket
      real(real64) :: scalar
   end type pair_entry
   type(pair_entry), allocatable :: entries(:)
   integer, allocatable :: buckets(:), term_labels(:,:)
   integer(int64), allocatable :: keys(:), coupling_hash(:)
   real(real64), allocatable :: values(:)
   integer :: entry_count, key_count, term_count, report_unit = 0
   integer(int64) :: byte_limit = 0, reserved_bytes, fixed_bytes, entry_bytes
   integer(int64) :: pairs, hits, misses, bypasses, not_stored, probes
   logical :: active = .false., report_open = .false.
   logical, protected :: reuse_requested = .false.
contains
   subroutine configure_angular_reuse(memory_mib, filename)
      integer, intent(in) :: memory_mib
      character(len=*), intent(in) :: filename
      integer :: status
      byte_limit = int(memory_mib,int64)*1048576_int64
      reuse_requested = memory_mib > 0
      if (len_trim(filename) == 0) return
      if (.not. reuse_requested) error stop "--reuse-stats requires --reuse-cache-mib > 0"
      open(newunit=report_unit, file=filename, status="new", action="write", iostat=status)
      if (status /= 0) error stop "Cannot create reuse statistics file (it must not already exist)"
      report_open = .true.
      write(report_unit,"(A)") "block,pairs,hits,misses,bypasses,not_stored,entries,probes," // &
         "cache_reserved_bytes,cache_used_bytes,cache_limit_bytes"
   end subroutine configure_angular_reuse

   ! All arithmetic stays below 2**48; no signed integer overflow is required.
   ! Hashes select buckets only. Both the complete key and coupling arrays are
   ! compared on every hit, including a collision in either coupling hash.
   pure function mix(hash, word) result(next)
      integer(int64), intent(in) :: hash, word
      integer(int64) :: next
      next = iand(hash*65599_int64 + ibits(word,0,32), 2147483647_int64)
      next = iand(next*65599_int64 + ibits(word,32,32), 2147483647_int64)
   end function mix

   subroutine begin_reuse_block()
      integer :: entry_capacity, table_capacity, key_capacity, term_capacity, status, j, i
      integer(int64) :: available, remaining, hash
      type(pair_entry) :: sample
      pairs = 0; hits = 0; misses = 0; bypasses = 0; not_stored = 0; probes = 0
      entry_count = 0; key_count = 0; term_count = 0
      reserved_bytes = 0; fixed_bytes = 0
      entry_bytes = storage_size(sample,kind=int64)/8
      active = .false.
      if (.not. reuse_requested) return
      ! Debug printouts require execution of the original kernels for each pair.
      if (ibug1 /= 0 .or. ibug2 /= 0 .or. ibug3 /= 0 .or. any(ldbpa) .or. any(ldbpg)) return
      available = byte_limit - 8_int64*csf_count
      if (available < 1024) return
      entry_capacity = int(min(available/512, 1048576_int64))
      table_capacity = 2
      do while (table_capacity < 2*entry_capacity)
         table_capacity = 2*table_capacity
      end do
      fixed_bytes = 8_int64*csf_count + 4_int64*table_capacity
      remaining = byte_limit - fixed_bytes - entry_bytes*entry_capacity
      key_capacity = int(remaining/2/8)
      term_capacity = int((remaining - 8_int64*key_capacity)/28)
      allocate(coupling_hash(csf_count), buckets(table_capacity), entries(entry_capacity), &
         keys(key_capacity), term_labels(5,term_capacity), values(term_capacity), stat=status)
      if (status /= 0) then
         call release_cache()
         fixed_bytes = 0
         return
      end if
      reserved_bytes = fixed_bytes + entry_bytes*entry_capacity + 8_int64*key_capacity + 28_int64*term_capacity
      buckets = 0
      do j = 1, csf_count
         hash = 1
         do i = 1, max_orbitals
            hash = mix(hash, int(csf_coupling_quantum_codes(i,j),int64))
         end do
         coupling_hash(j) = hash
      end do
      active = .true.
   end subroutine begin_reuse_block

   pure function shell_word(orbital, csf) result(word)
      integer, intent(in) :: orbital, csf
      integer(int64) :: word
      integer :: i
      word = iand(int(csf_subshell_occupations(orbital,csf),int64),255_int64)
      do i = 1, 3
         word = ior(word,shiftl(iand(int(csf_subshell_quantum_codes(orbital,i,csf),int64),255_int64),8*i))
      end do
   end function shell_word

   subroutine pair_key(ja, jb, key, length, orbital_map, inverse_map, hash)
      integer, intent(in) :: ja, jb
      integer(int64), intent(out) :: key(3+2*max_orbitals), hash
      integer, intent(out) :: length, orbital_map(0:max_orbitals), inverse_map(max_orbitals)
      integer :: orbital, slot, i
      key(1:3) = [coupling_hash(ja), coupling_hash(jb), int(merge(1,0,ja==jb),int64)]
      slot = 0
      orbital_map = 0
      inverse_map = 0
      length = 3
      do orbital = 1, orbital_count
         ! Removing only jointly empty shells preserves the union's ordering,
         ! bra/ket equality pattern, fermion phases and every coupling chain.
         if (csf_subshell_occupations(orbital,ja) == 0 .and. csf_subshell_occupations(orbital,jb) == 0) cycle
         slot = slot + 1
         orbital_map(slot) = orbital
         inverse_map(orbital) = slot
         key(length+1) = ior(shiftl(int(orbital_kappa(orbital),int64),32),shell_word(orbital,ja))
         key(length+2) = shell_word(orbital,jb)
         length = length + 2
      end do
      hash = 1
      do i = 1, length
         hash = mix(hash,key(i))
      end do
      hash = ieor(hash,shiftr(hash,16))
   end subroutine pair_key

   subroutine original_pair(ja, jb, first, second, scalar)
      integer, intent(in) :: ja, jb
      integer, intent(out) :: first, second
      real(real64), intent(out) :: scalar
      real(real64) :: tshell(max_orbitals)
      first = 0; second = 0; tshell = 0
      if (ja /= jb) call onescalar(ja,jb,first,second,tshell)
      ! ONESCALAR leaves its other outputs undefined on some early returns.
      scalar = 0
      if (first /= 0) then
         if (first /= second) scalar = tshell(1)
      else
         second = 0
      end if
      angular_coefficient_count = 0
      call rkco_gg(ja,jb,cord,0,1)
   end subroutine original_pair

   subroutine load_reused_pair(ja, jb, first, second, scalar)
      integer, intent(in) :: ja, jb
      integer, intent(out) :: first, second
      real(real64), intent(out) :: scalar
      integer(int64) :: key(3+2*max_orbitals), hash
      integer :: length, orbital_map(0:max_orbitals), inverse_map(max_orbitals)
      integer :: bucket, entry, start, count, i, j
      pairs = pairs + 1
      if (.not. active) then
         bypasses = bypasses + 1
         call original_pair(ja,jb,first,second,scalar)
         return
      end if
      call pair_key(ja,jb,key,length,orbital_map,inverse_map,hash)
      bucket = int(iand(hash,int(size(buckets)-1,int64))) + 1
      do
         entry = buckets(bucket)
         if (entry == 0) exit
         probes = probes + 1
         if (entries(entry)%hash == hash .and. entries(entry)%key_size == length) then
            start = entries(entry)%key_start
            if (all(keys(start:start+length-1) == key(1:length))) then
               if (all(csf_coupling_quantum_codes(:,ja) == csf_coupling_quantum_codes(:,entries(entry)%bra)) .and. &
                   all(csf_coupling_quantum_codes(:,jb) == csf_coupling_quantum_codes(:,entries(entry)%ket))) then
                  hits = hits + 1
                  first = orbital_map(entries(entry)%first)
                  second = orbital_map(entries(entry)%second)
                  scalar = entries(entry)%scalar
                  angular_coefficient_count = entries(entry)%term_count
                  do while (angular_buffer_capacity < angular_coefficient_count)
                     call alcbuf(2)
                  end do
                  start = entries(entry)%term_start
                  do i = 1, angular_coefficient_count
                     do j = 1, 4
                        angular_coefficient_labels(j,i) = orbital_map(term_labels(j,start+i-1))
                     end do
                     angular_coefficient_labels(5,i) = term_labels(5,start+i-1)
                     angular_coefficients(i) = values(start+i-1)
                  end do
                  return
               end if
            end if
         end if
         bucket = iand(bucket,size(buckets)-1) + 1
      end do
      misses = misses + 1
      call original_pair(ja,jb,first,second,scalar)
      count = angular_coefficient_count
      if (entry_count == size(entries) .or. key_count+length > size(keys) .or. &
          int(term_count,int64)+count > size(values,kind=int64)) then
         not_stored = not_stored + 1
         return
      end if
      ! A future kernel producing an orbital outside the occupied union falls
      ! back without caching. Never infer an extra mapping or recoupling phase.
      if (first < 0 .or. first > orbital_count .or. second < 0 .or. second > orbital_count) then
         not_stored = not_stored + 1
         return
      end if
      if (first > 0) then
         if (inverse_map(first) == 0) then
            not_stored = not_stored + 1
            return
         end if
      end if
      if (second > 0) then
         if (inverse_map(second) == 0) then
            not_stored = not_stored + 1
            return
         end if
      end if
      do i = 1, count
         do j = 1, 4
            if (angular_coefficient_labels(j,i) < 1 .or. angular_coefficient_labels(j,i) > orbital_count) then
               not_stored = not_stored + 1
               return
            end if
            if (inverse_map(angular_coefficient_labels(j,i)) == 0) then
               not_stored = not_stored + 1
               return
            end if
         end do
      end do
      entry_count = entry_count + 1
      entry = entry_count
      buckets(bucket) = entry
      entries(entry)%hash = hash
      entries(entry)%key_start = key_count+1
      entries(entry)%key_size = length
      entries(entry)%term_start = term_count+1
      entries(entry)%term_count = count
      entries(entry)%bra = ja
      entries(entry)%ket = jb
      entries(entry)%first = 0
      entries(entry)%second = 0
      if (first > 0) entries(entry)%first = inverse_map(first)
      if (second > 0) entries(entry)%second = inverse_map(second)
      entries(entry)%scalar = scalar
      keys(key_count+1:key_count+length) = key(1:length)
      key_count = key_count + length
      do i = 1, count
         do j = 1, 4
            term_labels(j,term_count+i) = inverse_map(angular_coefficient_labels(j,i))
         end do
         term_labels(5,term_count+i) = angular_coefficient_labels(5,i)
         values(term_count+i) = angular_coefficients(i)
      end do
      term_count = term_count + count
   end subroutine load_reused_pair

   subroutine release_cache()
      if (allocated(entries)) deallocate(entries)
      if (allocated(buckets)) deallocate(buckets)
      if (allocated(keys)) deallocate(keys)
      if (allocated(coupling_hash)) deallocate(coupling_hash)
      if (allocated(term_labels)) deallocate(term_labels)
      if (allocated(values)) deallocate(values)
      active = .false.
   end subroutine release_cache

   subroutine end_reuse_block(block_index)
      integer, intent(in) :: block_index
      integer(int64) :: used_bytes
      used_bytes = fixed_bytes + entry_bytes*entry_count + 8_int64*key_count + 28_int64*term_count
      print "(A,I0,A,I0,A,I0,A,I0)", "Angular reuse block ",block_index,": hits=",hits, &
         ", pairs=",pairs,", cache reserved bytes=",reserved_bytes
      if (report_open) then
         write(report_unit,'(I0,10(",",I0))') block_index,pairs,hits,misses,bypasses,not_stored, &
            entry_count,probes,reserved_bytes,used_bytes,byte_limit
         flush(report_unit)
      end if
      call release_cache()
   end subroutine end_reuse_block
end module rangular_reuse
