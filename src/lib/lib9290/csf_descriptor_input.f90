! Versioned descriptor transport and ordinary-CSF semantic adapter.
! The physical state table remains owned by terms_C; no table is duplicated.
module csf_descriptor_input
   use iso_fortran_env, only: int8, int64, error_unit, file_storage_size
   use parameter_def, only: orbital_capacity => NNNW
   use terms_C, only: state_counts => ITAB, state_offsets => JTAB, state_values => NTAB
   implicit none
   private
   public :: descriptor_source, open_descriptor_source, close_descriptor_source
   public :: read_descriptor_block, install_descriptor_header, descriptor_block_label

   type :: descriptor_source
      integer :: unit_number = 0
      logical :: is_open = .false.
      integer :: core_count = 0, peel_count = 0, orbital_count = 0, block_count = 0
      integer :: row_width = 0, total_csf_count = 0, electron_count = 0, next_block_index = 1
      integer, allocatable :: principal_numbers(:), relativistic_kappas(:)
      integer, allocatable :: orbital_momenta(:), single_electron_two_j(:)
      character(len=2), allocatable :: orbital_labels(:)
      integer, allocatable :: block_lengths(:), block_two_j(:), block_parities(:)
      character(len=64) :: header_sha256 = ''
   end type descriptor_source

contains

   subroutine require_descriptor(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(error_unit, '(A)') 'CSF descriptor input: '//message
         error stop 1
      end if
   end subroutine require_descriptor

   function read_unsigned_integer(unit_number, byte_count) result(value)
      integer, intent(in) :: unit_number, byte_count
      integer(int64) :: value
      integer(int8) :: encoded_bytes(8)
      integer :: byte_index, read_status
      read(unit_number, iostat=read_status) encoded_bytes(:byte_count)
      call require_descriptor(read_status == 0, 'truncated binary header')
      value = 0_int64
      do byte_index = 1, byte_count
         value = ior(value, ishft(iand(int(encoded_bytes(byte_index), int64), 255_int64), 8*(byte_index-1)))
      end do
      call require_descriptor(value >= 0_int64, 'unsigned count exceeds signed int64 range')
   end function read_unsigned_integer

   function read_default_integer(unit_number, signed_value) result(value)
      integer, intent(in) :: unit_number
      logical, intent(in), optional :: signed_value
      integer :: value
      integer(int64) :: wide_value
      wide_value = read_unsigned_integer(unit_number, 4)
      if (present(signed_value)) then
         if (signed_value .and. wide_value >= 2147483648_int64) wide_value = wide_value - 4294967296_int64
      end if
      call require_descriptor(abs(wide_value) <= int(huge(value), int64), 'integer exceeds suite default-integer range')
      value = int(wide_value)
   end function read_default_integer

   subroutine open_descriptor_source(source, filename, maximum_blocks)
      type(descriptor_source), intent(inout) :: source
      character(len=*), intent(in) :: filename
      integer, intent(in) :: maximum_blocks
      character(len=8) :: file_magic
      character(len=*), parameter :: orbital_letters = 'spdfghiklm'
      integer :: read_status, orbital_index, previous_index, block_index, file_version, byte_marker
      integer :: allocation_status, relativistic_kappa, orbital_momentum, hash_index
      integer(int64) :: total_rows, row_count, counted_rows, file_bytes, expected_bytes

      call close_descriptor_source(source)
      source = descriptor_source()
      call require_descriptor(file_storage_size == 8 .and. storage_size(0_int8) == 8, &
                              'this transport requires eight-bit file storage units and int8')
      open(newunit=source%unit_number, file=filename, access='stream', form='unformatted', &
           status='old', action='read', iostat=read_status)
      call require_descriptor(read_status == 0, 'cannot open '//filename)
      source%is_open = .true.
      inquire(unit=source%unit_number, size=file_bytes)
      read(source%unit_number, iostat=read_status) file_magic
      call require_descriptor(read_status == 0, 'missing binary magic')
      call require_descriptor(file_magic == 'GCSDESC1', 'unrecognized binary magic')
      file_version = read_default_integer(source%unit_number)
      byte_marker = read_default_integer(source%unit_number)
      call require_descriptor(file_version == 1, 'unsupported binary version')
      call require_descriptor(byte_marker == 16909060, 'invalid little-endian marker')
      source%core_count = read_default_integer(source%unit_number)
      source%peel_count = read_default_integer(source%unit_number)
      call require_descriptor(source%core_count <= orbital_capacity .and. source%peel_count <= orbital_capacity, &
                              'orbital table exceeds configured capacity')
      source%orbital_count = source%core_count + source%peel_count
      call require_descriptor(source%orbital_count > 0 .and. source%orbital_count <= orbital_capacity, &
                              'orbital table exceeds configured capacity or is empty')
      source%block_count = read_default_integer(source%unit_number)
      call require_descriptor(source%block_count > 0 .and. source%block_count <= maximum_blocks, &
                              'block count exceeds caller capacity or is empty')
      source%row_width = read_default_integer(source%unit_number)
      call require_descriptor(source%row_width == 4*source%peel_count + 2, 'incorrect descriptor row width')
      total_rows = read_unsigned_integer(source%unit_number, 8)
      call require_descriptor(total_rows > 0 .and. total_rows <= int(huge(source%total_csf_count), int64), &
                              'CSF count exceeds current suite integer range or is empty')
      source%total_csf_count = int(total_rows)
      source%electron_count = read_default_integer(source%unit_number)
      call require_descriptor(source%electron_count > 0, 'electron count must be positive')
      read(source%unit_number, iostat=read_status) source%header_sha256
      call require_descriptor(read_status == 0, 'missing source header digest')
      do hash_index = 1, len(source%header_sha256)
         call require_descriptor(index('0123456789abcdef', source%header_sha256(hash_index:hash_index)) > 0, &
                                 'invalid source header digest')
      end do
      expected_bytes = 108_int64 + 8_int64*source%orbital_count + 16_int64*source%block_count &
                       + int(source%row_width, int64)*total_rows
      call require_descriptor(file_bytes == expected_bytes, 'binary length mismatch (truncation or trailing data)')
      allocate(source%principal_numbers(source%orbital_count), source%relativistic_kappas(source%orbital_count), &
               source%orbital_momenta(source%orbital_count), source%single_electron_two_j(source%orbital_count), &
               source%orbital_labels(source%orbital_count), source%block_lengths(source%block_count), &
               source%block_two_j(source%block_count), source%block_parities(source%block_count), stat=allocation_status)
      call require_descriptor(allocation_status == 0, 'cannot allocate header arrays')
      do orbital_index = 1, source%orbital_count
         source%principal_numbers(orbital_index) = read_default_integer(source%unit_number)
         relativistic_kappa = read_default_integer(source%unit_number, .true.)
         call require_descriptor(relativistic_kappa /= 0 .and. abs(relativistic_kappa) <= 10, &
                                 'unsupported relativistic kappa')
         orbital_momentum = relativistic_kappa
         if (relativistic_kappa < 0) orbital_momentum = -relativistic_kappa - 1
         call require_descriptor(orbital_momentum < len(orbital_letters), 'unsupported orbital letter')
         call require_descriptor(source%principal_numbers(orbital_index) > orbital_momentum, &
                                 'invalid principal quantum number')
         source%relativistic_kappas(orbital_index) = relativistic_kappa
         source%orbital_momenta(orbital_index) = orbital_momentum
         source%single_electron_two_j(orbital_index) = 2*abs(relativistic_kappa) - 1
         source%orbital_labels(orbital_index) = orbital_letters(orbital_momentum+1:orbital_momentum+1)//' '
         if (relativistic_kappa > 0) source%orbital_labels(orbital_index)(2:2) = '-'
         do previous_index = 1, orbital_index - 1
            call require_descriptor(source%principal_numbers(previous_index) /= source%principal_numbers(orbital_index) &
                                    .or. source%relativistic_kappas(previous_index) /= relativistic_kappa, &
                                    'duplicate orbital in core/peel table')
         end do
      end do
      counted_rows = 0_int64
      do block_index = 1, source%block_count
         row_count = read_unsigned_integer(source%unit_number, 8)
         call require_descriptor(row_count > 0 .and. row_count <= total_rows - counted_rows, 'invalid block length')
         counted_rows = counted_rows + row_count
         source%block_lengths(block_index) = int(row_count)
         source%block_two_j(block_index) = read_default_integer(source%unit_number)
         source%block_parities(block_index) = read_default_integer(source%unit_number, .true.)
         call require_descriptor(source%block_two_j(block_index) <= 126, 'block angular momentum exceeds byte limit')
         call require_descriptor(abs(source%block_parities(block_index)) == 1, 'invalid block parity')
      end do
      call require_descriptor(counted_rows == total_rows, 'block lengths do not cover all rows')
   end subroutine open_descriptor_source

   subroutine close_descriptor_source(source)
      type(descriptor_source), intent(inout) :: source
      if (source%is_open) close(source%unit_number)
      source%is_open = .false.
   end subroutine close_descriptor_source

   function descriptor_block_label(total_two_j, parity) result(label)
      integer, intent(in) :: total_two_j, parity
      character(len=8) :: label
      character(len=16) :: angular_text
      character(len=20) :: padded_text
      character(len=1) :: parity_symbol
      if (mod(total_two_j, 2) == 0) then
         write(angular_text, '(I0)') total_two_j/2
      else
         write(angular_text, '(I0,A)') total_two_j, '/2'
      end if
      parity_symbol = '+'
      if (parity < 0) parity_symbol = '-'
      padded_text = trim(angular_text)//parity_symbol
      padded_text = adjustr(padded_text)
      ! SETCSLL takes the final five characters, followed by three blanks.
      label = padded_text(16:20)
   end function descriptor_block_label

   subroutine install_descriptor_header(source, core_count, block_labels)
      use orb_C, only: current_csf_count => NCF, current_orbital_count => NW, &
                       principal_numbers => NP, relativistic_kappas => NAK, orbital_momenta => NKL, &
                       single_electron_two_j => NKJ, orbital_labels => NH
      use def_C, only: electron_count => NELEC
      use hblock_C, only: block_count => NBLOCK, block_lengths => NCFBLK
      type(descriptor_source), intent(in) :: source
      integer, intent(out) :: core_count
      character(len=8), intent(out) :: block_labels(:)
      integer :: block_index, allocation_status
      call require_descriptor(size(block_labels) >= source%block_count, 'insufficient block label capacity')
      core_count = source%core_count
      current_orbital_count = source%orbital_count
      current_csf_count = source%total_csf_count
      electron_count = source%electron_count
      block_count = source%block_count
      principal_numbers(:current_orbital_count) = source%principal_numbers
      relativistic_kappas(:current_orbital_count) = source%relativistic_kappas
      orbital_momenta(:current_orbital_count) = source%orbital_momenta
      single_electron_two_j(:current_orbital_count) = source%single_electron_two_j
      orbital_labels(:current_orbital_count) = source%orbital_labels
      allocate(block_lengths(block_count+1), stat=allocation_status)
      call require_descriptor(allocation_status == 0, 'cannot allocate block lengths')
      block_lengths = 0
      block_lengths(:block_count) = source%block_lengths
      block_labels = ''
      do block_index = 1, block_count
         block_labels(block_index) = descriptor_block_label(source%block_two_j(block_index), &
                                                           source%block_parities(block_index))
      end do
   end subroutine install_descriptor_header

   function resolve_subshell_state(single_two_j, occupation, physical_seniority, printed_two_j, block_index, row_index) &
      result(subshell_two_j)
      integer, intent(in) :: single_two_j, occupation, physical_seniority, printed_two_j
      integer, intent(in) :: block_index, row_index
      integer :: subshell_two_j
      integer :: capacity, population, table_row, table_entry, table_offset, match_count, matching_j_count
      integer :: candidate_two_j, candidate_seniority, candidate_label
      capacity = single_two_j + 1
      population = min(occupation, capacity - occupation)
      if (population == 0) then
         call require_state(physical_seniority == 0 .and. (printed_two_j == -1 .or. printed_two_j == 0), &
                                 'empty/closed shell must have seniority zero and J=0 or omitted J')
         subshell_two_j = 0
         return
      end if
      call require_state(single_two_j <= 17, 'unsupported open-shell state table for single-electron 2j > 17')
      if (single_two_j > 9) call require_state(population <= 2, 'unsupported high-j occupation')
      table_row = (capacity-2)/2
      table_row = table_row*(table_row+1)/2 + population + 1
      call require_state(table_row >= 1 .and. table_row <= size(state_counts), 'unsupported subshell table row')
      match_count = 0
      subshell_two_j = -1
      do table_entry = 1, state_counts(table_row)
         table_offset = state_offsets(table_row) + 3*(table_entry-1)
         candidate_seniority = state_values(table_offset+1)
         candidate_label = state_values(table_offset+2)
         candidate_two_j = state_values(table_offset+3) - 1
         if (candidate_seniority /= physical_seniority) cycle
         if (printed_two_j /= -1 .and. printed_two_j /= candidate_two_j) cycle
         match_count = match_count + 1
         call require_state(candidate_label == 0, 'unsupported state requires an additional w label')
         subshell_two_j = candidate_two_j
      end do
      call require_state(match_count == 1, 'invalid or ambiguous subshell (occupation, seniority, J)')
      ! Ordinary LODCSH2 only distinguishes repeated J by v for j=7/2,q=4.
      matching_j_count = 0
      do table_entry = 1, state_counts(table_row)
         table_offset = state_offsets(table_row) + 3*(table_entry-1)
         if (state_values(table_offset+3) == subshell_two_j+1) matching_j_count = matching_j_count + 1
      end do
      call require_state(matching_j_count == 1 .or. (single_two_j == 7 .and. occupation == 4), &
                              'unsupported state is ambiguous in the ordinary text reader')
   contains
      subroutine require_state(condition, message)
         logical, intent(in) :: condition
         character(len=*), intent(in) :: message
         if (.not. condition) then
            write(error_unit, '(A,I0,A,I0,2A)') 'CSF descriptor input: block ', block_index, &
                                                ', CSF ', row_index, ': ', message
            error stop 1
         end if
      end subroutine require_state
   end function resolve_subshell_state

   subroutine decode_descriptor_row(source, block_index, row_index, row, occupations, subshell_quantum_numbers, couplings)
      type(descriptor_source), intent(in) :: source
      integer, intent(in) :: block_index, row_index
      integer(int8), intent(in) :: row(:)
      integer(int8), intent(out) :: occupations(:), subshell_quantum_numbers(:,:), couplings(:)
      integer :: shell_two_j(orbital_capacity), printed_two_k(orbital_capacity)
      integer :: peel_index, orbital_index, channel_offset, occupation, seniority, printed_two_j
      integer :: occupied_count, occupied_position, open_count, capacity, previous_two_j, coupled_two_j
      integer :: last_nonzero_shell, parity_exponent, total_two_j, parity, electron_count

      occupations = 0_int8
      subshell_quantum_numbers = 0_int8
      couplings = 0_int8
      shell_two_j = 0
      printed_two_k = -1
      electron_count = 0
      parity_exponent = 0
      last_nonzero_shell = 0
      do orbital_index = 1, source%core_count
         capacity = source%single_electron_two_j(orbital_index) + 1
         occupations(orbital_index) = int(capacity, int8)
         subshell_quantum_numbers(orbital_index,3) = 1_int8
         electron_count = electron_count + capacity
      end do
      do peel_index = 1, source%peel_count
         orbital_index = source%core_count + peel_index
         channel_offset = 4*(peel_index-1)
         occupation = int(row(channel_offset+1))
         printed_two_j = int(row(channel_offset+2))
         seniority = int(row(channel_offset+3))
         printed_two_k(orbital_index) = int(row(channel_offset+4))
         capacity = source%single_electron_two_j(orbital_index) + 1
         call require_row(occupation >= 0 .and. occupation <= capacity, 'occupation outside shell capacity')
         call require_row(seniority >= 0, 'negative physical seniority')
         call require_row(printed_two_j >= -1 .and. printed_two_j <= 126, 'invalid subshell 2J')
         call require_row(printed_two_k(orbital_index) >= -1 .and. printed_two_k(orbital_index) <= 126, &
                                 'invalid intermediate 2K')
         if (occupation == 0) then
            call require_row(printed_two_j == -1 .and. seniority == 0 .and. &
                                    printed_two_k(orbital_index) == -1, 'noncanonical empty shell descriptor')
         end if
         shell_two_j(orbital_index) = resolve_subshell_state(source%single_electron_two_j(orbital_index), &
                                                            occupation, seniority, printed_two_j, block_index, row_index)
         occupations(orbital_index) = int(occupation, int8)
         subshell_quantum_numbers(orbital_index,1) = int(seniority, int8)
         subshell_quantum_numbers(orbital_index,3) = int(shell_two_j(orbital_index)+1, int8)
         if (shell_two_j(orbital_index) /= 0) last_nonzero_shell = orbital_index
         electron_count = electron_count + occupation
         parity_exponent = parity_exponent + occupation*source%orbital_momenta(orbital_index)
      end do
      total_two_j = int(row(source%row_width-1))
      parity = int(row(source%row_width))
      call require_row(total_two_j == source%block_two_j(block_index) .and. &
                              parity == source%block_parities(block_index), 'row/block symmetry mismatch')
      call require_row(parity == (-1)**parity_exponent, 'incorrect row parity')
      call require_row(electron_count == source%electron_count, 'incorrect row electron count')
      occupied_count = count(occupations(source%core_count+1:source%orbital_count) > 0)
      occupied_position = 0
      open_count = 0
      previous_two_j = 0
      do orbital_index = source%core_count+1, source%orbital_count
         occupation = int(occupations(orbital_index))
         if (occupation == 0) cycle
         occupied_position = occupied_position + 1
         capacity = source%single_electron_two_j(orbital_index) + 1
         coupled_two_j = printed_two_k(orbital_index)
         if (coupled_two_j >= 0) then
            call require_row(occupied_position >= 2 .and. occupied_position < occupied_count, &
                                    'intermediate coupling at invalid occupied-shell boundary')
         else if (orbital_index >= last_nonzero_shell) then
            coupled_two_j = total_two_j
         else if (shell_two_j(orbital_index) == 0) then
            coupled_two_j = previous_two_j
         else if (previous_two_j == 0) then
            coupled_two_j = shell_two_j(orbital_index)
         else
            call require_row(.false., 'missing nontrivial intermediate coupling')
         end if
         call require_row(coupled_two_j >= abs(previous_two_j-shell_two_j(orbital_index)) .and. &
                                 coupled_two_j <= previous_two_j+shell_two_j(orbital_index) .and. &
                                 mod(previous_two_j+shell_two_j(orbital_index)-coupled_two_j, 2) == 0, &
                                 'coupling violates angular momentum triangle/parity')
         previous_two_j = coupled_two_j
         if (occupation < capacity) then
            open_count = open_count + 1
            if (open_count > 1) couplings(open_count-1) = int(coupled_two_j+1, int8)
         end if
      end do
      call require_row(previous_two_j == total_two_j, 'incorrect final angular momentum')
      couplings(orbital_capacity) = int((total_two_j+1)*parity, int8)
   contains
      subroutine require_row(condition, message)
         logical, intent(in) :: condition
         character(len=*), intent(in) :: message
         if (.not. condition) then
            write(error_unit, '(A,I0,A,I0,2A)') 'CSF descriptor input: block ', block_index, &
                                                ', CSF ', row_index, ': ', message
            error stop 1
         end if
      end subroutine require_row
   end subroutine decode_descriptor_row

   subroutine read_descriptor_block(source, block_index, occupations, subshell_quantum_numbers, couplings)
      type(descriptor_source), intent(inout) :: source
      integer, intent(in) :: block_index
      integer(int8), intent(out) :: occupations(:,:), subshell_quantum_numbers(:,:,:), couplings(:,:)
      integer(int8), allocatable :: batch_values(:,:)
      integer :: batch_capacity, batch_count, row_offset, batch_row, read_status, allocation_status
      call require_descriptor(source%is_open, 'descriptor source is not open')
      call require_descriptor(block_index == source%next_block_index .and. block_index <= source%block_count, &
                              'blocks must be loaded once in original order')
      call require_descriptor(size(occupations,1) == orbital_capacity .and. &
                              size(occupations,2) == source%block_lengths(block_index), 'occupation array shape mismatch')
      call require_descriptor(all(shape(couplings) == shape(occupations)), 'coupling array shape mismatch')
      call require_descriptor(size(subshell_quantum_numbers,1) == orbital_capacity .and. &
                              size(subshell_quantum_numbers,2) == 3 .and. &
                              size(subshell_quantum_numbers,3) == size(occupations,2), 'quantum number array shape mismatch')
      batch_capacity = min(8192, source%block_lengths(block_index))
      allocate(batch_values(source%row_width,batch_capacity), stat=allocation_status)
      call require_descriptor(allocation_status == 0, 'cannot allocate descriptor batch')
      row_offset = 0
      do while (row_offset < source%block_lengths(block_index))
         batch_count = min(batch_capacity, source%block_lengths(block_index)-row_offset)
         read(source%unit_number, iostat=read_status) batch_values(:,:batch_count)
         call require_descriptor(read_status == 0, 'truncated descriptor payload')
         do batch_row = 1, batch_count
            call decode_descriptor_row(source, block_index, row_offset+batch_row, batch_values(:,batch_row), &
                                       occupations(:,row_offset+batch_row), &
                                       subshell_quantum_numbers(:,:,row_offset+batch_row), couplings(:,row_offset+batch_row))
         end do
         row_offset = row_offset + batch_count
      end do
      source%next_block_index = source%next_block_index + 1
   end subroutine read_descriptor_block
end module csf_descriptor_input
