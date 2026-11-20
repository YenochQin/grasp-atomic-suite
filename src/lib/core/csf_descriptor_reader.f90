! Read and validate CSF descriptor headers and blocks, then decode the rows into Fortran basis arrays.
module csf_descriptor_reader
    ! Module dependencies
    use iso_fortran_env, only: int8, int64, error_unit
    use iso_c_binding, only: c_ptr, c_null_ptr, c_null_char, c_char, c_size_t, c_int32_t, c_int64_t
    use csf_descriptor_c_bindings
    use grasp_parameters, only: orbital_capacity => max_orbitals
    use subshell_term_table, only: state_counts => itab, state_offsets => jtab, state_values => ntab
    implicit none (type, external)

    ! Module interface
    private
    public :: descriptor_source, open_descriptor_source, close_descriptor_source
    public :: read_descriptor_block, install_descriptor_header, descriptor_block_label
    public :: read_descriptor_zero_order_counts

    ! Derived types
    ! Owns a forward-only Rust stream. Close the source explicitly and avoid copying a live handle.
    type :: descriptor_source
        ! Components
        type(c_ptr) :: handle = c_null_ptr
        logical :: is_open = .false., header_loaded = .false.
        integer :: core_count = 0, peel_count = 0, orbital_count = 0, block_count = 0
        integer :: row_width = 0, total_csf_count = 0, electron_count = 0, next_block_index = 1
        integer, allocatable :: principal_numbers(:), relativistic_kappas(:)
        integer, allocatable :: orbital_momenta(:), single_electron_two_j(:)
        character(len=2), allocatable :: orbital_labels(:)
        integer, allocatable :: block_lengths(:), block_two_j(:), block_parities(:)
    end type descriptor_source

contains

    logical function failed(status)
        ! Arguments
        integer, intent(in), optional :: status
        failed = .false.
        if (present(status)) failed = status /= 0
    end function failed

    ! Return a failure status when requested. Otherwise report the error and stop.
    logical function require_descriptor(condition, message, status) result(ok)
        ! Arguments
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        integer, intent(out), optional :: status
        ok = condition
        if (present(status)) status = 0
        if (.not. condition) then
            write(error_unit, "(A)") "CSF descriptor input: "//message
            if (present(status)) then
                status = 1
            else
                error stop 1
            end if
        end if
    end function require_descriptor

    subroutine open_descriptor_source(source, filename, maximum_blocks, header_filename, status)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        character(len=*), intent(in) :: filename
        integer, intent(in) :: maximum_blocks
        character(len=*), intent(in), optional :: header_filename
        integer, intent(out), optional :: status

        ! Local declarations
        integer :: open_status
        call open_descriptor_source_impl(source, filename, maximum_blocks, header_filename, open_status)
        if (open_status /= 0) then
            call close_descriptor_source(source)
            source = descriptor_source()
        end if
        if (present(status)) then
            status = open_status
        else if (open_status /= 0) then
            error stop "Cannot open descriptor source"
        end if
    end subroutine open_descriptor_source

    subroutine open_descriptor_source_impl(source, filename, maximum_blocks, header_filename, status)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        character(len=*), intent(in) :: filename
        integer, intent(in) :: maximum_blocks
        character(len=*), intent(in), optional :: header_filename
        integer, intent(out), optional :: status

        ! Local declarations
        character(len=:), allocatable :: header_path, descriptor_stem
        character(kind=c_char,len=1024) :: reader_message
        integer(c_int64_t) :: header_counts(5)
        integer(c_int64_t), allocatable :: block_lengths(:)
        integer(c_int32_t), allocatable :: orbital_table(:,:), symmetry_table(:,:)
        integer :: read_status, orbital_index, kappa, orbital_l_value, path_length, allocation_status
        character(len=*), parameter :: orbital_letters = "spdfghiklm"
        reader_message = c_null_char
        call close_descriptor_source(source)
        source = descriptor_source()
        if (.not. require_descriptor(rust_descriptor_version() == 1, "unsupported Rust ABI", status)) return
        path_length = len_trim(filename)
        if (.not. require_descriptor(path_length > 8, "expected descriptor Parquet path", status)) return
        if (.not. require_descriptor(filename(path_length-7:path_length) == ".parquet", &
            & "only descriptor Parquet input is supported", status)) return
        descriptor_stem = filename(:path_length-8)
        path_length = len(descriptor_stem)
        if (path_length > 12) then
            if (descriptor_stem(path_length-11:) == "_descriptors") descriptor_stem = descriptor_stem(:path_length-12)
        end if
        ! Derive the paired header path unless the caller supplied one.
        header_path = descriptor_stem//"_header.toml"
        if (present(header_filename)) then
            if (len_trim(header_filename) > 0) header_path = header_filename
        end if
        read_status = rust_descriptor_open(trim(filename)//c_null_char, header_path//c_null_char, &
            source%handle, reader_message, 1024_c_size_t)
        if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), &
            & status)) return
        source%is_open = .true.
        read_status = rust_descriptor_info(source%handle, header_counts, reader_message, 1024_c_size_t)
        if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), &
            & status)) return
        if (.not. require_descriptor(all(header_counts >= 0) .and. all(header_counts <= int(huge(0),int64)), &
            & "header integer overflow",  &
            & status)) return
        source%core_count = int(header_counts(1))
        source%peel_count = int(header_counts(2))
        if (.not. require_descriptor(header_counts(1)+header_counts(2) <= orbital_capacity, "orbital capacity exceeded", &
            & status)) return
        source%orbital_count = int(header_counts(1)+header_counts(2))
        source%block_count = int(header_counts(3))
        source%total_csf_count = int(header_counts(4))
        source%electron_count = int(header_counts(5))
        source%row_width = 4*source%peel_count+2
        if (.not. require_descriptor(source%block_count > 0 .and. source%block_count < huge(0) .and. &
            source%block_count <= maximum_blocks, &
            "invalid block count", status)) return
        allocate(orbital_table(2,source%orbital_count), block_lengths(source%block_count), &
            symmetry_table(2,source%block_count), &
            source%principal_numbers(source%orbital_count), source%relativistic_kappas(source%orbital_count), &
            source%orbital_momenta(source%orbital_count), source%single_electron_two_j(source%orbital_count), &
            source%orbital_labels(source%orbital_count), source%block_lengths(source%block_count), &
            source%block_two_j(source%block_count), source%block_parities(source%block_count), stat=allocation_status)
        if (.not. require_descriptor(allocation_status == 0, "cannot allocate header arrays", status)) return
        ! The leading dimension of 2 matches the interleaved (n, kappa) and (2J, parity) ABI pairs.
        read_status = rust_descriptor_tables(source%handle, orbital_table, int(source%orbital_count,c_size_t), &
            block_lengths, symmetry_table, int(source%block_count,c_size_t), reader_message, 1024_c_size_t)
        if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), &
            & status)) return
        if (.not. require_descriptor(all(block_lengths > 0) .and. all(block_lengths <= huge(0)), "invalid block lengths", &
            & status)) return
        source%block_lengths = int(block_lengths)
        source%block_two_j = int(symmetry_table(1,:))
        source%block_parities = int(symmetry_table(2,:))
        source%principal_numbers = int(orbital_table(1,:))
        source%relativistic_kappas = int(orbital_table(2,:))
        do orbital_index = 1, source%orbital_count
            kappa = source%relativistic_kappas(orbital_index)
            orbital_l_value = kappa
            if (kappa < 0) orbital_l_value = -kappa-1
            if (.not. require_descriptor(kappa /= 0 .and. orbital_l_value >= 0 .and. orbital_l_value < len(orbital_letters), &
                & "unsupported orbital", status)) return
            source%orbital_momenta(orbital_index) = orbital_l_value
            source%single_electron_two_j(orbital_index) = 2*abs(kappa)-1
            source%orbital_labels(orbital_index) = orbital_letters(orbital_l_value+1:orbital_l_value+1)//" "
            if (kappa > 0) source%orbital_labels(orbital_index)(2:2) = "-"
        end do
        source%header_loaded = .true.
    end subroutine open_descriptor_source_impl

    subroutine read_descriptor_zero_order_counts(source, counts, status)
        ! Arguments
        type(descriptor_source), intent(in) :: source
        integer, intent(out) :: counts(:), status

        ! Local declarations
        integer(c_int64_t) :: raw_zero_order_counts(size(counts))
        character(kind=c_char,len=1024) :: reader_message
        integer :: read_status
        if (.not. require_descriptor(source%is_open .and. size(counts) == source%block_count, &
            "invalid source or buffer for zero-order counts", status)) return
        read_status = rust_descriptor_zero_order_counts(source%handle, raw_zero_order_counts, int(size(counts),c_size_t), &
            reader_message, 1024_c_size_t)
        if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), &
            & status)) return
        counts = int(raw_zero_order_counts)
    end subroutine read_descriptor_zero_order_counts

    subroutine close_descriptor_source(source)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        call rust_descriptor_close(source%handle)
        source%is_open = .false.
    end subroutine close_descriptor_source

    ! Right-align J and parity in the first five characters of the eight-character block label.
    function descriptor_block_label(total_two_j, parity) result(label)
        ! Arguments
        integer, intent(in) :: total_two_j, parity
        character(len=8) :: label

        ! Local declarations
        character(len=16) :: angular_text
        character(len=20) :: padded_text
        character(len=1) :: parity_symbol
        if (mod(total_two_j, 2) == 0) then
            write(angular_text, "(I0)") total_two_j/2
        else
            write(angular_text, "(I0,A)") total_two_j, "/2"
        end if
        parity_symbol = "+"
        if (parity < 0) parity_symbol = "-"
        padded_text = trim(angular_text)//parity_symbol
        padded_text = adjustr(padded_text)
        label = padded_text(16:20)
    end function descriptor_block_label

    subroutine install_descriptor_header(source, core_count, block_labels, status)
        ! Module dependencies
        use orbital_csf_basis, only: current_csf_count => csf_count, current_orbital_count => orbital_count, &
            principal_numbers => orbital_principal_n, relativistic_kappas => orbital_kappa, orbital_momenta &
            & => orbital_l, &
            single_electron_two_j => orbital_two_j, orbital_labels => orbital_symmetry_labels
        use atomic_calculation_state, only: electron_count => electron_count
        use hamiltonian_symmetry_blocks, only: block_count => symmetry_block_count, block_lengths => csf_counts_by_block

        ! Arguments
        type(descriptor_source), intent(in) :: source
        integer, intent(out) :: core_count
        character(len=8), intent(out) :: block_labels(:)
        integer, intent(out), optional :: status

        ! Local declarations
        integer :: block_index, allocation_status
        integer, pointer :: prepared_lengths(:)
        if (.not. require_descriptor(source%header_loaded, "descriptor header is not loaded", status)) return
        if (.not. require_descriptor(size(block_labels) >= source%block_count, "insufficient block label capacity", &
            & status)) return
        allocate(prepared_lengths(source%block_count+1), stat=allocation_status)
        if (.not. require_descriptor(allocation_status == 0, "cannot allocate block lengths", status)) return
        prepared_lengths=0
        prepared_lengths(:source%block_count)=source%block_lengths
        if (associated(block_lengths)) deallocate(block_lengths)
        block_lengths=>prepared_lengths
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
        block_lengths = 0
        block_lengths(:block_count) = source%block_lengths
        block_labels = ""
        do block_index = 1, block_count
            block_labels(block_index) = descriptor_block_label(source%block_two_j(block_index), &
                source%block_parities(block_index))
        end do
    end subroutine install_descriptor_header

    ! Look up the subshell state using occupation, seniority and the angular momentum available in the
    ! descriptor.
    function resolve_subshell_state(single_two_j, occupation, physical_seniority, printed_two_j, block_index, row_index, &
        & status) &
        result(subshell_two_j)
        ! Arguments
        integer, intent(in) :: single_two_j, occupation, physical_seniority, printed_two_j
        integer, intent(in) :: block_index, row_index
        integer, intent(out), optional :: status
        integer :: subshell_two_j

        ! Local declarations
        integer :: subshell_capacity, population, table_row, table_entry, table_offset, match_count, matching_j_count
        integer :: candidate_two_j, candidate_seniority, candidate_label
        subshell_two_j = -1
        subshell_capacity = single_two_j + 1
        population = min(occupation, subshell_capacity - occupation)
        if (population == 0) then
            if (.not. require_state(physical_seniority == 0 .and. (printed_two_j == -1 .or. printed_two_j == 0), &
                "empty/closed shell must have seniority zero and J=0 or omitted J")) return
            subshell_two_j = 0
            return
        end if
        if (.not. require_state(single_two_j <= 17, "unsupported open-shell state table for single-electron 2j > 17")) return
        if (.not. require_state(single_two_j <= 9 .or. population <= 2, "unsupported high-j occupation")) return
        table_row = (subshell_capacity-2)/2
        table_row = table_row*(table_row+1)/2 + population + 1
        if (.not. require_state(table_row >= 1 .and. table_row <= size(state_counts), "unsupported subshell table row")) return
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
            if (.not. require_state(candidate_label == 0, "unsupported state requires an additional w label")) return
            subshell_two_j = candidate_two_j
        end do
        if (.not. require_state(match_count == 1, "invalid or ambiguous subshell (occupation, seniority, J)")) return
        matching_j_count = 0
        do table_entry = 1, state_counts(table_row)
            table_offset = state_offsets(table_row) + 3*(table_entry-1)
            if (state_values(table_offset+3) == subshell_two_j+1) matching_j_count = matching_j_count + 1
        end do
        ! For j=7/2 with four electrons, seniority distinguishes states with the same J.
        if (.not. require_state(matching_j_count == 1 .or. (single_two_j == 7 .and. occupation == 4), &
            "unsupported state is ambiguous in the ordinary text reader")) return
    contains
        logical function require_state(condition, message) result(ok)
            ! Arguments
            logical, intent(in) :: condition
            character(len=*), intent(in) :: message
            ok = condition
            if (present(status)) status = 0
            if (.not. condition) then
                write(error_unit, "(A,I0,A,I0,2A)") "CSF descriptor input: block ", block_index, &
                    ", CSF ", row_index, ": ", message
                if (present(status)) then
                    status = 1
                else
                    error stop 1
                end if
            end if
        end function require_state
    end function resolve_subshell_state

    ! Insert full core shells and validate the active-shell row before packing the coupling codes.
    subroutine decode_descriptor_row(source, block_index, row_index, packed_row, occupations, subshell_quantum_numbers, &
        & couplings, status)
        ! Arguments
        type(descriptor_source), intent(in) :: source
        integer, intent(in) :: block_index, row_index
        integer(int8), intent(in) :: packed_row(:)
        integer(int8), intent(out) :: occupations(:), subshell_quantum_numbers(:,:), couplings(:)
        integer, intent(out), optional :: status

        ! Local declarations
        integer :: shell_two_j(orbital_capacity), printed_two_k(orbital_capacity)
        integer :: peel_index, orbital_index, channel_offset, occupation, seniority, printed_two_j
        integer :: occupied_count, occupied_position, open_count, subshell_capacity, previous_two_j, coupled_two_j
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
            subshell_capacity = source%single_electron_two_j(orbital_index) + 1
            occupations(orbital_index) = int(subshell_capacity, int8)
            subshell_quantum_numbers(orbital_index,3) = 1_int8
            electron_count = electron_count + subshell_capacity
        end do
        ! V2 stores (occupation, printed 2J, physical seniority, printed 2K) for each peel shell.
        do peel_index = 1, source%peel_count
            orbital_index = source%core_count + peel_index
            channel_offset = 4*(peel_index-1)
            occupation = int(packed_row(channel_offset+1))
            printed_two_j = int(packed_row(channel_offset+2))
            seniority = int(packed_row(channel_offset+3))
            printed_two_k(orbital_index) = int(packed_row(channel_offset+4))
            subshell_capacity = source%single_electron_two_j(orbital_index) + 1
            if (.not. require_row(occupation >= 0 .and. occupation <= subshell_capacity, &
                & "occupation outside shell capacity")) return
            if (.not. require_row(seniority >= 0, "negative physical seniority")) return
            if (.not. require_row(printed_two_j >= -1 .and. printed_two_j <= 126, "invalid subshell 2J")) return
            if (.not. require_row(printed_two_k(orbital_index) >= -1 .and. printed_two_k(orbital_index) <= 126, &
                "invalid intermediate 2K")) return
            if (occupation == 0) then
                if (.not. require_row(printed_two_j == -1 .and. seniority == 0 .and. &
                    printed_two_k(orbital_index) == -1, "noncanonical empty shell descriptor")) return
            end if
            shell_two_j(orbital_index) = resolve_subshell_state(source%single_electron_two_j(orbital_index), &
                occupation, seniority, printed_two_j, block_index, row_index, &
                & status)
            if (failed(status)) return
            occupations(orbital_index) = int(occupation, int8)
            subshell_quantum_numbers(orbital_index,1) = int(seniority, int8)
            subshell_quantum_numbers(orbital_index,3) = int(shell_two_j(orbital_index)+1, int8)
            if (shell_two_j(orbital_index) /= 0) last_nonzero_shell = orbital_index
            electron_count = electron_count + occupation
            parity_exponent = parity_exponent + occupation*source%orbital_momenta(orbital_index)
        end do
        total_two_j = int(packed_row(source%row_width-1))
        parity = int(packed_row(source%row_width))
        if (.not. require_row(total_two_j == source%block_two_j(block_index) .and. &
            parity == source%block_parities(block_index), "row/block symmetry mismatch")) return
        if (.not. require_row(parity == (-1)**parity_exponent, "incorrect row parity")) return
        if (.not. require_row(electron_count == source%electron_count, "incorrect row electron count")) return
        ! Infer only trivial or final couplings. Intermediate nontrivial values must be supplied by the
        ! descriptor.
        occupied_count = count(occupations(source%core_count+1:source%orbital_count) > 0)
        occupied_position = 0
        open_count = 0
        previous_two_j = 0
        do orbital_index = source%core_count+1, source%orbital_count
            occupation = int(occupations(orbital_index))
            if (occupation == 0) cycle
            occupied_position = occupied_position + 1
            subshell_capacity = source%single_electron_two_j(orbital_index) + 1
            coupled_two_j = printed_two_k(orbital_index)
            if (coupled_two_j >= 0) then
                if (.not. require_row(occupied_position >= 2 .and. occupied_position < occupied_count, &
                    "intermediate coupling at invalid occupied-shell boundary")) return
            else if (orbital_index >= last_nonzero_shell) then
                coupled_two_j = total_two_j
            else if (shell_two_j(orbital_index) == 0) then
                coupled_two_j = previous_two_j
            else if (previous_two_j == 0) then
                coupled_two_j = shell_two_j(orbital_index)
            else
                if (.not. require_row(.false., "missing nontrivial intermediate coupling")) return
            end if
            if (.not. require_row(coupled_two_j >= abs(previous_two_j-shell_two_j(orbital_index)) .and. &
                coupled_two_j <= previous_two_j+shell_two_j(orbital_index) .and. &
                mod(previous_two_j+shell_two_j(orbital_index)-coupled_two_j, 2) == 0, &
                "coupling violates angular momentum triangle/parity")) return
            previous_two_j = coupled_two_j
            if (occupation < subshell_capacity) then
                open_count = open_count + 1
                if (open_count > 1) couplings(open_count-1) = int(coupled_two_j+1, int8)
            end if
        end do
        if (.not. require_row(previous_two_j == total_two_j, "incorrect final angular momentum")) return
        ! Store the total 2J+1 and parity together in the final coupling slot.
        couplings(orbital_capacity) = int((total_two_j+1)*parity, int8)
    contains
        logical function require_row(condition, message) result(ok)
            ! Arguments
            logical, intent(in) :: condition
            character(len=*), intent(in) :: message
            ok = condition
            if (present(status)) status = 0
            if (.not. condition) then
                write(error_unit, "(A,I0,A,I0,2A)") "CSF descriptor input: block ", block_index, &
                    ", CSF ", row_index, ": ", message
                if (present(status)) then
                    status = 1
                else
                    error stop 1
                end if
            end if
        end function require_row
    end subroutine decode_descriptor_row

    subroutine read_descriptor_block(source, block_index, occupations, subshell_quantum_numbers, couplings, status)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        integer, intent(in) :: block_index
        integer(int8), intent(out) :: occupations(:,:)
        integer(int8), intent(out), optional :: subshell_quantum_numbers(:,:,:), couplings(:,:)
        integer, intent(out), optional :: status

        ! Local declarations
        integer :: block_status
        call read_descriptor_block_impl(source,block_index,occupations,subshell_quantum_numbers,couplings,block_status)
        if (block_status /= 0) then
            occupations=0_int8
            if (present(subshell_quantum_numbers)) subshell_quantum_numbers=0_int8
            if (present(couplings)) couplings=0_int8
            call close_descriptor_source(source)
        end if
        if (present(status)) then
            status=block_status
        else if (block_status /= 0) then
            error stop "Descriptor block validation failed"
        end if
    end subroutine read_descriptor_block

    subroutine read_descriptor_block_impl(source, block_index, occupations, subshell_quantum_numbers, couplings, status)
        ! Arguments and local declarations
        type(descriptor_source), intent(inout) :: source
        integer, intent(in) :: block_index
        integer(int8), intent(out) :: occupations(:,:)
        integer(int8), intent(out), optional :: subshell_quantum_numbers(:,:,:), couplings(:,:)
        integer(int8) :: row_subshell_quantum_numbers(orbital_capacity,3), row_couplings(orbital_capacity)
        integer, intent(out), optional :: status
        integer(int8), allocatable :: batch_values(:,:)
        character(kind=c_char,len=1024) :: reader_message
        integer :: batch_capacity, batch_row_count, row_offset, batch_row_index, read_status, allocation_status
        reader_message = c_null_char
        if (.not. require_descriptor(source%is_open .and. source%header_loaded, "descriptor source is not open", status)) return
        if (.not. require_descriptor(block_index == source%next_block_index .and. block_index <= source%block_count, &
            "blocks must be loaded once in original order", status)) return
        if (.not. require_descriptor(size(occupations,1) == orbital_capacity .and. &
            size(occupations,2) == source%block_lengths(block_index),  &
            & "occupation array shape mismatch", status)) return
        if (present(couplings)) then
            if (.not. require_descriptor(all(shape(couplings) == shape(occupations)), "coupling shape", status)) return
        end if
        if (present(subshell_quantum_numbers)) then
            if (.not. require_descriptor(size(subshell_quantum_numbers,1) == orbital_capacity .and. &
                size(subshell_quantum_numbers,2) == 3 .and. &
                size(subshell_quantum_numbers,3) == size(occupations,2), "quantum shape", status)) return
        end if
        ! Read bounded batches from the current block, in the same order as the file.
        batch_capacity = min(8192, source%block_lengths(block_index))
        ! Fortran column j contains all channels of descriptor row j, contiguous across the ABI.
        allocate(batch_values(source%row_width,batch_capacity), stat=allocation_status)
        if (.not. require_descriptor(allocation_status == 0, "cannot allocate descriptor batch", status)) return
        row_offset = 0
        do while (row_offset < source%block_lengths(block_index))
            batch_row_count = min(batch_capacity, source%block_lengths(block_index)-row_offset)
            read_status = rust_descriptor_read(source%handle, batch_values, &
                int(source%row_width,c_size_t)*int(batch_row_count,c_size_t), int(batch_row_count,c_size_t), &
                reader_message, 1024_c_size_t)
            if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), &
                & status)) return
            do batch_row_index = 1, batch_row_count
                call decode_descriptor_row(source, block_index, row_offset+batch_row_index, batch_values(:,batch_row_index), &
                    occupations(:,row_offset+batch_row_index), row_subshell_quantum_numbers, &
                    & row_couplings, status)
                if (failed(status)) return
                if (present(subshell_quantum_numbers)) subshell_quantum_numbers(:,:, &
                    & row_offset+batch_row_index) = row_subshell_quantum_numbers
                if (present(couplings)) couplings(:,row_offset+batch_row_index) = row_couplings
            end do
            row_offset = row_offset + batch_row_count
        end do
        source%next_block_index = source%next_block_index + 1
    end subroutine read_descriptor_block_impl
end module csf_descriptor_reader
