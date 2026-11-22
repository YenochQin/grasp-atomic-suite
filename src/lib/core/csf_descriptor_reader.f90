! Read and validate CSF descriptor headers and blocks, then decode the rows into Fortran basis arrays.
module csf_descriptor_reader
    ! Module dependencies
    use iso_fortran_env, only: int8, int64, error_unit
    use iso_c_binding, only: c_ptr, c_null_ptr, c_null_char, c_char, c_size_t, c_int32_t, c_int64_t
    use csf_descriptor_c_bindings
    use grasp_parameters, only: max_orbitals, csf_orbital_letters
    use angular_symmetry_labels, only: format_j_label
    use csf_coupling_data, only: csf_angular_code_count, csf_total_j_slot, csf_parity_slot
    implicit none (type, external)

    ! Module interface
    private
    public :: descriptor_source, open_descriptor_source, close_descriptor_source
    public :: read_descriptor_block, install_descriptor_header, descriptor_block_label
    public :: read_descriptor_zero_order_counts

    ! Derived types
    ! Owns an indexed Rust descriptor reader. Close the source explicitly and avoid copying a live handle.
    type :: descriptor_source
        ! Components
        type(c_ptr) :: handle = c_null_ptr
        logical :: is_open = .false., header_loaded = .false.
        integer :: core_count = 0, peel_count = 0, orbital_count = 0, block_count = 0
        integer :: format_version = 0
        character(len=64) :: basis_id = ""
        integer :: row_width = 0, total_csf_count = 0, electron_count = 0
        integer, allocatable :: principal_numbers(:), relativistic_kappas(:)
        integer, allocatable :: orbital_momenta(:), single_electron_two_j(:)
        character(len=2), allocatable :: orbital_labels(:)
        integer, allocatable :: block_lengths(:), block_two_j(:), block_parities(:)
    end type descriptor_source

contains

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

    subroutine open_descriptor_source(source, filename, maximum_blocks, header_filename, status, spin_angular)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        character(len=*), intent(in) :: filename
        integer, intent(in) :: maximum_blocks
        character(len=*), intent(in), optional :: header_filename
        integer, intent(out), optional :: status
        logical, intent(in), optional :: spin_angular

        ! Local declarations
        integer :: open_status
        character(kind=c_char,len=1024) :: message
        call open_descriptor_source_impl(source, filename, maximum_blocks, header_filename, open_status)
        if (open_status == 0 .and. present(spin_angular)) then
            if (spin_angular) then
                open_status = rust_descriptor_check_spin_angular(source%handle, message, len(message,kind=c_size_t))
                if (open_status /= 0) write(error_unit, "(2A)") "CSF descriptor input: ", message(:index(message,c_null_char)-1)
            end if
        end if
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
        reader_message = c_null_char
        call close_descriptor_source(source)
        source = descriptor_source()
        if (.not. require_descriptor(rust_descriptor_abi_version() == 3, "unsupported Rust ABI", status)) return
        if (.not. require_descriptor(rust_descriptor_format_version() == 3, "unsupported descriptor format", status)) return
        path_length = len_trim(filename)
        if (.not. require_descriptor(path_length > 8, "expected descriptor Parquet path", status)) return
        if (.not. require_descriptor(filename(path_length-7:path_length) == ".parquet", &
            & "only descriptor Parquet input is supported", status)) return
        descriptor_stem = filename(:path_length-8)
        path_length = len(descriptor_stem)
        if (path_length >= len("_desc")) then
            if (descriptor_stem(path_length-len("_desc")+1:) == "_desc") &
                descriptor_stem = descriptor_stem(:path_length-len("_desc"))
        end if
        if (path_length >= len("_descriptors") .and. len(descriptor_stem) == path_length) then
            if (descriptor_stem(path_length-len("_descriptors")+1:) == "_descriptors") &
                descriptor_stem = descriptor_stem(:path_length-len("_descriptors"))
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
        call rust_descriptor_basis_id(source%handle, source%basis_id)
        source%format_version = int(rust_descriptor_format_version())
        read_status = rust_descriptor_info(source%handle, header_counts, reader_message, 1024_c_size_t)
        if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), &
            & status)) return
        if (.not. require_descriptor(all(header_counts >= 0) .and. all(header_counts <= int(huge(0),int64)), &
            & "header integer overflow",  &
            & status)) return
        source%core_count = int(header_counts(1))
        source%peel_count = int(header_counts(2))
        if (.not. require_descriptor(header_counts(1)+header_counts(2) <= max_orbitals, "orbital capacity exceeded", &
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
            if (.not. require_descriptor(kappa /= 0 .and. orbital_l_value >= 0 .and. orbital_l_value < len(csf_orbital_letters), &
                & "unsupported orbital", status)) return
            source%orbital_momenta(orbital_index) = orbital_l_value
            source%single_electron_two_j(orbital_index) = 2*abs(kappa)-1
            source%orbital_labels(orbital_index) = csf_orbital_letters(orbital_l_value+1:orbital_l_value+1)//" "
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
        if (.not. require_descriptor(all(raw_zero_order_counts > 0), &
            "using zero-first requires positive zero_order_counts in every block", status)) return
        counts = int(raw_zero_order_counts)
    end subroutine read_descriptor_zero_order_counts

    subroutine close_descriptor_source(source)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        call rust_descriptor_close(source%handle)
        source%is_open = .false.
        source%header_loaded = .false.
    end subroutine close_descriptor_source

    ! Keep the legacy four-column J field, expanding it for three-digit half-integers.
    function descriptor_block_label(total_two_j, parity) result(label)
        ! Arguments
        integer, intent(in) :: total_two_j, parity
        character(len=8) :: label

        label = format_j_label(total_two_j+1)//merge("+","-",parity > 0)
    end function descriptor_block_label

    subroutine install_descriptor_header(source, core_count, block_labels, status)
        ! Module dependencies
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
            orbital_symmetry_labels
        use atomic_calculation_state, only: electron_count
        use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block

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
        if (associated(csf_counts_by_block)) deallocate(csf_counts_by_block)
        csf_counts_by_block=>prepared_lengths
        core_count = source%core_count
        orbital_count = source%orbital_count
        csf_count = source%total_csf_count
        electron_count = source%electron_count
        symmetry_block_count = source%block_count
        orbital_principal_n(:orbital_count) = source%principal_numbers
        orbital_kappa(:orbital_count) = source%relativistic_kappas
        orbital_l(:orbital_count) = source%orbital_momenta
        orbital_two_j(:orbital_count) = source%single_electron_two_j
        orbital_symmetry_labels(:orbital_count) = source%orbital_labels
        csf_counts_by_block = 0
        csf_counts_by_block(:symmetry_block_count) = source%block_lengths
        block_labels = ""
        do block_index = 1, symmetry_block_count
            block_labels(block_index) = descriptor_block_label(source%block_two_j(block_index), &
                source%block_parities(block_index))
        end do
    end subroutine install_descriptor_header

    ! Rust has validated each V3 row. Install the core and copy the complete peel codes.
    subroutine install_descriptor_rows(source, rows, occupations, subshell_quantum_numbers, angular_codes)
        type(descriptor_source), intent(in) :: source
        integer(int8), intent(in) :: rows(:,:)
        integer(int8), intent(out) :: occupations(:,:)
        integer(int8), intent(out), optional :: subshell_quantum_numbers(:,:,:), angular_codes(:,:)
        integer :: csf_index, core_end, orbital_end

        core_end = source%core_count
        orbital_end = source%orbital_count
        do csf_index = 1, size(rows,2)
            occupations(:core_end,csf_index) = int(source%single_electron_two_j(:core_end)+1,int8)
            occupations(core_end+1:orbital_end,csf_index) = rows(1:source%row_width-2:4,csf_index)
            occupations(orbital_end+1:,csf_index) = 0_int8
            if (present(subshell_quantum_numbers)) then
                subshell_quantum_numbers(:,:,csf_index) = 0_int8
                subshell_quantum_numbers(:core_end,3,csf_index) = 1_int8
                subshell_quantum_numbers(core_end+1:orbital_end,1,csf_index) = rows(2:source%row_width-2:4,csf_index)
                subshell_quantum_numbers(core_end+1:orbital_end,3,csf_index) = rows(3:source%row_width-2:4,csf_index)
            end if
            if (present(angular_codes)) then
                angular_codes(:core_end,csf_index) = 1_int8
                angular_codes(core_end+1:orbital_end,csf_index) = rows(4:source%row_width-2:4,csf_index)
                angular_codes(orbital_end+1:max_orbitals,csf_index) = 0_int8
                angular_codes(csf_total_j_slot,csf_index) = rows(source%row_width-1,csf_index)
                angular_codes(csf_parity_slot,csf_index) = rows(source%row_width,csf_index)
            end if
        end do
    end subroutine install_descriptor_rows

    subroutine read_descriptor_block(source, block_index, occupations, subshell_quantum_numbers, angular_codes, status)
        ! Arguments
        type(descriptor_source), intent(inout) :: source
        integer, intent(in) :: block_index
        integer(int8), intent(out) :: occupations(:,:)
        integer(int8), intent(out), optional :: subshell_quantum_numbers(:,:,:), angular_codes(:,:)
        integer, intent(out), optional :: status

        ! Local declarations
        integer :: block_status
        call read_descriptor_block_impl(source,block_index,occupations,subshell_quantum_numbers,angular_codes,block_status)
        if (block_status /= 0) then
            occupations=0_int8
            if (present(subshell_quantum_numbers)) subshell_quantum_numbers=0_int8
            if (present(angular_codes)) angular_codes=0_int8
            call close_descriptor_source(source)
        end if
        if (present(status)) then
            status=block_status
        else if (block_status /= 0) then
            error stop "Descriptor block validation failed"
        end if
    end subroutine read_descriptor_block

    subroutine read_descriptor_block_impl(source, block_index, occupations, subshell_quantum_numbers, angular_codes, status)
        ! Arguments and local declarations
        type(descriptor_source), intent(inout) :: source
        integer, intent(in) :: block_index
        integer(int8), intent(out) :: occupations(:,:)
        integer(int8), intent(out), optional, target :: subshell_quantum_numbers(:,:,:), angular_codes(:,:)
        integer, intent(out), optional :: status
        integer(int8), allocatable :: batch_values(:,:)
        integer(int8), pointer :: batch_quantum(:,:,:), batch_angular(:,:)
        character(kind=c_char,len=1024) :: reader_message
        integer :: batch_capacity, batch_row_count, row_offset, last_row, read_status, allocation_status
        reader_message = c_null_char
        if (.not. require_descriptor(source%is_open .and. source%header_loaded, "descriptor source is not open", status)) return
        if (.not. require_descriptor(block_index >= 1 .and. block_index <= source%block_count, &
            "block index out of range", status)) return
        if (.not. require_descriptor(size(occupations,1) == max_orbitals .and. &
            size(occupations,2) == source%block_lengths(block_index),  &
            & "occupation array shape mismatch", status)) return
        if (present(angular_codes)) then
            if (.not. require_descriptor(all(shape(angular_codes) == [csf_angular_code_count,size(occupations,2)]), &
                "angular code shape", status)) return
        end if
        if (present(subshell_quantum_numbers)) then
            if (.not. require_descriptor(size(subshell_quantum_numbers,1) == max_orbitals .and. &
                size(subshell_quantum_numbers,2) == 3 .and. &
                size(subshell_quantum_numbers,3) == size(occupations,2), "quantum shape", status)) return
        end if
        read_status = rust_descriptor_select_block(source%handle, int(block_index-1,c_size_t), reader_message, 1024_c_size_t)
        if (.not. require_descriptor(read_status == 0, reader_message(:max(0,index(reader_message,c_null_char)-1)), status)) return
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
            last_row = row_offset + batch_row_count
            nullify(batch_quantum,batch_angular)
            if (present(subshell_quantum_numbers)) batch_quantum => subshell_quantum_numbers(:,:,row_offset+1:last_row)
            if (present(angular_codes)) batch_angular => angular_codes(:,row_offset+1:last_row)
            call install_descriptor_rows(source,batch_values(:,:batch_row_count),occupations(:,row_offset+1:last_row), &
                batch_quantum,batch_angular)
            row_offset = last_row
        end do
    end subroutine read_descriptor_block_impl
end module csf_descriptor_reader
