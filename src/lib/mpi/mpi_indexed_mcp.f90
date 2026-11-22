! Rank-local MCP stream format. All file positions are one-based byte positions.
module mpi_indexed_mcp
    use iso_c_binding, only: c_char, c_int, c_null_char
    use iso_fortran_env, only: int32, int64, real64, file_storage_size
    use ieee_arithmetic, only: ieee_support_datatype
    use grasp_parameters, only: orbital_index_radix
    use mpi_state, only: mpi_rank, mpi_process_count, mpi_bcast, mpi_comm_world, mpi_integer8, mpi_error_code, &
        mpi_barrier, mpi_character
    use mpi_io_errors, only: fail_input_io_mpi
    implicit none (type, external)
    private

    public :: indexed_mcp, open_indexed_mcp_writer, finish_indexed_mcp, open_indexed_mcp_reader, close_indexed_mcp
    public :: write_indexed_sparse, write_indexed_group, start_indexed_channel, read_indexed_sparse
    public :: read_indexed_group, read_indexed_contributions, indexed_column_count, decode_mcp_label
    public :: mcp_segment_size, mcp_group, contribution_bytes, mcp_group_directory_bytes
    public :: read_indexed_directory, validate_indexed_group, check_sparse_structure

    character(len=8), parameter :: magic = "GMCPIDX1", complete_magic = "MCPDONE1"
    integer(int32), parameter :: layout(8) = [1_int32, int(z'01020304',int32), 4_int32, 8_int32, 8_int32, &
        int(orbital_index_radix,int32), 1_int32, 1_int32]
    integer(int64), parameter :: header_bytes = 8 + 4*size(layout) + 8*12 + 64
    integer(int64), parameter :: block_entry_bytes = 8 + 8*5, channel_entry_bytes = 8*3
    integer(int64), parameter :: footer_bytes = 8 + 8*3, mcp_group_directory_bytes = 8*4
    integer(int64), parameter :: contribution_bytes = 16
    integer, parameter :: mcp_segment_size = 65536

    type :: mcp_group
        integer(int64) :: position, bytes, label, count
    end type
    type :: indexed_mcp
        integer :: unit = -1, scratch = -1
        integer :: core_count, block_count, max_rank, electron_count, csf_count, orbital_count
        integer :: interaction_flags
        character(len=64) :: basis_id
        character(:), allocatable :: output_path, staging_path
        character(len=8), allocatable :: block_labels(:)
        integer(int64), allocatable :: block_sizes(:), cutoffs(:), sparse_positions(:), row_counts(:), column_counts(:)
        ! Channel 1 is T; channels 2 and above are V(0), V(1), ...
        integer(int64), allocatable :: group_counts(:,:), first_groups(:,:), channel_counts(:,:)
        integer(int64) :: group_count = 0, directory, generation, next_position
        integer :: block = 0, channel = 0
    end type

    interface
        function rename_shard(source, destination) bind(C, name="rename") result(status)
            import c_char, c_int
            character(kind=c_char), intent(in) :: source(*), destination(*)
            integer(c_int) :: status
        end function rename_shard
    end interface

contains

    subroutine check_format(unit)
        integer, intent(in) :: unit
        if (file_storage_size /= 8 .or. storage_size(0) /= 32 .or. .not. ieee_support_datatype(0.0_real64)) &
            call fail_input_io_mpi(0, unit, "MCP indexed unsupported integer, real, or file storage representation")
    end subroutine check_format

    subroutine verify(condition, file, context, status)
        logical, intent(in) :: condition
        type(indexed_mcp), intent(in) :: file
        character(*), intent(in) :: context
        integer, intent(in), optional :: status
        if (condition) return
        if (present(status)) call fail_input_io_mpi(status, file%unit, "MCP indexed "//context)
        call fail_input_io_mpi(0, file%unit, "MCP indexed "//context)
    end subroutine verify

    pure integer(int64) function table_end(file) result(position)
        type(indexed_mcp), intent(in) :: file
        position = header_bytes + file%block_count*(block_entry_bytes + (file%max_rank+2)*channel_entry_bytes) + 1
    end function table_end

    subroutine allocate_block_tables(file)
        type(indexed_mcp), intent(inout) :: file
        integer :: blocks, channels, status
        blocks = file%block_count
        channels = file%max_rank + 2
        allocate(file%block_labels(blocks), file%block_sizes(blocks), file%cutoffs(blocks), file%sparse_positions(blocks), &
            file%row_counts(blocks), file%column_counts(blocks), file%group_counts(channels,blocks), &
            file%first_groups(channels,blocks), file%channel_counts(channels,blocks), stat=status)
        call verify(status == 0, file, "allocate block directory", status)
    end subroutine allocate_block_tables

    subroutine write_header(file)
        type(indexed_mcp), intent(in) :: file
        integer :: status, block
        write(file%unit, pos=1, iostat=status) magic, layout, int(mpi_rank,int64), int(mpi_process_count,int64), &
            int(file%core_count,int64), int(file%block_count,int64), int(file%max_rank,int64), int(file%electron_count,int64), &
            int(file%csf_count,int64), int(file%orbital_count,int64), int(file%interaction_flags,int64), &
            file%generation, file%directory, file%group_count, file%basis_id
        call verify(status == 0, file, "write header", status)
        do block = 1, file%block_count
            write(file%unit, iostat=status) file%block_labels(block), file%block_sizes(block), file%cutoffs(block), &
                file%sparse_positions(block), file%column_counts(block), file%row_counts(block), &
                file%first_groups(:,block), file%group_counts(:,block), file%channel_counts(:,block)
            call verify(status == 0, file, "write block directory", status)
        end do
    end subroutine write_header

    subroutine open_indexed_mcp_writer(file, prefix, core, max_rank, electrons, orbitals, sizes, labels, cutoffs, flags, basis_id)
        type(indexed_mcp), intent(out) :: file
        character(*), intent(in) :: prefix, labels(:)
        character(len=64), intent(in) :: basis_id
        integer, intent(in) :: core, max_rank, electrons, orbitals, sizes(:), cutoffs(:), flags
        integer :: status
        integer(int64) :: total_csf
        character(len=32) :: suffix
        call check_format(-1)
        total_csf = sum(int(sizes,int64))
        call verify(size(sizes) > 0 .and. size(labels) == size(sizes) .and. size(cutoffs) == size(sizes) .and. &
            all(sizes > 0) .and. total_csf < huge(0), file, "invalid writer basis dimensions")
        if (mpi_rank == 0) call system_clock(file%generation)
        call mpi_bcast(file%generation, 1, mpi_integer8, 0, mpi_comm_world, mpi_error_code)
        write(suffix,"(I0)") file%generation
        file%output_path = trim(prefix)//".mcpidx"
        file%staging_path = file%output_path//".part-"//trim(suffix)
        open(newunit=file%unit, file=file%staging_path, form="unformatted", access="stream", &
            status="new", action="readwrite", iostat=status)
        call verify(status == 0, file, "create "//file%staging_path, status)
        open(newunit=file%scratch, form="unformatted", access="stream", status="scratch", iostat=status)
        call verify(status == 0, file, "open temporary group directory", status)
        file%core_count = core
        file%max_rank = max_rank
        file%electron_count = electrons
        file%orbital_count = orbitals
        file%block_count = size(sizes)
        file%csf_count = int(total_csf)
        file%interaction_flags = flags
        file%basis_id = basis_id
        call allocate_block_tables(file)
        file%block_labels = labels
        file%block_sizes = sizes
        file%cutoffs = cutoffs
        file%sparse_positions = 0
        file%column_counts = 0
        file%row_counts = 0
        file%first_groups = 0
        file%group_counts = 0
        file%channel_counts = 0
        file%directory = 0
        call write_header(file)
        file%next_position = table_end(file)
    end subroutine open_indexed_mcp_writer

    pure integer function indexed_column_count(csf_count, rank, ranks) result(columns)
        integer, intent(in) :: csf_count, rank, ranks
        columns = 0
        if (csf_count > rank) columns = 1 + (csf_count-1-rank)/ranks
    end function indexed_column_count

    subroutine write_indexed_sparse(file, block, ends, rows)
        type(indexed_mcp), intent(inout) :: file
        integer, intent(in) :: block, ends(:), rows(:)
        integer :: status
        call verify(block == file%block+1 .and. block <= file%block_count, file, "out-of-order sparse block")
        if (file%block > 0) call verify(file%channel == file%max_rank+2, file, "incomplete previous block")
        call verify(size(ends) == indexed_column_count(int(file%block_sizes(block)),mpi_rank,mpi_process_count), &
            file, "invalid compact column count")
        call verify(size(rows,kind=int64) < huge(0), file, "sparse positions exceed 32-bit capacity")
        call check_sparse_structure(ends,rows,file%unit,"MCP indexed write")
        file%block = block
        file%channel = 0
        file%sparse_positions(block) = file%next_position
        file%column_counts(block) = size(ends,kind=int64)
        file%row_counts(block) = size(rows,kind=int64)
        write(file%unit, pos=file%next_position, iostat=status) ends, rows
        call verify(status == 0, file, "write sparse structure", status)
        file%next_position = file%next_position + 4*(file%column_counts(block)+file%row_counts(block))
    end subroutine write_indexed_sparse

    subroutine start_indexed_channel(file, channel)
        type(indexed_mcp), intent(inout) :: file
        integer, intent(in) :: channel
        call verify(file%block > 0 .and. channel == file%channel+1 .and. channel <= file%max_rank+2, &
            file, "out-of-order channel")
        file%channel = channel
        file%first_groups(channel,file%block) = file%group_count + 1
    end subroutine start_indexed_channel

    subroutine write_indexed_group(file, label, columns, positions, coefficients)
        type(indexed_mcp), intent(inout) :: file
        integer, intent(in) :: label, columns(:), positions(:)
        real(real64), intent(in) :: coefficients(:)
        integer :: status, channel, block, i
        integer(int64) :: count, bytes, first, last
        call verify(file%block > 0 .and. file%channel > 0, file, "group written before channel")
        count = size(coefficients,kind=int64)
        call verify(count > 0 .and. count < huge(0) .and. size(columns,kind=int64) == count .and. &
            size(positions,kind=int64) == count, file, "invalid group size for in-memory sorter")
        bytes = contribution_bytes*count
        call verify(bytes <= huge(bytes)-file%next_position, file, "contribution byte offset overflow")
        write(file%scratch, iostat=status) file%next_position, bytes, int(label,int64), count
        call verify(status == 0, file, "write group directory", status)
        do first = 1, count, mcp_segment_size
            last = min(count,first+mcp_segment_size-1)
            write(file%unit, pos=file%next_position, iostat=status) (columns(i),positions(i),coefficients(i),i=int(first),int(last))
            call verify(status == 0, file, "write contribution segment", status)
            file%next_position = file%next_position + contribution_bytes*(last-first+1)
        end do
        block = file%block
        channel = file%channel
        file%group_count = file%group_count + 1
        file%group_counts(channel,block) = file%group_counts(channel,block) + 1
        file%channel_counts(channel,block) = file%channel_counts(channel,block) + count
    end subroutine write_indexed_group

    subroutine finish_indexed_mcp(file)
        type(indexed_mcp), intent(inout) :: file
        integer(int64) :: first
        integer :: status, entries
        integer(int64), allocatable :: buffer(:,:)
        call verify(file%block == file%block_count .and. file%channel == file%max_rank+2, file, "incomplete output")
        allocate(buffer(4,mcp_segment_size))
        file%directory = file%next_position
        call verify(file%group_count <= (huge(first)-file%directory-footer_bytes)/mcp_group_directory_bytes, &
            file, "directory byte offset overflow")
        do first = 1, file%group_count, mcp_segment_size
            entries = int(min(int(mcp_segment_size,int64),file%group_count-first+1))
            read(file%scratch, pos=1+mcp_group_directory_bytes*(first-1), iostat=status) buffer(:,:entries)
            call verify(status == 0, file, "read temporary directory", status)
            write(file%unit, pos=file%next_position, iostat=status) buffer(:,:entries)
            call verify(status == 0, file, "append group directory", status)
            file%next_position = file%next_position + mcp_group_directory_bytes*entries
        end do
        close(file%scratch, iostat=status)
        call verify(status == 0, file, "close temporary directory", status)
        call write_header(file)
        flush(file%unit, iostat=status)
        call verify(status == 0, file, "flush shard before completion", status)
        call mpi_barrier(mpi_comm_world, mpi_error_code)
        write(file%unit, pos=file%next_position, iostat=status) complete_magic, file%generation, file%directory, file%group_count
        call verify(status == 0, file, "write completion footer", status)
        call close_indexed_mcp(file)
        call mpi_barrier(mpi_comm_world, mpi_error_code)
        status = int(rename_shard(file%staging_path//c_null_char,file%output_path//c_null_char))
        call verify(status == 0, file, "publish "//file%output_path, status)
        call mpi_barrier(mpi_comm_world, mpi_error_code)
    end subroutine finish_indexed_mcp

    subroutine close_indexed_mcp(file)
        type(indexed_mcp), intent(inout) :: file
        integer :: status
        close(file%unit, iostat=status)
        call verify(status == 0, file, "close shard", status)
        file%unit = -1
    end subroutine close_indexed_mcp

    subroutine open_indexed_mcp_reader(file, prefix, maximum_blocks)
        type(indexed_mcp), intent(out) :: file
        character(*), intent(in) :: prefix
        integer, intent(in) :: maximum_blocks
        character(len=8) :: tag
        integer(int32) :: file_layout(size(layout))
        integer(int64) :: fields(9), file_size, footer(3), identity(8), expected_group, expected_position, remaining
        integer(int64), allocatable :: root_partitions(:,:)
        character(len=64) :: root_basis_id
        character(len=8), allocatable :: root_labels(:)
        integer :: status, block, channel
        call check_format(-1)
        open(newunit=file%unit, file=trim(prefix)//".mcpidx", form="unformatted", access="stream", &
            status="old", action="read", iostat=status)
        call verify(status == 0, file, "open "//trim(prefix)//".mcpidx", status)
        read(file%unit, iostat=status) tag, file_layout, fields, file%generation, file%directory, file%group_count, file%basis_id
        call verify(status == 0, file, "read header", status)
        call verify(tag == magic .and. all(file_layout == layout), file, "incompatible format, endianness, or label radix")
        call verify(fields(1) == mpi_rank .and. fields(2) == mpi_process_count, file, "incompatible rank layout")
        call verify(all(fields(3:) >= 0) .and. all(fields(3:) < huge(0)), file, "invalid header dimensions")
        call verify(fields(4) >= 1 .and. fields(4) <= maximum_blocks .and. fields(5) < orbital_index_radix .and. &
            fields(8) > 0 .and. fields(8) < orbital_index_radix .and. fields(3) <= fields(8) .and. fields(9) <= 3, &
            file, "invalid basis or interaction metadata")
        file%core_count = int(fields(3))
        file%block_count = int(fields(4))
        file%max_rank = int(fields(5))
        file%electron_count = int(fields(6))
        file%csf_count = int(fields(7))
        file%orbital_count = int(fields(8))
        file%interaction_flags = int(fields(9))
        inquire(unit=file%unit, size=file_size, iostat=status)
        call verify(status == 0, file, "query shard size", status)
        call verify(file%directory >= table_end(file) .and. file%directory <= file_size-footer_bytes+1, &
            file, "truncated shard or invalid directory position")
        remaining = file_size-footer_bytes-file%directory+1
        call verify(mod(remaining,mcp_group_directory_bytes) == 0 .and. &
            file%group_count == remaining/mcp_group_directory_bytes, file, "truncated shard or invalid group directory length")
        call allocate_block_tables(file)
        expected_group = 1
        expected_position = table_end(file)
        do block = 1, file%block_count
            read(file%unit, iostat=status) file%block_labels(block), file%block_sizes(block), file%cutoffs(block), &
                file%sparse_positions(block), file%column_counts(block), file%row_counts(block), &
                file%first_groups(:,block), file%group_counts(:,block), file%channel_counts(:,block)
            call verify(status == 0, file, "read block directory", status)
            call verify(file%block_sizes(block) > 0 .and. file%block_sizes(block) < huge(0) .and. &
                file%cutoffs(block) >= 0 .and. file%cutoffs(block) <= file%block_sizes(block), file, "invalid block metadata")
            call verify(file%column_counts(block) == indexed_column_count(int(file%block_sizes(block)),mpi_rank,mpi_process_count) &
                .and. file%row_counts(block) >= file%column_counts(block) .and. file%row_counts(block) < huge(0) &
                .and. file%sparse_positions(block) == expected_position, file, "invalid sparse block dimensions or offset")
            call verify(file%column_counts(block)+file%row_counts(block) <= (file%directory-expected_position)/4, &
                file, "sparse block exceeds data region")
            expected_position = expected_position + 4*(file%column_counts(block)+file%row_counts(block))
            do channel = 1, file%max_rank+2
                call verify(file%first_groups(channel,block) == expected_group .and. &
                    file%group_counts(channel,block) >= 0 .and. &
                    file%group_counts(channel,block) <= file%group_count-expected_group+1 .and. &
                    file%channel_counts(channel,block) >= file%group_counts(channel,block) .and. &
                    file%channel_counts(channel,block) <= (file%directory-expected_position)/contribution_bytes, &
                    file, "invalid channel dimensions or offsets")
                expected_group = expected_group + file%group_counts(channel,block)
                expected_position = expected_position + contribution_bytes*file%channel_counts(channel,block)
            end do
        end do
        call verify(expected_group == file%group_count+1 .and. expected_position == file%directory .and. &
            sum(file%block_sizes) == file%csf_count, file, "directory totals do not match shard dimensions")
        read(file%unit, pos=file_size-footer_bytes+1, iostat=status) tag, footer
        call verify(status == 0, file, "read completion footer", status)
        call verify(tag == complete_magic .and. all(footer == [file%generation,file%directory,file%group_count]), &
            file, "incomplete shard or invalid footer")
        identity = [file%generation,fields(3:)]
        call mpi_bcast(identity,size(identity),mpi_integer8,0,mpi_comm_world,mpi_error_code)
        call verify(identity(1) == file%generation, file, "shards generated by different runs")
        call verify(all(identity(2:) == fields(3:)), file, "inconsistent metadata across shards")
        root_basis_id = file%basis_id
        call mpi_bcast(root_basis_id,len(root_basis_id),mpi_character,0,mpi_comm_world,mpi_error_code)
        call verify(root_basis_id == file%basis_id, file, "inconsistent basis identifier across shards")
        allocate(root_partitions(2,file%block_count),root_labels(file%block_count),stat=status)
        call verify(status == 0, file, "allocate shard identity check", status)
        root_partitions(1,:) = file%block_sizes
        root_partitions(2,:) = file%cutoffs
        root_labels = file%block_labels
        call mpi_bcast(root_partitions,size(root_partitions),mpi_integer8,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast(root_labels,8*size(root_labels),mpi_character,0,mpi_comm_world,mpi_error_code)
        call verify(all(root_partitions(1,:) == file%block_sizes) .and. all(root_partitions(2,:) == file%cutoffs) .and. &
            all(root_labels == file%block_labels), file, "inconsistent block partitions across shards")
    end subroutine open_indexed_mcp_reader

    subroutine read_indexed_sparse(file, block, ends, rows)
        type(indexed_mcp), intent(in) :: file
        integer, intent(in) :: block
        integer, intent(out) :: ends(:), rows(:)
        integer :: status
        call verify(size(ends,kind=int64) == file%column_counts(block) .and. size(rows,kind=int64) == file%row_counts(block), &
            file, "invalid sparse destination dimensions")
        read(file%unit, pos=file%sparse_positions(block), iostat=status) ends, rows
        call verify(status == 0, file, "read sparse structure", status)
        call check_sparse_structure(ends,rows,file%unit,"MCP indexed")
    end subroutine read_indexed_sparse

    subroutine check_sparse_structure(ends, rows, unit, context)
        integer, intent(in) :: ends(:), rows(:), unit
        character(*), intent(in) :: context
        integer :: i, previous, column
        previous = 0
        do i = 1, size(ends)
            column = mpi_rank + 1 + (i-1)*mpi_process_count
            if (ends(i) <= previous .or. ends(i) > size(rows)) &
                call fail_input_io_mpi(0,unit,context//" invalid sparse column end")
            if (rows(ends(i)) /= column .or. any(rows(previous+1:ends(i)) < 1) .or. &
                any(rows(previous+1:ends(i)) > column)) call fail_input_io_mpi(0,unit,context//" invalid sparse row index")
            previous = ends(i)
        end do
        if (previous /= size(rows)) call fail_input_io_mpi(0,unit,context//" sparse column ends do not cover all rows")
    end subroutine check_sparse_structure

    subroutine read_indexed_directory(file, index, entries)
        type(indexed_mcp), intent(in) :: file
        integer(int64), intent(in) :: index
        type(mcp_group), intent(out) :: entries(:)
        integer :: status
        call verify(index >= 1 .and. index <= file%group_count .and. size(entries) <= file%group_count-index+1, &
            file, "group index outside directory")
        read(file%unit, pos=file%directory+mcp_group_directory_bytes*(index-1), iostat=status) entries
        call verify(status == 0, file, "read group directory entries", status)
    end subroutine read_indexed_directory

    subroutine read_indexed_group(file, index, position, available, previous_label, group, one_electron)
        type(indexed_mcp), intent(in) :: file
        integer(int64), intent(in) :: index, position, available, previous_label
        type(mcp_group), intent(out) :: group
        logical, intent(in) :: one_electron
        type(mcp_group) :: entry(1)
        call read_indexed_directory(file,index,entry)
        group = entry(1)
        call validate_indexed_group(file,position,available,previous_label,group,one_electron)
    end subroutine read_indexed_group

    subroutine validate_indexed_group(file, position, available, previous_label, group, one_electron)
        type(indexed_mcp), intent(in) :: file
        integer(int64), intent(in) :: position, available, previous_label
        type(mcp_group), intent(in) :: group
        logical, intent(in) :: one_electron
        integer :: orbitals(4)
        call verify(group%position == position .and. group%count > 0 .and. group%count <= available .and. &
            group%label > previous_label .and. group%label <= huge(0), file, "invalid group position, label, or count")
        call verify(group%bytes == contribution_bytes*group%count, file, "invalid group length")
        call decode_mcp_label(int(group%label),one_electron,orbitals)
        call verify(all(orbitals >= 1) .and. all(orbitals <= file%orbital_count), file, "orbital label outside basis")
    end subroutine validate_indexed_group

    subroutine read_indexed_contributions(file, group, first, columns, positions, coefficients)
        type(indexed_mcp), intent(in) :: file
        type(mcp_group), intent(in) :: group
        integer(int64), intent(in) :: first
        integer, intent(out) :: columns(:), positions(:)
        real(real64), intent(out) :: coefficients(:)
        integer :: status, i
        call verify(first >= 1 .and. size(columns) <= group%count-first+1 .and. &
            size(columns) == size(positions) .and. size(columns) == size(coefficients), file, "invalid contribution segment")
        read(file%unit, pos=group%position+contribution_bytes*(first-1), iostat=status) &
            (columns(i),positions(i),coefficients(i),i=1,size(columns))
        call verify(status == 0, file, "read contribution segment", status)
    end subroutine read_indexed_contributions

    pure subroutine decode_mcp_label(label, one_electron, orbitals)
        integer, intent(in) :: label
        logical, intent(in) :: one_electron
        integer, intent(out) :: orbitals(4)
        integer :: high
        if (one_electron) then
            orbitals(1:2) = [mod(label,orbital_index_radix),label/orbital_index_radix]
            orbitals(3:4) = orbitals(1:2)
        else
            high = label/orbital_index_radix
            orbitals(4) = mod(label,orbital_index_radix)
            orbitals(2) = mod(high,orbital_index_radix)
            high = high/orbital_index_radix
            orbitals(3) = mod(high,orbital_index_radix)
            orbitals(1) = high/orbital_index_radix
        end if
    end subroutine decode_mcp_label
end module mpi_indexed_mcp
