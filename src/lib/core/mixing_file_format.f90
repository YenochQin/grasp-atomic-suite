! Indexed little-endian ASF mixing streams.
module mixing_file_format
    use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real64, file_storage_size, iostat_end, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_support_datatype, ieee_is_finite
    use input_error_handling, only: fail_input_io, abort_input
    implicit none (type, external)
    private
    public :: mixing_file, open_mixing_input, next_mixing_block, read_mixing_levels, read_mixing_vectors
    public :: open_mixing_stream, open_transformed_mixing_stream, write_mixing_block, close_mixing_stream

    character(8), parameter :: stream_signature = "CSF--MIX"
    character(8), parameter :: incomplete_signature = "GRMIXTMP"
    integer(int64), parameter :: header_bytes = 112, block_bytes = 48, level_bytes = 16
    integer, parameter :: buffer_values = 131072

    type :: mixing_file
        integer :: unit = -1
        integer :: electron_count = 0, orbital_count = 0, csf_count = 0, level_count = 0, block_count = 0
        integer(int64) :: coefficient_count = 0
        character(64) :: basis_id = ""
        integer(int64), allocatable :: directory(:,:)
        integer :: block = 0, active_csfs = 0, active_levels = 0
    end type

    ! Each program writes one mixing stream at a time. MPI writers use this on rank zero.
    type(mixing_file), save :: output
    character(:), allocatable, save :: output_path
contains
    logical function supports_stream_encoding()
        supports_stream_encoding = file_storage_size == 8 .and. ieee_support_datatype(0.0_real64) .and. &
            all(transfer(1_int32, [0_int8]) == [1_int8,0_int8,0_int8,0_int8])
    end function

    subroutine open_mixing_stream(path, unit, status)
        use atomic_calculation_state, only: electron_count
        use orbital_csf_basis, only: orbital_count
        use csf_descriptor_reader, only: loaded_csf_basis_id
        use hamiltonian_symmetry_blocks, only: symmetry_block_count, csf_counts_by_block, eigenvector_counts_by_block
        use level_symmetry_metadata, only: packed_block_symmetries
        character(*), intent(in) :: path
        integer, intent(out) :: unit, status
        integer :: block, symmetry
        integer(int64) :: level_start, offset, csfs, levels, values, total_csfs, total_levels

        unit = -1
        status = 1
        if (output%unit /= -1 .or. symmetry_block_count < 1) return
        if (.not. supports_stream_encoding() .or. len_trim(loaded_csf_basis_id) /= 64) return
        total_csfs = sum(int(csf_counts_by_block(1:symmetry_block_count),int64))
        total_levels = sum(int(eigenvector_counts_by_block(:symmetry_block_count),int64))
        if (min(total_csfs,total_levels) < 1 .or. max(total_csfs,total_levels) > huge(0)) return
        if (allocated(output%directory)) deallocate(output%directory)
        allocate(output%directory(6,symmetry_block_count), stat=status)
        if (status /= 0) return
        output%electron_count = electron_count
        output%orbital_count = orbital_count
        output%basis_id = loaded_csf_basis_id
        output%block_count = symmetry_block_count
        output%csf_count = int(total_csfs)
        output%level_count = int(total_levels)
        offset = header_bytes + block_bytes*symmetry_block_count + level_bytes*output%level_count
        level_start = 0
        output%coefficient_count = 0
        do block = 1, symmetry_block_count
            symmetry = packed_block_symmetries(block)
            csfs = csf_counts_by_block(block)
            levels = eigenvector_counts_by_block(block)
            values = csfs*levels
            if (symmetry == 0 .or. csfs < 1 .or. levels < 1 .or. levels > csfs .or. &
                values > (huge(offset)-offset)/8) then
                status = 1
                return
            end if
            output%directory(:,block) = [int(abs(symmetry)-1,int64), int(sign(1,symmetry),int64), &
                csfs, levels, level_start, offset]
            level_start = level_start + levels
            offset = offset + 8*values
            output%coefficient_count = output%coefficient_count + values
        end do
        call open_mixing_output(path, unit, status)
        if (status == 0) output_path = path
    end subroutine

    subroutine open_transformed_mixing_stream(input_file, path, unit, status)
        type(mixing_file), intent(in) :: input_file
        character(*), intent(in) :: path
        integer, intent(out) :: unit, status

        unit = -1
        status = 1
        if (output%unit /= -1) return
        ! The transformation changes coefficients, not the source block layout or basis digest.
        output = input_file
        output%unit = -1
        output%block = 0
        call open_mixing_output(path, unit, status)
        if (status == 0) output_path = path
    end subroutine

    subroutine open_mixing_output(path, unit, status)
        character(*), intent(in) :: path
        integer, intent(out) :: unit, status
        integer :: close_status

        unit = -1
        open(newunit=unit, file=path, status="replace", action="write", access="stream", form="unformatted", iostat=status)
        if (status /= 0) return
        write(unit,iostat=status) incomplete_signature, transfer([int(output%electron_count,int64), &
            int(output%orbital_count,int64), int(output%block_count,int64), int(output%csf_count,int64), &
            int(output%level_count,int64)], [0_int8]), output%basis_id, transfer(output%directory,[0_int8])
        if (status == 0) flush(unit,iostat=status)
        if (status /= 0) then
            close(unit,iostat=close_status)
            unit = -1
            return
        end if
        output%unit = unit
    end subroutine

    subroutine write_mixing_block(unit, block, roots, energies, coefficients, abort_handler)
        integer, intent(inout) :: unit
        integer, intent(in) :: block, roots(:)
        real(real64), intent(in) :: energies(:), coefficients(:)
        procedure(abort_input), optional :: abort_handler
        integer(int64) :: position, first, last, actual_size
        integer :: root, status
        integer(int8) :: table(16,size(roots))

        if (unit == -1 .or. unit /= output%unit .or. block < 1 .or. block > output%block_count) &
            call fail_mixing_output(0, "Invalid mixing output unit or block", abort_handler)
        if (block == 1 .and. output%block == output%block_count) then
            ! Start the new pass with an empty payload. A short rewrite must not
            ! inherit bytes from the preceding pass and satisfy the size check.
            close(unit,iostat=status)
            if (status /= 0) call fail_mixing_output(status, "Cannot close preceding mixing pass", abort_handler)
            call open_mixing_output(output_path,unit,status)
            if (status /= 0) call fail_mixing_output(status, "Cannot start new mixing pass", abort_handler)
            output%block = 0
        end if
        if (block /= output%block+1) call fail_mixing_output(0, "Out-of-order mixing block", abort_handler)
        if (size(roots) /= output%directory(4,block) .or. size(energies) /= size(roots) .or. &
            size(coefficients,kind=int64) /= output%directory(3,block)*output%directory(4,block)) &
            call fail_mixing_output(0, "Mixing block dimensions differ from the directory", abort_handler)
        do root = 1, size(roots)
            table(:8,root) = transfer(int(roots(root),int64), [0_int8])
            table(9:,root) = transfer(energies(root), [0_int8])
        end do
        if (block == 1) then
            write(unit,pos=1,iostat=status) incomplete_signature
            if (status /= 0) call fail_mixing_output(status, "Cannot begin mixing coefficient update", abort_handler)
            flush(unit,iostat=status)
            if (status /= 0) call fail_mixing_output(status, "Cannot flush incomplete mixing signature", abort_handler)
        end if
        position = header_bytes + block_bytes*output%block_count + level_bytes*output%directory(5,block) + 1
        write(unit,pos=position,iostat=status) table
        if (status /= 0) call fail_mixing_output(status, "Cannot write ASF energies", abort_handler)
        position = output%directory(6,block) + 1
        do first = 1_int64, size(coefficients,kind=int64), buffer_values
            last = min(first+buffer_values-1, size(coefficients,kind=int64))
            write(unit,pos=position,iostat=status) transfer(coefficients(first:last), [0_int8])
            if (status /= 0) call fail_mixing_output(status, "Cannot write ASF mixing coefficients", abort_handler)
            position = position + 8*(last-first+1)
        end do
        if (block == output%block_count) then
            flush(unit,iostat=status)
            if (status /= 0) call fail_mixing_output(status, "Cannot flush mixing payload", abort_handler)
            ! Check the disconnected file's size. A runtime may report the logical
            ! position on a connected unit even if the underlying write was short.
            close(unit,iostat=status)
            if (status /= 0) call fail_mixing_output(status, "Cannot close the mixing payload", abort_handler)
            actual_size = -1_int64
            inquire(file=output_path,size=actual_size,iostat=status)
            if (status /= 0 .or. actual_size /= position-1) &
                call fail_mixing_output(status, "Incomplete mixing payload after flush", abort_handler)
            open(newunit=unit,file=output_path,status="old",action="write",access="stream",form="unformatted",iostat=status)
            if (status /= 0) call fail_mixing_output(status, "Cannot reopen mixing payload", abort_handler)
            output%unit = unit
            write(unit,pos=1,iostat=status) stream_signature
            if (status /= 0) call fail_mixing_output(status, "Cannot finalize mixing stream", abort_handler)
            flush(unit,iostat=status)
            if (status /= 0) call fail_mixing_output(status, "Cannot flush complete mixing signature", abort_handler)
        end if
        output%block = block
    end subroutine

    subroutine close_mixing_stream(unit, abort_handler)
        integer, intent(in) :: unit
        procedure(abort_input), optional :: abort_handler
        integer :: status

        if (unit == -1 .or. unit /= output%unit) &
            call fail_mixing_output(0, "Invalid mixing output unit on close", abort_handler)
        if (output%block /= output%block_count) &
            call fail_mixing_output(0, "Cannot close an incomplete mixing stream", abort_handler)
        close(unit,iostat=status)
        if (status /= 0) call fail_mixing_output(status, "Cannot close mixing stream", abort_handler)
        output = mixing_file()
        deallocate(output_path)
    end subroutine close_mixing_stream

    subroutine fail_mixing_output(status, context, abort_handler)
        integer, intent(in) :: status
        character(*), intent(in) :: context
        procedure(abort_input), optional :: abort_handler

        write(error_unit,"(a,a,a,i0,a)") "Mixing output failure: ", context, " (iostat ", status, ")"
        if (allocated(output_path)) write(error_unit,"(a,a)") "Mixing output file: ", output_path
        flush(error_unit)
        if (present(abort_handler)) call abort_handler()
        error stop 1
    end subroutine fail_mixing_output

    subroutine open_mixing_input(file, path, basis_id)
        type(mixing_file), intent(out) :: file
        character(*), intent(in) :: path
        character(*), intent(in), optional :: basis_id
        character(8) :: signature
        integer :: status, block
        integer(int64) :: fields(5), file_size, offset, levels, csfs, values
        integer(int8) :: header(40)
        integer(int8), allocatable :: directory(:)

        open(newunit=file%unit,file=path,status="old",action="read",access="stream",form="unformatted",iostat=status)
        if (status /= 0) call fail_input_io(status,-1,"cannot open mixing file: "//trim(path))
        read(file%unit,iostat=status) signature
        if (status /= 0) call fail_input_io(status,file%unit,"Cannot read mixing stream signature: "//trim(path))
        if (signature == incomplete_signature) call fail_input_io(0,file%unit,"Incomplete mixing stream")
        if (signature /= stream_signature) &
            call fail_input_io(0,file%unit,"Expected CSF--MIX signature: "//trim(path))
        if (.not. supports_stream_encoding()) call fail_input_io(0,file%unit,"Mixing stream requires little-endian IEEE binary64")
        read(file%unit,iostat=status) header, file%basis_id
        if (status /= 0) call fail_input_io(status,file%unit,"Cannot read mixing stream header")
        fields = transfer(header, fields)
        if (minval(fields) < 1 .or. maxval(fields) > huge(0) .or. fields(5) > fields(4) .or. fields(3) > fields(5)) &
            call fail_input_io(0,file%unit,"Invalid mixing stream counts")
        if (verify(file%basis_id,"0123456789abcdef") /= 0) call fail_input_io(0,file%unit,"Invalid mixing basis digest")
        if (present(basis_id)) then
            if (basis_id /= "" .and. basis_id /= file%basis_id) &
                call fail_input_io(0,file%unit,"Mixing file does not match the CSF descriptor SHA-256")
        end if
        file%electron_count = int(fields(1))
        file%orbital_count = int(fields(2))
        file%block_count = int(fields(3))
        file%csf_count = int(fields(4))
        file%level_count = int(fields(5))
        inquire(unit=file%unit,size=file_size,iostat=status)
        if (status /= 0) call fail_input_io(status,file%unit,"Cannot determine mixing stream size")
        offset = header_bytes + block_bytes*fields(3) + level_bytes*fields(5)
        if (offset > file_size) call fail_input_io(0,file%unit,"Truncated mixing block directory or level table")
        allocate(directory(block_bytes*fields(3)),file%directory(6,file%block_count),stat=status)
        if (status /= 0) call fail_input_io(status,file%unit,"Cannot allocate the mixing block directory")
        read(file%unit,iostat=status) directory
        if (status /= 0) call fail_input_io(status,file%unit,"Cannot read mixing block directory")
        file%directory = reshape(transfer(directory,[0_int64]),[6,file%block_count])
        levels = 0
        csfs = 0
        do block = 1, file%block_count
            associate(row => file%directory(:,block))
                if (row(1) < 0 .or. row(1) >= huge(0) .or. (row(2) /= -1 .and. row(2) /= 1) .or. &
                    row(3) < 1 .or. row(3) > file%csf_count-csfs .or. &
                    row(4) < 1 .or. row(4) > min(row(3),file%level_count-levels) .or. &
                    row(5) /= levels .or. row(6) /= offset) &
                    call fail_input_io(0,file%unit,"Invalid mixing block directory entry")
                values = row(3)*row(4)
                if (values > (file_size-offset)/8) call fail_input_io(0,file%unit,"Truncated mixing coefficient matrix")
                csfs = csfs + row(3)
                levels = levels + row(4)
                offset = offset + 8*values
                file%coefficient_count = file%coefficient_count + values
            end associate
        end do
        if (offset /= file_size .or. levels /= file%level_count .or. csfs /= file%csf_count) &
            call fail_input_io(0,file%unit,"Mixing stream totals disagree with the directory")
    end subroutine

    subroutine next_mixing_block(file, block, csfs, levels, j_code, parity, status)
        type(mixing_file), intent(inout) :: file
        integer, intent(out) :: block, csfs, levels, j_code, parity
        integer, optional, intent(out) :: status
        integer :: io_status

        io_status = 0
        if (file%block == file%block_count) then
            io_status = iostat_end
        else
            file%block = file%block + 1
            block = file%block
            csfs = int(file%directory(3,block))
            levels = int(file%directory(4,block))
            j_code = int(file%directory(1,block))+1
            parity = int(file%directory(2,block))
        end if
        if (present(status)) then
            status = io_status
        else if (io_status /= 0) then
            call fail_input_io(io_status,file%unit,"Cannot read mixing block")
        end if
        if (io_status == 0) then
            file%active_csfs = csfs
            file%active_levels = levels
        end if
    end subroutine

    subroutine read_mixing_levels(file, roots, shift, energies)
        type(mixing_file), intent(in) :: file
        integer, intent(out) :: roots(:)
        real(real64), intent(out) :: shift, energies(:)
        integer :: root, status
        integer(int64) :: position, root_id
        integer(int8) :: table(16,size(roots))

        position = header_bytes + block_bytes*file%block_count + level_bytes*file%directory(5,file%block)+1
        read(file%unit,pos=position,iostat=status) table
        if (status /= 0) call fail_input_io(status,file%unit,"Cannot read mixing stream levels")
        shift = 0.0_real64
        do root = 1, size(roots)
            root_id = transfer(table(:8,root),root_id)
            if (root_id < 1 .or. root_id > file%active_csfs) call fail_input_io(0,file%unit,"Invalid ASF root index")
            roots(root) = int(root_id)
            energies(root) = transfer(table(9:,root),shift)
        end do
        if (.not. ieee_is_finite(shift) .or. .not. all(ieee_is_finite(energies))) &
            call fail_input_io(0,file%unit,"Nonfinite ASF energy")
    end subroutine

    ! stride permits a block to be read into an existing global-CSF matrix.
    subroutine read_mixing_vectors(file, coefficients, stride)
        type(mixing_file), intent(in) :: file
        real(real64), intent(inout) :: coefficients(:)
        integer, intent(in) :: stride
        integer :: status, level, count
        integer(int64) :: row, offset, first, last
        integer(int8), allocatable :: buffer(:)

        allocate(buffer(8*min(buffer_values,file%active_csfs)))
        offset = file%directory(6,file%block)+1
        do level = 1, file%active_levels
            row = int(level-1,int64)*stride
            do first = 1, file%active_csfs, buffer_values
                last = min(first+buffer_values-1,int(file%active_csfs,int64))
                count = int(last-first+1)
                read(file%unit,pos=offset,iostat=status) buffer(:8*count)
                if (status /= 0) call fail_input_io(status,file%unit,"Cannot read mixing coefficient matrix")
                coefficients(row+first:row+last) = transfer(buffer(:8*count),[0.0_real64])
                offset = offset + 8_int64*count
            end do
        end do
    end subroutine
end module mixing_file_format
