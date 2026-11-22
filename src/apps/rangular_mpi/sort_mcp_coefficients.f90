module sort_mcp_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: sort_mcp_coefficients
contains

    ! This routine sorts lists
    !
    !                   (ICLMN,INDEX,LABEL,COEFF)
    ! into lists
    !                   (LABEL,ICLMN,INDEX,COEFF)
    !
    ! using Heapsort. File NFILE is closed by this routine.  NCOEFF is
    ! the number of triads (INDEX, ...) and is an input. NTGRAL is the
    ! number of different values of LABEL, and is an output. If LPRINT
    ! is  .TRUE. , the contents of the sorted file are interpreted and
    ! printed to unit 99.
    !
    ! Written by Farid A. Parpia. Last revision: 21 Dec 1992
    subroutine sort_mcp_coefficients(mcp_channel_index, coefficient_count, distinct_integral_count, print_details, &
        & symmetry_block_index, &
        indexed_output)
        ! Module dependencies
        use mpi_indexed_mcp, only: indexed_mcp, start_indexed_channel, write_indexed_group
        use mcp_generation_options, only: mcp_temporary_units
        use, intrinsic :: iso_fortran_env, only: output_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc
        use orbital_csf_basis,           only: orbital_principal_n, orbital_symmetry_labels
        use standard_io_units, only: stderr_unit
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Local declarations
        integer :: heap_parent_index, heap_child_index

        ! Arguments
        integer, intent(in) :: mcp_channel_index
        integer, intent(in) :: coefficient_count
        integer, intent(out)  :: distinct_integral_count
        integer, intent(in)   :: symmetry_block_index
        logical, intent(in)   :: print_details
        type(indexed_mcp), intent(inout) :: indexed_output

        ! Local constants

        ! Local declarations
        integer, pointer :: csf_column_indices(:), matrix_element_positions(:), packed_orbital_labels(:)
        real(real64), pointer :: angular_coefficients(:)
        integer :: text_length, channel_number_text_length, file_block_index, i, heap_root_index, heap_end_index, &
            & saved_csf_column, &
            & saved_matrix_position, saved_orbital_label, &
            & j,&
            group_orbital_label, first_group_index, last_group_index, contribution_count, ia, ib, k, id, ic
        real(real64) :: saved_coefficient
        character :: number_text*20, mcp_channel_label*3, channel_number_text*2

        call format_integer (coefficient_count, number_text, text_length)
        if (mcp_channel_index > 31) then
            call format_integer (mcp_channel_index - 32, channel_number_text, channel_number_text_length)
            write (output_unit, *) &
                & "Sorting "//number_text(1:text_length)//" V(k="//channel_number_text(1:channel_number_text_length)//&
                ") coefficients ...", mcp_channel_index
        else
            write (output_unit, *) "Sorting "//number_text(1:text_length)//" T coefficients ...", mcp_channel_index
        end if
        call start_indexed_channel(indexed_output,mcp_channel_index-30)
        read (mcp_temporary_units(mcp_channel_index)) mcp_channel_label, file_block_index
        if (symmetry_block_index /= file_block_index) then
            write (stderr_unit, *) "sort_mcp_coefficients: nb = ", symmetry_block_index, ".NE. mb (=", file_block_index, ")"
            stop
        end if
        ! Sort the list
        if (coefficient_count > 0) then
            ! Allocate storage for all required arrays
            call alloc (angular_coefficients, coefficient_count, "angular_coefficients", "sort_mcp_coefficients")
            call alloc (csf_column_indices, coefficient_count, "csf_column_indices", "sort_mcp_coefficients")
            call alloc (matrix_element_positions, coefficient_count, "matrix_element_positions", "sort_mcp_coefficients")
            call alloc (packed_orbital_labels, coefficient_count, "packed_orbital_labels", "sort_mcp_coefficients")
            ! Read arrays into memory from NFILE
            do i = 1, coefficient_count
                read (mcp_temporary_units(mcp_channel_index)) csf_column_indices(i), matrix_element_positions(i), &
                    & packed_orbital_labels(i), angular_coefficients(i)
            end do
        end if
        ! Sort LABEL into ascending order using the heapsort algorithm;
        ! move the associated members of COEFF and INDEX in the same
        ! manner; the code below is adapted from Press et al.
        if (coefficient_count > 1) then
            heap_root_index = coefficient_count/2 + 1
            heap_end_index = coefficient_count
            do
                if (heap_root_index > 1) then
                    heap_root_index = heap_root_index - 1
                    saved_coefficient = angular_coefficients(heap_root_index)
                    saved_csf_column = csf_column_indices(heap_root_index)
                    saved_matrix_position = matrix_element_positions(heap_root_index)
                    saved_orbital_label = packed_orbital_labels(heap_root_index)
                else
                    saved_coefficient = angular_coefficients(heap_end_index)
                    saved_csf_column = csf_column_indices(heap_end_index)
                    saved_matrix_position = matrix_element_positions(heap_end_index)
                    saved_orbital_label = packed_orbital_labels(heap_end_index)
                    angular_coefficients(heap_end_index) = angular_coefficients(1)
                    csf_column_indices(heap_end_index) = csf_column_indices(1)
                    matrix_element_positions(heap_end_index) = matrix_element_positions(1)
                    packed_orbital_labels(heap_end_index) = packed_orbital_labels(1)
                    heap_end_index = heap_end_index - 1
                    if (heap_end_index == 1) then
                        angular_coefficients(1) = saved_coefficient
                        csf_column_indices(1) = saved_csf_column
                        matrix_element_positions(1) = saved_matrix_position
                        packed_orbital_labels(1) = saved_orbital_label
                        exit
                    end if
                end if
                heap_parent_index = heap_root_index
                heap_child_index = heap_root_index + heap_root_index
                do while (heap_child_index <= heap_end_index)
                    if (heap_child_index < heap_end_index) then
                        if (packed_orbital_labels(heap_child_index) < packed_orbital_labels(heap_child_index+1)) &
                            & heap_child_index = heap_child_index + 1
                    end if
                    if (saved_orbital_label < packed_orbital_labels(heap_child_index)) then
                        angular_coefficients(heap_parent_index) = angular_coefficients(heap_child_index)
                        csf_column_indices(heap_parent_index) = csf_column_indices(heap_child_index)
                        matrix_element_positions(heap_parent_index) = matrix_element_positions(heap_child_index)
                        packed_orbital_labels(heap_parent_index) = packed_orbital_labels(heap_child_index)
                        heap_parent_index = heap_child_index
                        if (heap_parent_index > heap_end_index/2) exit
                        heap_child_index = 2*heap_parent_index
                    else
                        exit
                    end if
                end do
                angular_coefficients(heap_parent_index) = saved_coefficient
                csf_column_indices(heap_parent_index) = saved_csf_column
                matrix_element_positions(heap_parent_index) = saved_matrix_position
                packed_orbital_labels(heap_parent_index) = saved_orbital_label
            end do

        end if
        first_group_index = 1
        distinct_integral_count = 0
        do i = 2, coefficient_count
            if (packed_orbital_labels(i) == packed_orbital_labels(first_group_index)) cycle
            call write_group(first_group_index,i-1)
            first_group_index = i
        end do
        if (coefficient_count > 0) call write_group(first_group_index,coefficient_count)
        ! Completion message
        if (print_details) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From sort_mcp_coefficients:")')
            write (output_unit, '(/,"From sort_mcp_coefficients:")')
            if (coefficient_count > 0) then
                group_orbital_label = packed_orbital_labels(1)
                first_group_index = 1
                last_group_index = 1
                do i = 2, coefficient_count
                    if (packed_orbital_labels(i) == group_orbital_label) last_group_index = last_group_index + 1
                    do
                        if (packed_orbital_labels(i)==group_orbital_label .and. i/=coefficient_count) exit
                        saved_orbital_label = group_orbital_label
                        contribution_count = last_group_index - first_group_index + 1
                        if (mcp_channel_index == 31) then
                            ia = mod(saved_orbital_label,orbital_index_radix)
                            ib = saved_orbital_label/orbital_index_radix
                            call ensure_debug_output()
                            write (angular_debug_unit, '(" I(",1i2,1a2,",",1i2,1a2,"):")') &
                                orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ib), &
                                & orbital_symmetry_labels(ib)
                            do j = first_group_index, last_group_index
                                call ensure_debug_output()
                                write (angular_debug_unit, '("  T_[",1i2,",",1i4,"] = ",1p,d19.12)') csf_column_indices(j), &
                                & matrix_element_positions(j), angular_coefficients(j)
                            end do
                        else
                            k = mcp_channel_index - 32
                            id = mod(saved_orbital_label,orbital_index_radix)
                            saved_orbital_label = saved_orbital_label/orbital_index_radix
                            ib = mod(saved_orbital_label,orbital_index_radix)
                            saved_orbital_label = saved_orbital_label/orbital_index_radix
                            ic = mod(saved_orbital_label,orbital_index_radix)
                            ia = saved_orbital_label/orbital_index_radix
                            call ensure_debug_output()
                            write (angular_debug_unit, '(" R^[(",1i2,")] (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,1a2,"):")') &
                                k,orbital_principal_n(ia),orbital_symmetry_labels(ia),orbital_principal_n(ib), &
                                & orbital_symmetry_labels(ib),orbital_principal_n(ic), &
                                orbital_symmetry_labels(ic),orbital_principal_n(id),orbital_symmetry_labels(id)
                            do j = first_group_index, last_group_index
                                call ensure_debug_output()
                                write (angular_debug_unit, &
                                & '("  V^[(",1i2,")]_[",1i8,",",1i8,"] = ",1p,d19.12)') &
                                k, csf_column_indices(j), matrix_element_positions( &
                                    & j),angular_coefficients(j)
                            end do
                        end if
                        group_orbital_label = packed_orbital_labels(i)
                        first_group_index = last_group_index + 1
                        last_group_index = first_group_index
                        if (last_group_index /= coefficient_count) exit
                    end do
                end do
            end if
            call ensure_debug_output()
            write (angular_debug_unit, '("  Number of integrals is ",1i4)') distinct_integral_count
        end if
        ! Deallocate storage
        if (coefficient_count > 0) then
            call dalloc (angular_coefficients,"angular_coefficients", "sort_mcp_coefficients")
            call dalloc (csf_column_indices, "csf_column_indices", "sort_mcp_coefficients")
            call dalloc (matrix_element_positions, "matrix_element_positions", "sort_mcp_coefficients")
            call dalloc (packed_orbital_labels, "packed_orbital_labels", "sort_mcp_coefficients")
        end if

    contains

        subroutine write_group(first, last)
            integer, intent(in) :: first, last
            call write_indexed_group(indexed_output,packed_orbital_labels(first),csf_column_indices(first:last), &
                matrix_element_positions(first:last),angular_coefficients(first:last))
            distinct_integral_count = distinct_integral_count + 1
        end subroutine write_group
    end subroutine sort_mcp_coefficients
end module sort_mcp_coefficients_interface
