! Declare the C ABI v3 used to call the Rust CSF descriptor reader.
module csf_descriptor_c_bindings
    ! Module dependencies
    use iso_c_binding, only: c_ptr, c_char, c_int32_t, c_int8_t, c_int64_t, c_size_t
    implicit none (type, external)

    ! Module interface
    private
    public :: rust_descriptor_abi_version, rust_descriptor_open, rust_descriptor_info
    public :: rust_descriptor_tables, rust_descriptor_read, rust_descriptor_close
    public :: rust_descriptor_zero_order_counts, rust_descriptor_format_version
    public :: rust_descriptor_read_radial_header, rust_descriptor_basis_id
    public :: rust_descriptor_select_block, rust_descriptor_check_spin_angular

    ! Procedure interfaces
    interface
        subroutine rust_descriptor_basis_id(handle, digest) bind(C, name="rcsfs_descriptor_basis_id")
            import c_ptr, c_char
            type(c_ptr), value :: handle
            character(kind=c_char), intent(out) :: digest(*)
        end subroutine rust_descriptor_basis_id
        function rust_descriptor_abi_version() bind(C, name="rcsfs_descriptor_abi_version") result(status)
            ! Module dependencies
            import c_int32_t

            ! Arguments
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_format_version() bind(C, name="rcsfs_descriptor_format_version") result(version)
            import c_int32_t
            integer(c_int32_t) :: version
        end function
        function rust_descriptor_read_radial_header(header_path, orbital_table, orbital_capacity, header_counts, &
            & error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_read_radial_header") result(status)
            ! Module dependencies
            import c_char, c_size_t, c_int32_t

            ! Arguments
            character(kind=c_char), intent(in) :: header_path(*)
            integer(c_int32_t), intent(out) :: orbital_table(*), header_counts(3)
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_size_t), value :: orbital_capacity, message_capacity
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_open(descriptor_path, header_path, handle, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_open") result(status)
            ! Module dependencies
            import c_ptr, c_char, c_size_t, c_int32_t

            ! Arguments
            character(kind=c_char), intent(in) :: descriptor_path(*), header_path(*)
            type(c_ptr), intent(inout) :: handle
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_size_t), value :: message_capacity
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_check_spin_angular(handle, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_check_spin_angular") result(status)
            import c_ptr, c_char, c_size_t, c_int32_t
            type(c_ptr), value :: handle
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_size_t), value :: message_capacity
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_info(handle, header_counts, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_info") result(status)
            ! Module dependencies
            import c_ptr, c_char, c_size_t, c_int32_t, c_int64_t

            ! Arguments
            type(c_ptr), value :: handle
            integer(c_int64_t), intent(out) :: header_counts(5)
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_size_t), value :: message_capacity
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_tables(handle, orbital_table, orbital_capacity, block_lengths, symmetry_table, &
            & block_capacity, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_tables") result(status)
            ! Module dependencies
            import c_ptr, c_char, c_size_t, c_int32_t, c_int64_t

            ! Arguments
            type(c_ptr), value :: handle
            integer(c_int32_t), intent(out) :: orbital_table(*), symmetry_table(*)
            integer(c_int64_t), intent(out) :: block_lengths(*)
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_size_t), value :: message_capacity, orbital_capacity, block_capacity
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_zero_order_counts(handle, zero_order_counts, block_capacity, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_zero_order_counts") result(status)
            ! Module dependencies
            import c_ptr, c_char, c_size_t, c_int32_t, c_int64_t

            ! Arguments
            type(c_ptr), value :: handle
            integer(c_int64_t), intent(out) :: zero_order_counts(*)
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_size_t), value :: message_capacity, block_capacity
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_select_block(handle, block, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_select_block") result(status)
            import c_ptr, c_char, c_size_t, c_int32_t
            type(c_ptr), value :: handle
            integer(c_size_t), value :: block, message_capacity
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_int32_t) :: status
        end function
        function rust_descriptor_read(handle, row_values, value_capacity, row_count, error_message, message_capacity) &
            bind(C, name="rcsfs_descriptor_read") result(status)
            ! Module dependencies
            import c_ptr, c_char, c_size_t, c_int32_t, c_int8_t

            ! Arguments
            type(c_ptr), value :: handle
            integer(c_int8_t), intent(out) :: row_values(*)
            integer(c_size_t), value :: value_capacity, row_count, message_capacity
            character(kind=c_char), intent(out) :: error_message(*)
            integer(c_int32_t) :: status
        end function
        subroutine rust_descriptor_close(handle) bind(C, name="rcsfs_descriptor_close")
            ! Module dependencies
            import c_ptr

            ! Arguments
            type(c_ptr), intent(inout) :: handle
        end subroutine
    end interface
end module csf_descriptor_c_bindings
