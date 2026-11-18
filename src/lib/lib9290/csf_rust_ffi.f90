! C ABI v1 for the native rCSFs descriptor reader. No MPI or global state.
module csf_rust_ffi
   use iso_c_binding, only: c_ptr, c_char, c_int32_t, c_int8_t, c_int64_t, c_size_t
   implicit none (type, external)
   private
   public :: rust_descriptor_version, rust_descriptor_open, rust_descriptor_info
   public :: rust_descriptor_tables, rust_descriptor_read, rust_descriptor_close
   public :: rust_descriptor_zero_order_counts
   public :: rust_descriptor_read_radial_header
   interface
      function rust_descriptor_version() bind(C, name="rcsfs_descriptor_abi_version") result(status)
         import c_int32_t
         integer(c_int32_t) :: status
      end function
      function rust_descriptor_read_radial_header(header, orbitals, orbital_capacity, info, message, capacity) &
         bind(C, name="rcsfs_descriptor_read_radial_header") result(status)
         import c_char, c_size_t, c_int32_t
         character(kind=c_char), intent(in) :: header(*)
         integer(c_int32_t), intent(out) :: orbitals(*), info(3)
         character(kind=c_char), intent(out) :: message(*)
         integer(c_size_t), value :: orbital_capacity, capacity
         integer(c_int32_t) :: status
      end function
      function rust_descriptor_open(path, header, handle, message, capacity) &
         bind(C, name="rcsfs_descriptor_open") result(status)
         import c_ptr, c_char, c_size_t, c_int32_t
         character(kind=c_char), intent(in) :: path(*), header(*)
         type(c_ptr), intent(inout) :: handle
         character(kind=c_char), intent(out) :: message(*)
         integer(c_size_t), value :: capacity
         integer(c_int32_t) :: status
      end function
      function rust_descriptor_info(handle, info, message, capacity) &
         bind(C, name="rcsfs_descriptor_info") result(status)
         import c_ptr, c_char, c_size_t, c_int32_t, c_int64_t
         type(c_ptr), value :: handle
         integer(c_int64_t), intent(out) :: info(5)
         character(kind=c_char), intent(out) :: message(*)
         integer(c_size_t), value :: capacity
         integer(c_int32_t) :: status
      end function
      function rust_descriptor_tables(handle, orbitals, orbital_capacity, lengths, symmetries, block_capacity, message, capacity) &
         bind(C, name="rcsfs_descriptor_tables") result(status)
         import c_ptr, c_char, c_size_t, c_int32_t, c_int64_t
         type(c_ptr), value :: handle
         integer(c_int32_t), intent(out) :: orbitals(*), symmetries(*)
         integer(c_int64_t), intent(out) :: lengths(*)
         character(kind=c_char), intent(out) :: message(*)
         integer(c_size_t), value :: capacity, orbital_capacity, block_capacity
         integer(c_int32_t) :: status
      end function
      function rust_descriptor_zero_order_counts(handle, counts, block_capacity, message, capacity) &
         bind(C, name="rcsfs_descriptor_zero_order_counts") result(status)
         import c_ptr, c_char, c_size_t, c_int32_t, c_int64_t
         type(c_ptr), value :: handle
         integer(c_int64_t), intent(out) :: counts(*)
         character(kind=c_char), intent(out) :: message(*)
         integer(c_size_t), value :: capacity, block_capacity
         integer(c_int32_t) :: status
      end function
      function rust_descriptor_read(handle, values, size, rows, message, capacity) &
         bind(C, name="rcsfs_descriptor_read") result(status)
         import c_ptr, c_char, c_size_t, c_int32_t, c_int8_t
         type(c_ptr), value :: handle
         integer(c_int8_t), intent(out) :: values(*)
         integer(c_size_t), value :: size, rows, capacity
         character(kind=c_char), intent(out) :: message(*)
         integer(c_int32_t) :: status
      end function
      subroutine rust_descriptor_close(handle) bind(C, name="rcsfs_descriptor_close")
         import c_ptr
         type(c_ptr), intent(inout) :: handle
      end subroutine
   end interface
end module csf_rust_ffi
