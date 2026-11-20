! Copy the requested contiguous prefix of an integer vector.
module integer_vector_copy
    implicit none (type, external)
    private
    public :: icopy
contains
      subroutine icopy(element_count, source_values, source_stride, destination_values, destination_stride)
      implicit none (type, external)
      integer, intent(in) :: element_count
      integer, intent(in) :: source_stride
      integer, intent(in) :: destination_stride
      integer, intent(in) :: source_values(*)
      integer, intent(out) :: destination_values(*)
      integer :: i, prefix_length, first_full_block_index
      ! The stride arguments are accepted for compatibility. This routine copies contiguous elements.
      if (element_count <= 0) return
      prefix_length = mod(element_count,7)
      if (prefix_length /= 0) then
         destination_values(:prefix_length) = source_values(:prefix_length)
         if (element_count < 7) return
      end if
      first_full_block_index = prefix_length + 1
      destination_values(first_full_block_index:((element_count-first_full_block_index+7)/7)*7-1+first_full_block_index) = &
         & source_values(first_full_block_index:((element_count-first_full_block_index+7)/7)*7-1+first_full_block_index)
      return
      end subroutine icopy
end module integer_vector_copy
