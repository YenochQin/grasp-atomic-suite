! Read the requested interaction partition from the bound descriptor header.
subroutine getinf(source, status)
   use csf_descriptor_input, only: descriptor_source, read_descriptor_zero_order_counts
   use default_c, only: ndef
   use hblock_c, only: symmetry_block_count
   use iccu_c, only: iccut
   use mcp_c, only: diagonal_only, use_first_order_restriction
   implicit none (type, external)
   type(descriptor_source), intent(in) :: source
   integer, intent(out) :: status

   status = 0
   diagonal_only = .false.
   use_first_order_restriction = ndef /= 0
   iccut = 0
   if (use_first_order_restriction) then
      call read_descriptor_zero_order_counts(source, iccut(:symmetry_block_count), status)
   end if
end subroutine getinf
