! A sequential block reader for the RCI matrix construction lifecycle.
module csf_block_input
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                  install_descriptor_header, read_descriptor_block, read_descriptor_zero_order_counts
   use rmcdhf_descriptor_input, only: descriptor_path
   use hblock_c, only: csf_counts_by_block
   use orb_c, only: csf_subshell_occupations
   use stat_c, only: csf_subshell_quantum_codes, csf_coupling_quantum_codes
   implicit none (type, external)
   private
   type(descriptor_source), save :: source
   public :: open_csf_blocks, load_csf_block
   public :: read_csf_zero_order_counts
contains
   subroutine open_csf_blocks(name,ncore,labels,status)
      character(len=*), intent(in) :: name
      integer, intent(out) :: ncore,status
      character(len=8), intent(out) :: labels(:)
      call open_descriptor_source(source,descriptor_path(name),size(labels),status=status)
      if (status /= 0) return
      call install_descriptor_header(source,ncore,labels,status)
      if (status == 0) then
         deallocate(csf_counts_by_block)
         allocate(csf_counts_by_block(0:source%block_count),stat=status)
         if (status == 0) then
            csf_counts_by_block(0)=0
            csf_counts_by_block(1:)=source%block_lengths
         end if
      end if
      if (status /= 0) call close_descriptor_source(source)
   end subroutine open_csf_blocks

   subroutine read_csf_zero_order_counts(counts,status)
      integer, intent(out) :: counts(:), status
      call read_descriptor_zero_order_counts(source,counts,status)
   end subroutine read_csf_zero_order_counts

   subroutine load_csf_block(block_index,status)
      integer, intent(in) :: block_index
      integer, intent(out) :: status
      call read_descriptor_block(source,block_index,csf_subshell_occupations, &
                                 csf_subshell_quantum_codes,csf_coupling_quantum_codes,status)
      if (status /= 0 .or. block_index == source%block_count) call close_descriptor_source(source)
   end subroutine load_csf_block
end module csf_block_input
