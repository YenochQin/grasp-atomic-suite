subroutine setcsl(name,ncore,nblkin,idblk)
   use rci_csf_block_reader, only: open_csf_blocks
   implicit none (type, external)
   character(len=*), intent(in) :: name
   integer, intent(out) :: ncore
   integer, intent(in) :: nblkin
   character(len=8), intent(out) :: idblk(*)
   integer :: status
   call open_csf_blocks(name,ncore,idblk(:nblkin),status)
   if (status /= 0) error stop "RCI descriptor header failed"
end subroutine setcsl
