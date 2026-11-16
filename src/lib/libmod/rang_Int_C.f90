      module rang_int_c
      implicit none (type, external)
      integer :: notfound, found_csf
      integer, pointer :: found(:)
      character(len=256), pointer :: c_shell(:), c_quant(:), c_coupl(:)
      integer :: num_in_blk(20)
      end module rang_int_c
