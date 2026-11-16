      module rang_int_c
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: notfound, found_csf
      integer, dimension(:), pointer :: found
      character(len=256), dimension(:), pointer :: c_shell, c_quant, c_coupl
      integer, dimension(20) :: num_in_blk
      end module rang_int_c
