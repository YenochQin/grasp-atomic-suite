      module mcpdata_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ncoeff, nintg
      integer, dimension(:), pointer :: jann, jbnn, intgrl, intptr
      real(double), dimension(:), pointer :: cnn
      end module mcpdata_c
