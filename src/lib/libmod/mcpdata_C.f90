      module mcpdata_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: ncoeff, nintg
      integer, pointer :: jann(:), jbnn(:), intgrl(:), intptr(:)
      real(real64), pointer :: cnn(:)
      end module mcpdata_c
