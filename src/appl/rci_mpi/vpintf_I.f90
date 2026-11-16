      module vpintf_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      real(real64)  function vpintf (ia,ib)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer , intent(in) :: ia, ib
      end function
      end interface
      end module
