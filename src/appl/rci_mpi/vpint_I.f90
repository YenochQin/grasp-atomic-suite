      module vpint_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine vpint (ia,ib,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer , intent(inout) :: ia, ib
      real(real64), intent(out) :: tegral
      end subroutine vpint
      end interface
      end module
