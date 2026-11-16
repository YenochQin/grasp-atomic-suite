      module brint2_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine brint2 (ia,ib,ic,id,k,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ia, ib, ic, id, k
      real(real64), intent(out) :: tegral
      end subroutine
      end interface
      end module
