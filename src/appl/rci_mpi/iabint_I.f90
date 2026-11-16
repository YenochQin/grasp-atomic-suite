      module iabint_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine iabint (ia,ib,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: ia, ib
     real(real64), intent(out) :: tegral
      end subroutine
      end interface
      end module
