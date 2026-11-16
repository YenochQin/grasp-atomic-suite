      module rinti_nms_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine rinti_nms (j, k, valnmsk123,valnmsk1)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(out) :: valnmsk123, valnmsk1
      integer, intent(in) :: j
      integer, intent(in) :: k
      end subroutine
      end interface
      end module
