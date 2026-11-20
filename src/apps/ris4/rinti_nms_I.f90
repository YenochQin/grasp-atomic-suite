      module rinti_nms_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine rinti_nms (first_orbital, second_orbital, relativistic_nms_integral,leading_nms_integral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(out) :: relativistic_nms_integral, leading_nms_integral
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      end subroutine rinti_nms
      end interface
      end module rinti_nms_i
