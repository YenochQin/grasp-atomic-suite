      module dampck_i
      implicit none (type, external)
      interface
      subroutine dampck (ipr, j, ed1, ed2)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: ipr
      integer, intent(in) :: j
      real(real64), intent(inout) :: ed1
      real(real64), intent(inout) :: ed2
      end subroutine
      end interface
      end module
