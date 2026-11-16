      module auxblk_i
      implicit none (type, external)
      interface
      subroutine auxblk (j2max, atwinv)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j2max
      real(real64), intent(out) :: atwinv
      end subroutine
      end interface
      end module
