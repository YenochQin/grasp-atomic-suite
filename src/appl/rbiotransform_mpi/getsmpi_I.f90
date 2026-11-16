      module gets_i
      implicit none (type, external)
      interface
      subroutine gets (s, nshli, nshlf)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nshli
      integer, intent(in) :: nshlf
      real(real64), intent(out) :: s(nshli, nshlf)
      end subroutine gets
      end interface
      end module gets_i
