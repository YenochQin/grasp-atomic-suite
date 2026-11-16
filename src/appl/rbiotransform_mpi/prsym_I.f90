      module prsym_i
      implicit none (type, external)
      interface
      subroutine prsym (a, matdim)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: a(1)
      integer, intent(in) :: matdim
!...This routine performs I/O.
      end subroutine
      end interface
      end module
