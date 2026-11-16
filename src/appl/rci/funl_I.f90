      module funl_i
      implicit none (type, external)
      interface
      real(real64) function funl (x, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: x
      integer, intent(in) :: k
!...This routine performs I/O.
      end function
      end interface
      end module
