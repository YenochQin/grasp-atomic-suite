      module funk_i
      implicit none (type, external)
      interface
      real(real64) function funk (x, n)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: x
      integer, intent(in) :: n
!...This routine performs I/O.
      end function funk
      end interface
      end module funk_i
