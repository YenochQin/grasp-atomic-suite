      module matelt_i
      implicit none (type, external)
      interface
      subroutine matelt (i1, k, i2, apart, gjpart, dgjpart)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: i1
      integer, intent(in) :: k
      integer, intent(in) :: i2
      real(real64), intent(out) :: apart
      real(real64), intent(out) :: gjpart
      real(real64), intent(out) :: dgjpart
      end subroutine
      end interface
      end module
