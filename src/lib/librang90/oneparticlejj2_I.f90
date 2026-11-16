      module oneparticlejj2_i
      implicit none (type, external)
      interface
!
      subroutine oneparticlejj2(ns,ka,ja,jb,coeff)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ns,ka,ja,jb
      real(real64), intent(out) :: coeff
      end subroutine
      end interface
      end module
