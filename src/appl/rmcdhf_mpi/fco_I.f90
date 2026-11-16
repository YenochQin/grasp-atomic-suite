      module fco_i
      implicit none (type, external)
      interface
      real(real64) function fco (k, ir, ia, ib)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: k
      integer :: ir
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      end function
      end interface
      end module
