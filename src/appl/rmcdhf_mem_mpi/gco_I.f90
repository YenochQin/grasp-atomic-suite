      module gco_i
      implicit none (type, external)
      interface
      real(real64) function gco (k, ir, ia, ib)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer :: k
      integer :: ir
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      end function gco
      end interface
      end module gco_i
