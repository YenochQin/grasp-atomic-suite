      module bndinv_i
      implicit none (type, external)
      interface
      subroutine bndinv (a, el, n, determ, epsil, itest, nsize)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nsize
      real(real64), intent(inout) :: a(nsize, nsize)
      real(real64), intent(inout) :: el(nsize, nsize)
      integer, intent(in) :: n
      real(real64), intent(out) :: determ
      real(real64), intent(in) :: epsil
      integer, intent(out) :: itest
      end subroutine bndinv
      end interface
      end module bndinv_i
