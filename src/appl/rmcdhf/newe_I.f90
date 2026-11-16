      module newe_i
      implicit none (type, external)
      interface
      subroutine newe (j, sgn, nprime, mx, deleps, fail, inv)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      real(real64), intent(in) :: sgn
      integer, intent(in) :: nprime
      integer, intent(in) :: mx
      real(real64), intent(inout) :: deleps
      logical, intent(out) :: fail
      integer, intent(out) :: inv
      end subroutine
      end interface
      end module
