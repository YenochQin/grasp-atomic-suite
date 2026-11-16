      module mgs_nrm_i
      implicit none (type, external)
      interface
      subroutine mgs_nrm (n, kp, new, scra, b)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: kp
      integer, intent(in) :: new
      real(real64), intent(inout) :: scra(new)
      real(real64), intent(inout) :: b((kp + new)*n)
!VAST...Calls: DGEMV, DDOT, DSCAL
      end subroutine
      end interface
      end module
