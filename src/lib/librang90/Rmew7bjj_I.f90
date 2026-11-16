      module rmew7bjj_i
      implicit none (type, external)
      interface
!
      subroutine rmew7bjj(j1,j2,k1,k2,coef)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)  :: j1, j2, k1, k2
      real(real64), intent(out) :: coef
      end subroutine
      end interface
      end module
