      module rmeajj_i
      implicit none (type, external)
      interface
!
      subroutine rmeajj(ll,it,lq,j,its,lqs,j1s,coef)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)  :: ll, it, lq, j, its, lqs, j1s
      real(real64), intent(out) :: coef
      end subroutine
      end interface
      end module
