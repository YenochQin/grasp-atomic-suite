      module suwjj_i
      implicit none (type, external)
      interface
!
      subroutine suwjj(k1,k2,ll,j1,j2,suw)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)  :: k1, k2, j1, j2
      integer,      intent(out) :: ll
      real(real64), intent(out) :: suw
      end subroutine
      end interface
      end module
