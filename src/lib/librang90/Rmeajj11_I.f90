      module rmeajj11_i
      implicit none (type, external)
      interface
!
      subroutine rmeajj11(j1,j2,ll,s)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)  :: j1, j2, ll
      real(real64), intent(out) :: s
      end subroutine
      end interface
      end module
