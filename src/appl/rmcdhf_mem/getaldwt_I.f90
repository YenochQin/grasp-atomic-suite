      module getaldwt_i
      implicit none (type, external)
      interface
      subroutine getaldwt (ncf, wt)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ncf
      real(real64), intent(inout) :: wt(ncf)
      end subroutine getaldwt
      end interface
      end module getaldwt_i
