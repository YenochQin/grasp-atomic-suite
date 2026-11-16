      module engoutgg_i
      implicit none (type, external)
      interface
      subroutine engoutgg(e,ilev,nn,mode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nn
      real(real64), intent(in) :: e(nn)
      integer, intent(in) :: ilev(nn)
      integer, intent(in) :: mode
      end subroutine
      end interface
      end module
