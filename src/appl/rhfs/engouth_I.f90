      module engouth_i
      implicit none (type, external)
      interface
      subroutine engouth (eav, e, jtot, ipar, ilev, nn, mode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: eav
      integer, intent(in) :: nn
      real(real64), intent(in) :: e(nn)
      integer, intent(in) :: jtot(nn)
      integer, intent(in) :: ipar(nn)
      integer, intent(in) :: ilev(nn)
      integer, intent(in) :: mode
      end subroutine
      end interface
      end module
