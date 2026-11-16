      module engout_i
      implicit none (type, external)
      interface
      subroutine engout (eav, e, jtot, ipar, ilev, nn, mode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: eav
      integer, intent(in) :: nn
      real(real64), intent(in) :: e(nn)
      integer, intent(in) :: jtot
      integer, intent(in) :: ipar
      integer, intent(in) :: ilev(nn)
      integer, intent(in) :: mode
      end subroutine engout
      end interface
      end module engout_i
