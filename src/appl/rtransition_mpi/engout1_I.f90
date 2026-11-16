      module engout1_i
      implicit none (type, external)
      interface
      subroutine engout1 (eav, e, jtot, ipar, ilev, nn, mode, k)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: eav
      integer, intent(in) :: nn
      real(real64), intent(in) :: e(nn)
      integer, intent(in) :: jtot(nn)
      integer, intent(in) :: ipar(nn)
      integer, intent(in) :: ilev(nn)
      integer :: mode
!VAST...Dummy argument MODE is not referenced in this routine.
      integer, intent(in) :: k
      end subroutine engout1
      end interface
      end module engout1_i
