      module count_i
      implicit none (type, external)
      interface
      subroutine count (fr, mtpfr, nncff, sgn)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: fr(max_radial_points)
      integer, intent(in) :: mtpfr
      integer, intent(out) :: nncff
      real(real64), intent(out) :: sgn
!VAST.../COUN/ THRESH(IN)
!VAST.../DEF4/ ACCY(IN)
      end subroutine
      end interface
      end module
