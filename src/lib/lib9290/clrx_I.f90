      module clrx_i
      implicit none (type, external)
      interface
      real(real64) function clrx (kappaa, k, kappab)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, parameter :: mfact = 500
      integer, intent(in) :: kappaa
      integer, intent(in) :: k
      integer, intent(in) :: kappab
!VAST.../FACTS/ GAM(IN)
      end function
      end interface
      end module
