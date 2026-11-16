      module cre_i
      implicit none (type, external)
      interface
      real(real64) function cre (kap1, k, kap2)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: kap1
      integer, intent(in) :: k
      integer, intent(in) :: kap2
!VAST...Calls: CLRX
      end function
      end interface
      end module
