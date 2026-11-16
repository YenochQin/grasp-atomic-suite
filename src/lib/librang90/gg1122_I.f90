      module gg1122_i
      implicit none (type, external)
      interface
!
      subroutine gg1122(k1,k2,qm1,qm2,qm3,qm4,aa)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)               :: k1, k2
      real(real64), intent(in)               :: qm1, qm2, qm3, qm4
      real(real64), intent(out)              :: aa
      end subroutine
      end interface
      end module
