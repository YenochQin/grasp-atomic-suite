      module ww1_i
      implicit none (type, external)
      interface
!
      subroutine ww1(ik,bk,id,bd,k2,qm1,qm2,qm3,qm4,ww)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)               :: k2
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: qm1, qm2, qm3, qm4
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: ww
      end subroutine
      end interface
      end module
