      module wj1_i
      implicit none (type, external)
      interface
!
      subroutine wj1(ik,bk,id,bd,k2,qm1,qm2,wj)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)               :: k2
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: qm1, qm2
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: wj
      end subroutine
      end interface
      end module
