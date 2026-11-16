      module w1jjg_i
      implicit none (type, external)
      interface
!
      subroutine w1jjg(k1,qm1,qm2,ik,bk,id,bd,ww)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: qm1, qm2
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: ww
      end subroutine
      end interface
      end module
