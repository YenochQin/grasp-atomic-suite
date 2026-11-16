      module gg12_i
      implicit none (type, external)
      interface
!
      subroutine gg12(ik1,ik2,bk1,bk2,id1,id2,bd1,bd2,qm1,qm2,ww)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in) :: ik1(7), ik2(7), id1(7), id2(7)
      real(real64), intent(in)               :: qm1, qm2
      real(real64), intent(in) :: bk1(3), bk2(3), bd1(3), bd2(3)
      real(real64), intent(out)              :: ww
      end subroutine
      end interface
      end module
