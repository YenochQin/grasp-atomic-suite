      module gg1234_i
      implicit none (type, external)
      interface
!
      subroutine gg1234(ik1,ik2,ik3,ik4,bk1,bk2,bk3,bk4,id1,id2,  &
      id3,id4,bd1,bd2,bd3,bd4,qm1,qm2,qm3,qm4,ww)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in) :: ik1(7),ik2(7),ik3(7),ik4(7),id1(7),id2(7),id3(7),id4(7)
      real(real64), intent(in)               :: qm1,qm2,qm3,qm4
      real(real64), intent(in) :: bk1(3),bk2(3),bk3(3),bk4(3),bd1(3),bd2(3),bd3(3),bd4(3)
      real(real64), intent(out)              :: ww
      end subroutine
      end interface
      end module
