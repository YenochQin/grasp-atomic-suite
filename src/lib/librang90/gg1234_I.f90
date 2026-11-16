      module gg1234_i
      interface
!
      subroutine gg1234(ik1,ik2,ik3,ik4,bk1,bk2,bk3,bk4,id1,id2,  &
      id3,id4,bd1,bd2,bd3,bd4,qm1,qm2,qm3,qm4,ww)
      use vast_kind_param, only:  double
      integer,      intent(in), dimension(7) :: ik1,ik2,ik3,ik4,id1,id2,id3,id4
      real(double), intent(in)               :: qm1,qm2,qm3,qm4
      real(double), intent(in), dimension(3) :: bk1,bk2,bk3,bk4,bd1,bd2,bd3,bd4
      real(double), intent(out)              :: ww
      end subroutine
      end interface
      end module
