      module gg1233_i
      interface
!
      subroutine gg1233(ik1,ik2,ik3,bk1,bk2,bk3,id1,id2,id3,bd1,  &
      bd2,bd3,k1,qm1,qm2,qm3,qm4,ww)
      use vast_kind_param, only:  double
      integer,      intent(in)               :: k1
      integer,      intent(in), dimension(7) :: ik1,ik2,ik3,id1,id2,id3
      real(double), intent(in)               :: qm1,qm2,qm3,qm4
      real(double), intent(in), dimension(3) :: bk1,bk2,bk3,bd1,bd2,bd3
      real(double), intent(out)              :: ww
      end subroutine
      end interface
      end module
