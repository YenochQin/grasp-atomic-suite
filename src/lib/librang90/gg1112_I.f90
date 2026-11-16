      module gg1112_i
      interface
!
      subroutine gg1112(ik1,ik2,bk1,bk2,id1,id2,bd1,bd2,k1,  &
                        qm1,qm2,qm3,qm4,ww)
      use vast_kind_param, only:  double
      integer,      intent(in)               :: k1
      integer,      intent(in), dimension(7) :: ik1, ik2, id1, id2
      real(double), intent(in)               :: qm1, qm2, qm3, qm4
      real(double), intent(in), dimension(3) :: bk1, bk2, bd1, bd2
      real(double), intent(out)              :: ww
      end subroutine
      end interface
      end module
