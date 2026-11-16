      module gg12_i
      interface
!
      subroutine gg12(ik1,ik2,bk1,bk2,id1,id2,bd1,bd2,qm1,qm2,ww)
      use vast_kind_param, only:  double
      integer,      intent(in), dimension(7) :: ik1, ik2, id1, id2
      real(double), intent(in)               :: qm1, qm2
      real(double), intent(in), dimension(3) :: bk1, bk2, bd1, bd2
      real(double), intent(out)              :: ww
      end subroutine
      end interface
      end module
