      module ww1_i
      interface
!
      subroutine ww1(ik,bk,id,bd,k2,qm1,qm2,qm3,qm4,ww)
      use vast_kind_param, only:  double
      integer,      intent(in)               :: k2
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in)               :: qm1, qm2, qm3, qm4
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(out)              :: ww
      end subroutine
      end interface
      end module
