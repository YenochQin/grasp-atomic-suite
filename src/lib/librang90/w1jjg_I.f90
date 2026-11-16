      module w1jjg_i
      interface
!
      subroutine w1jjg(k1,qm1,qm2,ik,bk,id,bd,ww)
      use vast_kind_param, only:  double
      integer,      intent(in)               :: k1
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in)               :: qm1, qm2
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(out)              :: ww
      end subroutine
      end interface
      end module
