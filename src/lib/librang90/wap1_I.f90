      module wap1_i
      interface
!
      subroutine wap1(ik,bk,id,bd,k1,bk2,qm1,qm2,qm3,wa)
      use vast_kind_param, only:  double
      integer,      intent(in)               :: k1
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in)               :: bk2, qm1, qm2, qm3
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(out)              :: wa
      end subroutine
      end interface
      end module
