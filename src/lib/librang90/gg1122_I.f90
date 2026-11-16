      module gg1122_i
      interface
!
      subroutine gg1122(k1,k2,qm1,qm2,qm3,qm4,aa)
      use vast_kind_param, only:  double
      integer,      intent(in)               :: k1, k2
      real(double), intent(in)               :: qm1, qm2, qm3, qm4
      real(double), intent(out)              :: aa
      end subroutine
      end interface
      end module
