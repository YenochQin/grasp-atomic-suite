      module rwjj_i
      interface
!
      subroutine rwjj(j,j1,j2,k1,k2,coef)
      use vast_kind_param, only:  double
      integer,      intent(in)  :: j, j1, j2, k1, k2
      real(double), intent(out) :: coef
      end subroutine
      end interface
      end module
