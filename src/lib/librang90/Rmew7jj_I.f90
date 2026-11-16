      module rmew7jj_i
      interface
!
      subroutine rmew7jj(j1,j2,k1,k2,coef)
      use vast_kind_param, only:  double
      integer,      intent(in)  :: j1, j2, k1, k2
      real(double), intent(out) :: coef
      end subroutine
      end interface
      end module
