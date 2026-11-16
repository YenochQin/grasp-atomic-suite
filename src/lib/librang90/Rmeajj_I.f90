      module rmeajj_i
      interface
!
      subroutine rmeajj(ll,it,lq,j,its,lqs,j1s,coef)
      use vast_kind_param, only:  double
      integer,      intent(in)  :: ll, it, lq, j, its, lqs, j1s
      real(double), intent(out) :: coef
      end subroutine
      end interface
      end module
