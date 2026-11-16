      module suwjj_i
      interface
!
      subroutine suwjj(k1,k2,ll,j1,j2,suw)
      use vast_kind_param, only:  double
      integer,      intent(in)  :: k1, k2, j1, j2
      integer,      intent(out) :: ll
      real(double), intent(out) :: suw
      end subroutine
      end interface
      end module
