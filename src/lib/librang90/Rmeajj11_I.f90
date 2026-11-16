      module rmeajj11_i
      interface
!
      subroutine rmeajj11(j1,j2,ll,s)
      use vast_kind_param, only:  double
      integer,      intent(in)  :: j1, j2, ll
      real(double), intent(out) :: s
      end subroutine
      end interface
      end module
