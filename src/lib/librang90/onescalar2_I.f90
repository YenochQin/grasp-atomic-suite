      module onescalar2_i
      interface
!
      subroutine onescalar2(jja,jjb,ja,jb,coeff)
      use vast_kind_param, only:  double
      integer, intent(in) :: jja,jjb,ja,jb
      real(double), intent(out) :: coeff
      end subroutine
      end interface
      end module
