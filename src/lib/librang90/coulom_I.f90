      module coulom_i
      interface
!
      subroutine coulom(j1,j2,j3,j4,l1,l2,l3,l4,k,aa)
      use vast_kind_param, only:  double
      integer, intent(in) :: j1, j2, j3, j4, l1, l2, l3, l4, k
      real(double),intent(out)  :: aa
      end subroutine
      end interface
      end module
