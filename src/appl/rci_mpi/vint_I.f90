      module vint_i
      interface
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine vint (ia,ib,tegral)
      use vast_kind_param, only:  double
      integer , intent(in) :: ia, ib
      real(double), intent(out) :: tegral
      end subroutine vint
      end interface
      end module
