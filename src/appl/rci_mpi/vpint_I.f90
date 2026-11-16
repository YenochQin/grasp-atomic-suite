      module vpint_i
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine vpint (ia,ib,tegral)
      use vast_kind_param, only:  double
      integer , intent(inout) :: ia, ib
      real(double), intent(out) :: tegral
      end subroutine vpint
      end interface
      end module
