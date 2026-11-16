      module iabint_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine iabint (ia,ib,tegral)
      use vast_kind_param, only:  double
      integer, intent(inout) :: ia, ib
     real(double), intent(out) :: tegral
      end subroutine
      end interface
      end module
