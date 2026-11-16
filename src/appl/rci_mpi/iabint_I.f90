      module iabint_i
      implicit none (type, external)
      interface
      subroutine iabint (ia,ib,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: ia, ib
     real(real64), intent(out) :: tegral
      end subroutine
      end interface
      end module
