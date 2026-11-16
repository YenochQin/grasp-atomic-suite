      module vint_i
      implicit none (type, external)
      interface
      subroutine vint (ia,ib,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ia, ib
      real(real64), intent(out) :: tegral
      end subroutine vint
      end interface
      end module
