      module brint3_i
      implicit none (type, external)
      interface
      subroutine brint3 (ia,ib,ic,id,nu,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ia, ib, ic, id, nu
      real(real64), intent(out) :: tegral
      end subroutine brint3
      end interface
      end module brint3_i
