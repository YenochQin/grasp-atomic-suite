      module brint4_i
      implicit none (type, external)
      interface
      subroutine brint4 (ia,ib,ic,id,nu,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ia, ib, ic, id, nu
      real(real64), intent(out) :: tegral
      end subroutine
      end interface
      end module
