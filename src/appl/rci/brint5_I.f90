      module brint5_i
      implicit none (type, external)
      interface
      subroutine brint5 (ia,ib,ic,id,nu,tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ia, ib, ic, id, nu
      real(real64), intent(out) :: tegral
      end subroutine brint5
      end interface
      end module brint5_i
