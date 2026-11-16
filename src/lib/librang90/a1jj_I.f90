      module a1jj_i
      implicit none (type, external)
      interface
!
      subroutine a1jj (ik, bk, id, bd, qm1, a)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(in)               :: qm1
      real(real64), intent(out)              :: a
      end subroutine
      end interface
      end module
