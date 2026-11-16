      module awp1_i
      implicit none (type, external)
      interface
!
      subroutine awp1( ik, bk, id, bd, k1, bk2, qm1, qm2, qm3, aw)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: bk2, qm1, qm2, qm3
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: aw
      end subroutine
      end interface
      end module
