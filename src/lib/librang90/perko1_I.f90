      module perko1_i
      implicit none (type, external)
      interface
!
      subroutine perko1(ja,bk,ik,bd,id)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer,      intent(in)                :: ja
      integer,      intent(out) :: ik(7), id(7)
      real(real64), intent(out) :: bk(3), bd(3)
      end subroutine
      end interface
      end module
