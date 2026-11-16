      module recop1_i
      implicit none (type, external)
      interface
!
      subroutine recop1(ns,ja1,ka,ire,iat,recc)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in)       :: ns, ja1, ka, ire
      integer, intent(out)      :: iat
      real(real64), intent(out) :: recc
      end subroutine
      end interface
      end module
