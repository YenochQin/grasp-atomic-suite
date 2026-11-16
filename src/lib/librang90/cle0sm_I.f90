      module cle0sm_i
      implicit none (type, external)
      interface
!
      subroutine cle0sm(q,qm,s,c,cm,a)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in)  :: q, qm, s, c, cm
      real(real64), intent(out) :: a
      end subroutine
      end interface
      end module
