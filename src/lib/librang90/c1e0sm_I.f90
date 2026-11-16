      module c1e0sm_i
      implicit none (type, external)
      interface
!
      subroutine c1e0sm(q,qm,c,cm,a)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in)  :: q, qm, c, cm
      real(real64), intent(out) :: a
      end subroutine
      end interface
      end module
