      module c1e1sm_i
      implicit none (type, external)
      interface
!
      subroutine c1e1sm(q,qm,sm,c,cm,a)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in)  :: q, qm, sm, c, cm
      real(real64), intent(out) :: a
      end subroutine
      end interface
      end module
