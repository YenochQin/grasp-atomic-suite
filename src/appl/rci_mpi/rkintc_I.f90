      module rkintc_i
      implicit none (type, external)
      interface
      subroutine rkintc (ia, ib, ic, id, k, tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!      INTEGER KMAX
!      PARAMETER (KMAX = 20)
      integer, intent(inout) :: ia
      integer, intent(inout) :: ib
      integer, intent(inout) :: ic
      integer, intent(inout) :: id
      integer, intent(inout) :: k
      real(real64), intent(out) :: tegral
      end subroutine rkintc
      end interface
      end module rkintc_i
