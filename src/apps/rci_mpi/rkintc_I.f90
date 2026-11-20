      module rkintc_i
      implicit none (type, external)
      interface
      subroutine rkintc (orbital_a, orbital_b, orbital_c, orbital_d, multipole_rank, tegral)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!      INTEGER KMAX
!      PARAMETER (KMAX = 20)
      integer, intent(inout) :: orbital_a
      integer, intent(inout) :: orbital_b
      integer, intent(inout) :: orbital_c
      integer, intent(inout) :: orbital_d
      integer, intent(inout) :: multipole_rank
      real(real64), intent(out) :: tegral
      end subroutine rkintc
      end interface
      end module rkintc_i
