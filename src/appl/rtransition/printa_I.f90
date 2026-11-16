      module printa_i
      implicit none (type, external)
      interface
      subroutine printa (asfa, asfb, i, j, omega, factor, lines, lsame)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      logical, intent(in) :: lsame
      real(real64), intent(in) :: asfa
      real(real64), intent(in) :: asfb
      integer :: i
      integer, intent(in) :: j
      real(real64), intent(in) :: omega
      real(real64), intent(in) :: factor
      integer, intent(inout) :: lines
      end subroutine
      end interface
      end module
