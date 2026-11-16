      module adds_i
      implicit none (type, external)
      interface
      subroutine adds (n, lim, kpass, nncv, basis, ab, s)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(in) :: lim
      integer, intent(in) :: kpass
      integer, intent(in) :: nncv
      real(real64) :: basis(n*lim)
      real(real64) :: ab(n*lim)
      real(real64), intent(inout) :: s(lim*(lim + 1)/2)
!VAST...Calls: DDOT
      end subroutine
      end interface
      end module
