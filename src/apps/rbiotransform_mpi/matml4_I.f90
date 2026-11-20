      module matml4_i
      implicit none (type, external)
      interface
      subroutine matml4 (c, a, b, ncrow, nccol, narow, nacol, nbrow, nbcol&
         , itrnsp)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ncrow
      integer, intent(in) :: nccol
      real(real64), intent(inout) :: c(ncrow, nccol)
      integer, intent(in) :: narow
      integer, intent(in) :: nacol
      real(real64), intent(in) :: a(narow, nacol)
      integer, intent(in) :: nbrow
      integer, intent(in) :: nbcol
      real(real64), intent(in) :: b(nbrow, nbcol)
      integer, intent(in) :: itrnsp
      end subroutine matml4
      end interface
      end module matml4_i
