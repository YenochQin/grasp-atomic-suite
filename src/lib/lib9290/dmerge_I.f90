      module dmerge_i
      implicit none (type, external)
      interface
      subroutine dmerge (n, db, dc, idy, da, dconst, dl)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      real(real64), intent(in) :: db(*)
      real(real64), intent(inout) :: dc(*)
      integer, intent(in) :: idy(n)
      real(real64), intent(in) :: da(n)
      real(real64), intent(in) :: dconst
      real(real64), intent(out) :: dl
      end subroutine
      end interface
      end module
