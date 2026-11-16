      module dmerge_dnicmv_i
      implicit none (type, external)
      interface
      subroutine dmerge_dnicmv (n, db, dc, da, dconst, dl)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n
      real(real64), intent(in) :: db(*)
      real(real64), intent(inout) :: dc(*)
      real(real64), intent(in) :: da(n)
      real(real64), intent(in) :: dconst
      real(real64), intent(out) :: dl
      end subroutine
      end interface
      end module
