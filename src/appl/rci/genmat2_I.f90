      module genmat2_i
      implicit none (type, external)
      interface
      subroutine genmat2 (irestart, nelmnt_a, elsto)
      use, intrinsic :: iso_fortran_env, only: real64, int64
      implicit none (type, external)
      integer, intent(in)        :: irestart
      integer(int64), intent(out) :: nelmnt_a
      real(real64), intent(in)   :: elsto
      end subroutine genmat2
      end interface
      end module genmat2_i
