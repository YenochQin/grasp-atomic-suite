      module cxk_i
      implicit none (type, external)
      interface
      subroutine cxk (s, is, kaps, nu, k, ibr, iex)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(inout) :: s(12)
      integer, intent(in) :: is(4)
      integer, intent(in) :: kaps(4)
      integer, intent(in) :: nu
      integer, intent(in) :: k
      integer, intent(in) :: ibr
      integer, intent(in) :: iex
      end subroutine cxk
      end interface
      end module cxk_i
