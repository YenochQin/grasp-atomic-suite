      module biotr1_i
      implicit none (type, external)
      interface
      subroutine biotr1 (pi, qi, nli, ninshli, pf, qf, nlf, ninshlf, ngrid, mxl&
         , scr, lscr, ntestg, cishl, cici, cfshl, cfci)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ngrid
      real(real64) :: pi(ngrid, *)
      real(real64) :: qi(ngrid, *)
      integer, intent(in) :: mxl
      integer, intent(in) :: nli(mxl)
      integer :: ninshli(mxl)
      real(real64) :: pf(ngrid, *)
      real(real64) :: qf(ngrid, *)
      integer, intent(in) :: nlf(mxl)
      integer :: ninshlf(mxl)
      real(real64), intent(inout) :: scr(*)
      integer, intent(in) :: lscr
      integer, intent(in) :: ntestg
      real(real64) :: cishl(*)
      real(real64) :: cici(*)
      real(real64) :: cfshl(*)
      real(real64) :: cfci(*)
      end subroutine
      end interface
      end module
