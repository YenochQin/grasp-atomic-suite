      module initdvd_i
      implicit none (type, external)
      interface
      subroutine initdvd (irc, irev, n, noc, niv, nume, lim, hiend, scra1&
         , orthobasis, basis, ab, s)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: irc
      integer, intent(out) :: irev(*)
      integer, intent(in) :: n
      integer, intent(in) :: noc
      integer, intent(inout) :: niv
      integer, intent(in) :: nume
      integer, intent(in) :: lim
      logical, intent(in) :: hiend
      real(real64) :: scra1(*)
      real(real64) :: orthobasis(n*(noc + lim))
      real(real64) :: basis(n*lim)
      real(real64) :: ab(n*lim)
      real(real64) :: s(*)
!VAST...Calls: DINIT, DCOPY, MGS_NRM, DSCAL, ADDS
      end subroutine
      end interface
      end module
