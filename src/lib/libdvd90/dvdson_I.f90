      module dvdson_i
      implicit none (type, external)
      interface
      subroutine dvdson (irc, irev, n, lim, noc, ilow, ihigh, iselec, niv&
         , mblock, crite, critc, critr, maxiter, work, iwrsz, iwork, iiwsz&
         , hiend, nloops, ierr)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: irc
      integer, intent(inout) :: irev(*)
      integer, intent(in) :: n
      integer, intent(in) :: lim
      integer, intent(in) :: noc
      integer, intent(inout) :: ilow
      integer, intent(inout) :: ihigh
      integer, intent(inout) :: iselec(lim)
      integer, intent(inout) :: niv
      integer, intent(in) :: mblock
      real(real64) :: crite
      real(real64) :: critc
      real(real64) :: critr
      integer :: maxiter
      integer, intent(in) :: iwrsz
      real(real64), intent(inout) :: work(iwrsz)
      integer, intent(in) :: iiwsz
      integer, intent(inout) :: iwork(iiwsz)
      logical, intent(inout) :: hiend
      integer, intent(inout) :: nloops
      integer, intent(inout) :: ierr
!VAST...Calls: INITDVD, DVDRVR, DSCAL, DCOPY
!...This routine performs I/O.
      end subroutine
      end interface
      end module
