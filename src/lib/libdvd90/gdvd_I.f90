      module gdvd_i
      implicit none (type, external)
      abstract interface
      subroutine dvd_operator(n, m, b, c)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: n, m
      real(real64), intent(in) :: b(n,m)
      real(real64), intent(inout) :: c(n,m)
      end subroutine
      end interface
      interface
      subroutine gdvd (matrix_vector_callback, n, lim, diag, ilow, ihigh, iselec, niv, mblock&
         , crite, critc, critr, ortho, maxiter, work, iwrsz, iwork, iiwsz&
         , hiend, nloops, nmv, ierr)
      use, intrinsic :: iso_fortran_env, only: real64
      import dvd_operator
      implicit none (type, external)
      procedure(dvd_operator) :: matrix_vector_callback
      integer, intent(in) :: n
      integer, intent(in) :: lim
      real(real64), intent(inout) :: diag(n)
      integer :: ilow
      integer :: ihigh
      integer :: iselec(lim)
      integer, intent(inout) :: niv
      integer :: mblock
      real(real64) :: crite
      real(real64) :: critc
      real(real64) :: critr
      real(real64) :: ortho
!VAST...Dummy argument ORTHO is not referenced in this routine.
      integer :: maxiter
      integer, intent(in) :: iwrsz
      real(real64), intent(inout) :: work(iwrsz)
      integer, intent(in) :: iiwsz
      integer, intent(inout) :: iwork(iiwsz)
      logical :: hiend
      integer :: nloops
      integer, intent(out) :: nmv
      integer, intent(out) :: ierr
!VAST...Calls: OP, DVDSON
!...This routine performs I/O.
      end subroutine
      end interface
      end module
