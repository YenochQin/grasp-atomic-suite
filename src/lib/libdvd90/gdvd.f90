
!***************************************************************************

      subroutine gdvd(matrix_vector_callback, n, lim, diag, ilow, ihigh, iselec, niv, mblock, crite&
         , critc, critr, ortho, maxiter, work, iwrsz, iwork, iiwsz, hiend, &
         nloops, nmv, ierr)
!      Written by M. Saparov
!
! Note:
! HIEND, ISELEC() not used outside dvdson
!***************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use gdvd_i, only: dvd_operator
      use dvdson_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      procedure(dvd_operator) :: matrix_vector_callback
      integer, intent(in) :: n
      integer, intent(in) :: lim
      integer  :: ilow
      integer  :: ihigh
      integer, intent(inout) :: niv
      integer  :: mblock
      integer  :: maxiter
      integer, intent(in) :: iwrsz
      integer, intent(in) :: iiwsz
      integer  :: nloops
      integer, intent(out) :: nmv
      integer  :: ierr
      real(real64)  :: crite
      real(real64)  :: critc
      real(real64)  :: critr
      real(real64)  :: ortho
      logical  :: hiend
      integer  :: iselec(lim)
      integer, intent(inout) :: iwork(iiwsz)
      real(real64), intent(inout) :: diag(n)
      real(real64), intent(inout) :: work(iwrsz)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: irev(7)
      integer :: noc, irc, nb, iw1, iw2, iw3, iiw, iin, iw4, icur, i, indx, j
      real(real64) :: value, ovalue, rnorm, epsil
!
!**********************************************************************
!  NOC = number of orthogonalization constraints
!  IRC = reverse communication switch

      noc = 0
      irc = 0

!***********************************************************************
! CALLING DAVIDSON with reverse communication
!***********************************************************************
! Initial estimates
!
      if (niv == 0) then
         write (6, *) 'GDVD Error : No initial estimate!!!'
         ierr = -1000
         return
      endif

      !ttt=etime_(tarray)
      nmv = 0
      !xhh print*, 'MBLOCK = ', mblock
      !xhh print *, ' gdvd:  niv = ', niv
   99 continue
      call dvdson (irc, irev, n, lim, noc, ilow, ihigh, iselec, niv, mblock, &
         crite, critc, critr, maxiter, work, iwrsz, iwork, iiwsz, hiend, nloops&
         , ierr)
!
! * * Start Reverse Communication * * * * * * * * * * * * * * * * * * *
!
      nb = irev(1)
      iw1 = irev(2)
      iw2 = irev(3)
      iw3 = irev(4)
      iiw = irev(5)
      iin = irev(6)
      iw4 = irev(7)

      if (irc == 1) then
!********* ..Preconditioning. Solve NB times(M work(iw2)=work(iw1))
!          ..Results always on work(iw2)

         icur = 0
         do i = 1, nb
            indx = iwork(iiw+i-1) - 1
            value = work(iw3+indx)
            ovalue = work(iw4+indx)
            rnorm = work(iin+indx)
            epsil = work(iin+lim+indx)
! The current approximation of the eigenvector is x=Bc
! If needed it should be saved e.g: call dcopy(N,work(iw2),1,curx,1)

!            write(*,11) nloops,value,rnorm
   11       format('It ',i4,'  Dl',d10.3,' Res:',d10.3)
!
!-------------
! (M-lI)      .. Compute temporarily (M-valueI) (M preconditioner)
!             .. Needed for (M-valueI)^-1 res
!             .. Here M is the DIAG
            diag(:n) = diag(:n) - value
!-------------
!* Choice of Diagonal preconditioning
            do j = 1, n
               if (abs(diag(j)) > 1.0e-6_real64) then
                  work(iw2+icur+j-1) = work(iw1+icur+j-1)/diag(j)
               else
                  work(iw2+icur+j-1) = work(iw1+icur+j-1)*1.0e6_real64
               endif
            end do
!*-------------
!* e.g: For No preconditioner: Lanczos
!*            call dcopy(N,WORK(iw1+icur),1,WORK(iw2+icur),1)
!*-------------
!* (M+lI)      .. Restore (M+valueI)
            diag(:n) = diag(:n) + value

            icur = icur + n
         end do

         go to 99

!*********
      else if (irc==2 .or. irc==3) then
!********* ..Matrix-vector multiply.
         call matrix_vector_callback (n, nb, work(iw1), work(iw2))
         nmv = nmv + nb

         go to 99
      endif
! * * * * End of Reverse Communication * * * * * * * * * * * * * * * * *

      return
      end subroutine gdvd
