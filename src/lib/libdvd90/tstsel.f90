module tstsel_i
    implicit none (type, external)
    private
    public :: tstsel
contains
      logical function tstsel (kpass, nume, neig, iselec, svec, eigval, icv, &
         crite, critc, rowlast, ind, oldval, nncv, incv)
!=======================================================================
!
!     Called by: DVDRVR
!
!     It first checks if the wanted eigenvalues have reached
!     convergence and updates OLDVAL. Second, for each wanted and non
!     converged eigenvector, it finds the largest absolute coefficient
!       of the NNCV last added vectors (from SVEC) and if not coverged,
!     places it in ROWLAST. IND has the corresponding indices.
!     Third, it sorts ROWLAST in decreasing order and places the
!     corresponding indices in the array INCV. The index array INCV
!       and the number of unconverged pairs NNCV, are passed to DVDRVR.
!      Later in NEWVEC only the first MBLOCK of NNCV pairs
!      will be targeted, since if ROWLAST(i)>ROWLAST(j)
!      then approximately RESIDUAL(i)>RESIDUAL(j)
!
!     Subroutines called
!     IDAMAX
!-----------------------------------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!     USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE idamax_I
      use suite_linear_algebra, only: idamax
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: kpass
      integer , intent(in) :: nume
      integer , intent(in) :: neig
      integer , intent(inout) :: nncv
      real(real64) , intent(in) :: crite
      real(real64) , intent(in) :: critc
      integer , intent(in) :: iselec(neig)
      integer, intent(inout) :: icv(nume)
      integer , intent(inout) :: ind(neig)
      integer , intent(out) :: incv(neig)
      real(real64) , intent(in) :: svec(kpass*nume)
      real(real64) , intent(in) :: eigval(nume)
      real(real64), intent(inout) :: rowlast(neig)
      real(real64) , intent(in) :: oldval(nume)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nnce, i, ival, icnt, icur, l, indx, itemp
      real(real64) :: tmax, temp
      logical :: done
!-----------------------------------------------
!-----------------------------------------------------------------------
!
!   on entry
!   -------
!   KPASS       current dimension of the expanding Basis
!   NUME        Largest index of the wanted eigenvalues.
!   NEIG        number of wanted eigenvalues of original matrix
!   ISELEC index array of the wanted eigenvalues.
!   SVEC        the eigenvectors of the small system
!   EIGVAL The NUME lowest eigenvalues of the small problem
!   ICV  Index of converged eigenpairs.ICV(i)=1 iff eigenpair i
!      has converged, and ICV(i)=0 if eigenpair i has not.
!   CRITE,CRITC Convergence thresholds for eigenvalues and coefficients
!   ROWLAST    scratch array, keeping the largest absolute coefficient
!      of the NNCV last rows of Svec.
!   IND  scratch array, temporary keeping the indices of Rowlast
!   OLDVAL The previous iteration's eigenvalues.
!
!   on exit
!   -------
!   NNCV   Number of non converged eigenvectors (to be targeted)
!   INCV  Index to these columns in decreasing order of magnitude
!   TSTSEL     true if convergence has been reached
!
!-----------------------------------------------------------------------
!       !include 'mpif.h'

      done = .false.
!
! Test all wanted eigenvalues for convergence under CRITE
!
      nnce = 0
      do i = 1, neig
         ival = iselec(i)
         if (abs(oldval(ival)-eigval(ival)) < crite) cycle
         nnce = nnce + 1
      end do
      if (nnce == 0) then
!        IF (MYID == 0) WRITE (6, *) 'converged by crite'
         tstsel = .true.
         return
      endif
!
! Find the maximum element of the last NNCV coefficients of unconverged
! eigenvectors. For those unconverged coefficients, put their indices
! to IND and find their number NNCV
!
      icnt = 0
      do i = 1, neig
!Rasa      IF (ICV(ISELEC(I)).EQ.0) THEN
!             ..Find coefficient and test for convergence
         icur = kpass*iselec(i)
         tmax = abs(svec(icur))
         do l = 1, nncv - 1
            tmax = max(tmax,abs(svec(icur-l)))
         end do
         if (tmax < critc) then
!                ..this  coefficient converged
            icv(iselec(i)) = 1
         else
!                ..Not converged. Add it to the list.
!Rasa -- start change
            icv(iselec(i)) = 0
!Rasa -- end change
            icnt = icnt + 1
            ind(icnt) = iselec(i)
            rowlast(icnt) = tmax
         endif
!Rasa      ENDIF
      end do

      nncv = icnt
      if (nncv == 0) then
         done = .true.
      !IF (NNCV.EQ.0) print*, 'converged by critc', kpass
!        WRITE (6, *) 'converged by critc', KPASS, MYID
      endif
!
! Sort the ROWLAST elements interchanging their indices as well
!
      do i = 1, nncv
         indx = idamax(nncv - i + 1,rowlast(i),1)
         incv(i) = ind(indx+i-1)

         temp = rowlast(indx+i-1)
         rowlast(indx+i-1) = rowlast(i)
         rowlast(i) = temp
         itemp = ind(indx+i-1)
         ind(indx+i-1) = ind(i)
         ind(i) = itemp
      end do

      tstsel = done
      return
      end function tstsel
end module tstsel_i
