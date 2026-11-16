module dvdrvr_i
    implicit none (type, external)
    private
    public :: dvdrvr
contains
      subroutine dvdrvr(irc, irev, n, hiend, lim, mblock, noc, nume, niv, neig&
         , iselec, crite, critc, critr, maxiter, eigval, basis, orthobasis, ab&
         , s, temps, svec, scra1, iscra2, incv, icv, oldval, nloops, ierr)
!=======================================================================
!     called by DVDSON
!
!     Driver routine implementing Davidson's main loop. On entry it
!     is given the Basis, the work matrix D=AB and the small symmetric
!     matrix to be solved, S=B^TAB (as found by ADDS). In each step
!     the small problem is solved by calling DSPEVX.
!     TSTSEL tests for eigenvalue convergence and selects the next
!     pairs to be considered for targeting (as a block).
!       NEWVEC computes the new vectors (block) to be added in the
!     expanding basis, and tests for residual convergence.
!     After Preconditioning ((M-l_i)d=res), the basis is orthogonalized.
!       Then the matrix-vector multiplication finds the new vectors of D
!       (Dnew=ABnew), and the new small problem S, is calculated.
!       The algorithm is repeated.
!     In case of a large expanding basis (KPASS=LIM) the Basis, AB,
!     SVEC and S are collapsed.
!     At the end the current eigenvector estimates are computed as
!     well as the residuals and eigenvalue differences.
!
!     Subroutines called:
!     DSPEVX, MULTBC, TSTSEL, OVFLOW, NEWVEC, ADDS,
!     DCOPY, DDOT, DAXPY
!-----------------------------------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!     USE MPI_C,           ONLY: MYID
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use suite_linear_algebra, only: slamch
      use ovflow_i
      use tstsel_i
      use newvec_i
      use mgs_nrm_i
      use adds_i
      use multbc_i
      use suite_linear_algebra, only: daxpy, dcopy, ddot, dscal, dspevx
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: irc
      integer, intent(in) :: n
      integer, intent(in) :: lim
      integer  :: mblock
      integer, intent(in) :: noc
      integer, intent(in) :: nume
      integer, intent(in) :: niv
      integer, intent(in) :: neig
      integer, intent(in) :: maxiter
      integer, intent(inout) :: nloops
      integer, intent(inout) :: ierr
      real(real64)  :: crite
      real(real64)  :: critc
      real(real64)  :: critr
      logical, intent(in) :: hiend
      integer, intent(out) :: irev(*)
      integer  :: iselec(neig)
      integer  :: iscra2(5*lim)
      integer, intent(inout) :: incv(lim)
      integer, intent(inout) :: icv(nume + 1)
      real(real64), intent(inout) :: eigval(lim)
      real(real64)  :: basis(n*lim)
      real(real64)  :: orthobasis(n*lim + noc*n)
      real(real64)  :: ab(n*lim)
      real(real64)  :: s(lim*(lim + 1)/2)
      real(real64)  :: temps(lim*(lim + 1)/2)
      real(real64)  :: svec(lim*lim)
      real(real64), intent(inout) :: scra1(8*lim)
      real(real64), intent(inout) :: oldval(nume + 1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: rest, i, kpass, nncv, ifind, nfound, info, newstart
      real(real64) :: tol
      logical :: first, done

      save first, done, rest, i, kpass, nncv, ifind, tol, nfound, info, &
         newstart
!-----------------------------------------------
!CC   DIMENSION SVEC(LIM*(NUME+1)),EIGVAL(LIM)
!
!       !include 'mpif.h'
!       (this is for the reverse communication)

!-----------------------------------------------------------------------
!
!   on entry
!   -------
!
!   N           The order of the matrix A
!   HIEND Logical. True only if the highest eigenpairs are needed.
!   LIM  The limit on the size of the expanding Basis
!   MBLOCK Number of vectors to be targeted in each iteration.
!   NOC  Number of orthogonalization constraints (vectors)
!   NUME        The largest index of the eigenvalues wanted.
!   NIV  Starting dimension of expanding basis.
!   NEIG  Number of eigenvalues wanted.
!   ISELEC Array containg the indices of those NEIG eigenpairs.
!   CRITE       Convergence thresholds for eigenvalues, coefficients
!   CRITC,CRITR and residuals.
!   BASIS  Array with the basis vectors.
!   AB        Array with the vectors D=AB
!   S  Array keeping the symmetric matrix of the small problem.
!   TEMPS scratch array
!   SVEC Array for holding the eigenvectors of S
!   SCRA1  Srcatch array used by DSPEVX. It also holds residuals
!      and eigenvalue shifts for preconditioning.
!   ISCRA2 Integer Srcatch array used by DSPEVX.
!   INCV Srcatch array used in DSPEVX. Also used in TSTSEL and
!      NEWVEC where it holds the Indices of uNConVerged pairs
!   ICV  It contains "1" to the locations of ConVerged eigenpairs
!   OLDVAL Array keeping the previous' step eigenvalue estimates.
!
!   on exit
!   -------
!
!   EIGVAL Array containing the NUME lowest eigenvalues of the
!      the matrix A (or -A if the highest are sought).
!   Basis  On exit Basis stores the NUME corresponding eigenvectors
!   OLDVAL  On exit it stores the final differences of eigenvalues.
!   SCRA1 On exit it stores the NUME corresponding residuals.
!   NLOOPS Number of loops taken by the algorithm
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Reverse Communication
      if (irc == 1) then
!          ..Came from preconditioning
         go to 100
      else if (irc == 3) then
!          ..came from matrix vector multiply
         go to 200
      endif
!-----------------------------------------------------------------------

      i = 1
      if (nume > 0) then
         eigval(:nume) = 1.e30_real64
         icv(:nume) = 0
         i = nume + 1
      endif
      first = .true.
      kpass = niv
      nncv = kpass
!
! Decide HERE how many to restart with, not more than LIM-3.
!
      rest = max(lim/2,2*nume)
      rest = min(rest,lim - 3)

   10 continue
      ifind = min(kpass,rest)
      tol = slamch('S')
      call dcopy (ifind, eigval, 1, oldval, 1)
      call dcopy ((kpass*(kpass + 1))/2, s, 1, temps, 1)
      call dspevx ('Vectors also', 'In a range', 'Upper triangular', kpass, &
         temps, -1._real64, -1._real64, 1, ifind, tol, nfound, eigval, svec, kpass, scra1, &
         iscra2, incv, info)
      ierr = -abs(info)
      if (ierr /= 0) go to 60
!
! TeST for convergence on the absolute difference of eigenvalues between
! successive steps. Also SELect the unconverged eigenpairs and sort them
! by the largest magnitude in the last added NNCV rows of Svec.
!
      done = tstsel(kpass,nume,neig,iselec,svec,eigval,icv,crite,critc,scra1,&
         iscra2,oldval,nncv,incv)
      if (done .or. kpass>=n) go to 30
!
! Maximum size for expanding basis. Truncate basis, D, and S, Svec
! Consider the basis vectors found in TSTSEL for the newvec. KPASS=NUME
!

!          Change suggested by Anreas, March 18, 1996
      if (kpass == lim) then
!      PRINT*,'collapsing the basis: lim,nume,rest',lim,nume,rest
!        PRINT*,'myid = ', myid, '  nprocs = ', nprocs
!        PRINT*,'collapsing the basis: lim,nume,rest',lim,nume,rest,myid
!23456789012345678901234567890123456789012345678901234567890123456789012

         call ovflow (n, rest, kpass, scra1, basis, ab, s, svec, eigval)
!             CALL OVFLOW(N,NUME,KPASS,SCRA1,BASIS,AB,S,SVEC,EIGVAL)
      endif
!
! Compute and add the new residuals. NNCV is set to the number of new
! vectors that have not converged. If none, DONE=true, exit.
!
      call newvec (n, nume, lim, mblock, kpass, critr, nncv, incv, svec, eigval&
         , oldval, ab, basis, icv, scra1, scra1(lim+1), done)

      if (done) go to 30
!-----------------------------------------------------------------------
! Preconditioning the NNCV (residuals-deps x) stored in AB(kpass+1).
!          ..Robust Preconditioning (Eigenvalue shift).
!          ..Look for lowest, so Li-|eps_i|
      do i = 1, nncv
         eigval(incv(i)) = eigval(incv(i)) - scra1(lim+incv(i))
      end do
!
!          ..Change sign of eigenvalues if HIEND.
      if (hiend) call dscal (nume, -1._real64, eigval, 1)

! Use of Reverse Communication.
!
      irev(1) = nncv
      irev(2) = kpass + 1
      irev(3) = kpass + 1
      irc = 1
      if (irc == 1) return
!          ..Continue from preconditioning
  100 continue
      if (hiend) call dscal (nume, -1._real64, eigval, 1)
!
!          ..Shift the eigenvalues back to what they were.
      do i = 1, nncv
         eigval(incv(i)) = eigval(incv(i)) + scra1(lim+incv(i))
      end do
!-----------------------------------------------------------------------
!
! Orthonormalization of the previous vectors to the Basis and to any
! orthogonalization constraints. The not-yet-filled
! spaces of AB (from NEWSATRT onwards) are used for scratch.
!
      newstart = kpass*n + 1

      call mgs_nrm (n, noc + kpass, nncv, scra1(lim+1), orthobasis)
!-----------------------------------------------------------------------
! Use of Reverse Communication. Add new columns in D through matrix vector
! multiplication CALL OP(N,NNCV,BASIS(NEWSTART),AB(NEWSTART))
!
      irev(1) = nncv
      irev(2) = kpass + 1
      irev(3) = kpass + 1
      irc = 3
      return

!          ..Continue from matrix-vector multiply
  200 continue
      if (hiend) call dscal (n*nncv, -1._real64, ab(newstart), 1)
!
! Add new column in S, from the NNCV new vectors.
!
      call adds (n, lim, kpass, nncv, basis, ab, s)

      kpass = kpass + nncv
      nloops = nloops + 1

      if (nloops <= maxiter) go to 10
      ierr = ierr + 2048
      nloops = nloops - 1
      kpass = kpass - nncv
   30 continue
      do i = 1, nume
         oldval(i) = abs(oldval(i)-eigval(i))
      end do

      call multbc (n, kpass, nume, svec, scra1, basis)
      call multbc (n, kpass, nume, svec, scra1, ab)
!
! i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi
! temporarily stored in AB(NUME*N+1)
!
      do i = 1, nume
         call dcopy (n, ab((i-1)*n+1), 1, ab(nume*n+1), 1)
         call daxpy (n, (-eigval(i)),basis((i-1)*n+1), 1, ab(nume*n+1), 1)
         scra1(i) = ddot(n,ab(nume*n+1),1,ab(nume*n+1),1)
         scra1(i) = sqrt(scra1(i))
      end do
!
! Set IRC=0 for normal exit with no reverse communication
!
!     IF (MYID == 0) WRITE (6, *) 'DVDSON::NLOOPS =', NLOOPS
   60 continue
      irc = 0
      return
      end subroutine dvdrvr
end module dvdrvr_i
