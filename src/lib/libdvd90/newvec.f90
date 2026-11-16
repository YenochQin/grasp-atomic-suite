module newvec_i
    implicit none (type, external)
    private
    public :: newvec
contains
      subroutine newvec(n, nume, lim, mblock, kpass, critr, nncv, incv, svec, &
         eigval, oldval, ab, basis, icv, scra1, epsil, done)
!=======================================================================
!
!     Called by: DVDRVR
!
!     It calculates the new expansion vectors of the basis.
!     For each one of the vectors in INCV starting with the largest
!     megnitude one, calculate its residual Ri= DCi-liBCi and check
!     the ||Ri|| for convergence. If it is converged do not add it
!     but look for the immediate largest coefficient and its vector.
!     The above procedure continues until MBLOCK vectors have been
!     added to the basis, or the upper limit has been encountered.
!     Thus only  the required MBLOCK residuals are computed.
!     For robust Preconditioning, each residual is then modified as
!               \delta_eps x - res
!     where \delta_eps is an estimate of the eigenvalue error.
!     A different estimate is used for the eigenvalue shift.
!
!     Subroutines called:
!     DNRM2, DGEMV
!
!-----------------------------------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!     USE MPI_C,           ONLY: MYID
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use suite_linear_algebra, only: daxpy, ddot, dgemv
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer , intent(in) :: nume
      integer , intent(in) :: lim
      integer , intent(in) :: mblock
      integer, intent(in) :: kpass
      integer , intent(inout) :: nncv
      real(real64) , intent(in) :: critr
      logical , intent(out) :: done
      integer , intent(inout) :: incv(nume)
      integer, intent(inout) :: icv(nume)
      real(real64)  :: svec(lim*nume)
      real(real64) , intent(in) :: eigval(lim)
      real(real64) , intent(in) :: oldval(nume)
      real(real64)  :: ab(n*lim)
      real(real64)  :: basis(n*lim)
      real(real64) , intent(inout) :: scra1(lim)
      real(real64) , intent(out) :: epsil(lim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: newstart, nadded, icvc, limadd, icur, i, indx
      real(real64) :: sqres, gaplow, gap, gapup, dl, res, rhseps
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   --------
!   N           The order of the matrix A
!   NUME        The largest index of the eigenvalues wanted.
!   LIM         The limit on the size of the expanding Basis
!   MBLOCK Maximum number of vectora to enter the basis
!   KPASS the current dimension of the expanding basis
!   CRITR Convergence threshold for residuals
!   NNCV Number of Non ConVerged pairs (MBLOCK will be targeted)
!   INCV        Index to the corresponding SVEC columns of these pairs.
!   SVEC,EIGVAL Arrays holding the eigenvectors and eigenvalues of S
!   AB  Array with the vectors D=AB
!   BASIS the expanding basis having kpass vectors
!   ICV  Index of converged eigenpairs (ICV(i)=1 <=>i converged)
!
!   on exit
!   -------
!   NNCV The number of vectors finally added to the basis.
!   BASIS  The eigenector estimated reside at the end.
!   AB          The right hand sides for preconditioning.
!   EPSIL       The eigenvalue shifts for preconditioning,
!   ICV  Index of converged eigenpairs (updated)
!   DONE  logical, if covergance has been reached.
!      the Basis should be collapsed to current approximations.
!-----------------------------------------------------------------------
!       !include 'mpif.h'
      done = .false.
      newstart = kpass*n + 1
      nadded = 0
      icvc = 0
      limadd = min(lim,mblock + kpass)
      icur = newstart
!
! Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors.
!
      do i = 1, nncv
         indx = incv(i)
!
!          ..Compute b = BASIS*Svec_indx
!        ..Compute d = AB*Svec_indx
!        ..Daxpy d'= d - eigval b  gives the residual

         call dgemv ('N', n, kpass, 1._real64, basis, n, svec((indx-1)*kpass+1) &
            , 1, 0._real64, basis(icur), 1)
         call dgemv ('N', n, kpass, 1._real64, ab, n, svec((indx-1)*kpass+1), 1, &
            0._real64, ab(icur), 1)
         call daxpy (n, (-eigval(indx)),basis(icur), 1, ab(icur), 1)
!
!        ..Compute the norm of the residual
!        ..and check for convergence
!
         sqres = ddot(n,ab(icur),1,ab(icur),1)
         scra1(indx) = sqrt(sqres)

!        IF (MYID == 0) WRITE (6, '(A11,F22.16,I2,A10,F19.16)') ' EIGVAL(i) ', &
!           EIGVAL(INDX), INDX, ' Res.Norm ', SCRA1(INDX)

         if (scra1(indx) < critr) then
!           ..Converged,do not add. Go for next non converged one
            icvc = icvc + 1
            icv(indx) = 1
            if (icvc < nncv) cycle
!           ..All have converged.
       !print*, 'converged by critr'
       !print*, 'converged by critr',myid
            done = .true.
            return
         else
!           ..Not converged. Consider it for preconditioning
!          ---------------  ROBUST MODIFICATION ---------------------
!           ..Daxpy d'= d'- \delta_eps b  gives the rhs for precond.
!           ..It is stored in AB.
!             ..Deps are also stored in EPSIL
!
            if (indx == 1) then
               gaplow = 1.0e99_real64
               gap = abs(eigval(2)-eigval(1))
!                 print*, 'h1,h2:',eigval(1),eigval(2)
            else
               gaplow = abs(eigval(indx)-eigval(indx-1))
               gapup = abs(eigval(indx+1)-eigval(indx))
               gap = min(gaplow,gapup)
            endif
            dl = abs(oldval(indx)-eigval(indx))
            res = scra1(indx)

            if (gap > res) then
!                 EPSIL(indx) = min(dl,res,gaplow)
               rhseps = sqrt(dl*res)
            else
!                EPSIL(indx) = min( res, gaplow)
               rhseps = min(dl,sqrt(dl))
            endif

            epsil(indx) = 0._real64
            call daxpy (n, (-rhseps), basis(icur), 1, ab(icur), 1)

!          -------------- END OF ROBUST MODIFICATIONS -----------------
!
            nadded = nadded + 1
            incv(nadded) = indx
            if (nadded + kpass == limadd) exit
!           ..More to be added in the block
            icur = icur + n
         endif
      end do

      nncv = nadded

      return
      end subroutine newvec
end module newvec_i
