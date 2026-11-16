module dvdson_i
    implicit none (type, external)
    private
    public :: dvdson
contains
      subroutine dvdson(irc, irev, n, lim, noc, ilow, ihigh, iselec, niv, &
         mblock, crite, critc, critr, maxiter, work, iwrsz, iwork, iiwsz, hiend&
         , nloops, ierr)
!=======================================================================
!
!     Author: Andreas Stathopoulos, Charlotte F. Fischer
!
!     Computer Science Department
!     Vanderbilt University
!     Nashville, TN 37212
!     andreas@vuse.vanderbilt.edu
!     cff@vuse.vanderbilt.edu
!
!     June 1995
!
!     Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer
!
!     References: A Stathopoulos and C F Fischer,
!          Comp. Phys. Commun. 79 (1994) 268-290.
!
!                   A Stathopoulos, Y. Saad and C F Fischer,
!                   J. Computational and applied Mathematics (in print).
!
!     DVDSON is a Fortran77 program that finds a few selected
!     eigenvalues and their eigenvectors at either end of spectrum of
!     a large, symmetric (and usually sparse) matrix, denoted as A.
!     The matrix A is only referenced externally through the user
!     supplied routine which implements a block matrix-vector
!     operation(see below) and the user supplied preconditioning
!       routine (also perfomed externally). Either the range of the
!       eigenvalues wanted or an array of the indices of selected ones
!       can be specified.
!     DVDSON is a front-end routine for setting up arrays, and initial
!     guess. It also performs detailed error checking.
!     DVDRVR is the driver routine that implements a version of the
!     Davidson algorithm. The characteristics of this version are:
!      o  Use of REVERSE COMMUNICATION to perform the preconditioning
!         and the matrix vector multiply.
!      o  All arrays used by the program are stored in MEMORY.
!      o  BLOCK method (many vectors may be targeted per iteration.)
!      o  Eigenvectors are targeted in an optimum way without
!         the need to compute all unconverged residuals,
!      o  Uses two Modified Gram-Schmidt passes for orthogonality loss.
!      o  Finds HIGHEST eigenpairs by using the negative of the A.
!      o  Finds SELECTED eigenpairs specified by the user.
!      o  If not enough initial vectors are available (NIV<NUME)
!           it builds the first NUME-NIV Lanczos from those given.
!      o  The user can provide STOPPING CRITERIA for eigenvalues,
!         and residuals and coefficients and he  can CONTROL block size.
!      o  On exit INFORMATION is given about the convergence status
!         of eigenpairs.
!
!     The program consists of the following routines:
!     DVDSON, INITDVD, DVDRVR, ADDS, TSTSEL, MULTBC, OVFLOW,
!       NEWVEC, MGS_NRM
!
!     It also calls some basic BLAS routines:
!     DCOPY, DSCAL, DDOT, DAXPY, IDAMAX, DGEMV, DINIT
!
!     For solving the small eigenproblem, the routine DSPEVX from
!     LAPACK is used. DSPEVX is obtainable from NETLIB, together
!     with a series of subroutines that it calls.
!
!     All the routines have IMPLICIT REAL*8          (A-H,O-Z)
!
!-----------------------------------------------------------------------
!   Editted by C. F. Fischer                                       5/10/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use initdvd_i
      use dvdrvr_i
      !USE dscal_I
      !USE dcopy_I
      use suite_linear_algebra, only: dcopy, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: irc
      integer, intent(in) :: n
      integer, intent(in) :: lim
      integer, intent(in) :: noc
      integer, intent(inout) :: ilow
      integer, intent(inout) :: ihigh
      integer, intent(inout) :: niv
      integer, intent(in) :: mblock
      integer  :: maxiter
      integer, intent(in) :: iwrsz
      integer, intent(in) :: iiwsz
      integer, intent(inout) :: nloops
!      INTEGER, INTENT(INOUT) :: IERR
      integer, intent(inout) :: ierr
      real(real64)  :: crite
      real(real64)  :: critc
      real(real64)  :: critr
      logical, intent(inout) :: hiend
      integer, intent(inout) :: irev(*)
      integer, intent(inout) :: iselec(lim)
      integer, intent(inout) :: iwork(iiwsz)
      real(real64), intent(inout) :: work(iwrsz)
      real(kind(0.0_real64)) :: ddot
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: neig, i, nume, iortho, ibasis, ieigval, iab, is, itemps, isvec&
         , iscra1, ioldval, iscra2, iscra3, iicv
      save
!      SAVE NEIG, I, NUME, IORTHO, IBASIS, IEIGVAL, IAB, IS, ITEMPS, ISVEC, &
!         ISCRA1, IOLDVAL, ISCRA2, ISCRA3, IICV
!-----------------------------------------------
!
!       (this is for the reverse communication)

!-----------------------------------------------------------------------
!  (Important to the following is the concept of NUME, the distance of
!   the index of the eigenpair wanted which is farthest from the
!   extremes,i.e.,
!      if  lowest  eigepairs i1<i2<...<ik are wanted, NUME=ik
!      if highest eigenpairs i1<i2<...<ik are wanted, NUME=N-i1+1
!   where i1,...,ik are the indices of the wanted eigenpairs.
!   Obviously, NUME.GE.(No. of EiGenpairs wanted). )

!   on entry
!   -------
!   IRC         REVERSE COMMUNICATION integer.
!               IRC=0   --->>  Exit
!               IRC=1   --->>  Preconditioning
!               IRC=2   --->>  Matrix vector for initial setups.
!               IRC=3   --->>  Matrix vector in DVDRVR
!        NOTE Usually in reverse communication the routine called
!        ---- gives out two arrays w1 and w2. The first is used
!      as input for calculation Aw1->w2, or solve Mw2=w1.
!      Then these are passed back to the routine.
!      TO avoid extra work space, AND to be able to
!      perform block operations (i.e., where w1 is m blocks)
!      we use the same WORK array, where BASIS and AB are
!      stored in. An integer array with 6 integers is passed out,
!      except IRC: IREV.
!      IREV(1)=NB Simply NB operations must be performed
!                 from either the preconditioner or the matvec.
!          Also the number of targeted eigenpairs.
!      IREV(2)=iw1
!      IREV(3)=iw2
!              A*WORK(iw1:iw1+NB*N) --> WORK(iw2:iw1+NB*N)
!              Solve NB rhsides:
!              M*WORK(iw2:iw1+NB*N) = WORK(iw1:iw1+NB*N)
!              Obviously: iw1=ibasis+("indexB"-1)*N
!                  iw2=iab   +("indexD"-1)*N
!      for the matvec, while for the preconditioner it is
!      in reverse order.
!      IREV(4)=ieigval  (To denote where the eigenvalues are)
!      IREV(5)=iscra3   (where the eigenvalue index starts)
!              Simply the targeted eigenvalues are:
!           do i=1,NB
!          indx = IWORK(iscra3+i-1)
!          L_i  = WORK(ieigval+indx-1)
!                          solve (M-L_i) w2 = w1
!       enddo
!               The following is additional information that advanced
!      users may need.
!      IREV(6)=iscra1   (where the residuals of the corresponding
!             eigenpairs are. According to the above index)
!      IREV(7)=iscra1+LIM  this will hold the \delta eps shifts
!      from robust preconditioning
!
!   N           Order of the matrix.
!   LIM  The upper limit on the dimension of the expanding basis.
!      NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
!      only for LIM=NUME=N. The choice of LIM depends on the
!      available workspace (see below). If the space is
!      available it is preferable to have a large LIM, but not
!      larger than NUME$+$40.
!   NOC  Number of Orthogonalization Constraints (vectors).
!             These NOC vectors are given on input in the first
!      N*NOC elements of WORK.
!   ILOW The index of the lowest eigepair to be computed. If
!      (ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs
!      to be computed should be contained in array ISELEC.
!      (Modified on exit).
!   IHIGH  The index of the highest eigenpair to be computed.
!      Considered ONLY when ILOW is in the range
!      (0.LT.ILOW.LE.N). (Modified on exit).
!   ISELEC Array of size LIM holding the user specified indices
!      for the eigenpairs to be computed. Considered only when
!      (ILOW.LE.0).or.(ILOW.GT.N). The indices are read from
!      the first position until a non positive integer is met.
!         Example: if N=500, ILOW=0, and ISELEC(1)=495,
!         ISELEC(2)=497, ISELEC(3)=-1, the program will find
!         2 of the highest eigenpairs, pairs 495 and 497.
!      Any order of indices is acceptable (Modified on exit).
!   NIV  Number of Initial Vector estimates provided by the user.
!      If LIM > NIV > 0 then the NIV columns of size N of WORK
!      coming after the NOC constraint vectors,  should contain
!               the estimates (see below).
!               If NIV = 0, the vector (1,1,1,...,1)T is chosen.
!   MBLOCK Number of vectors to be targeted in each iteration.
!      1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold.
!      Large block size reduces the number of iterations
!      (matrix acceses) but increases the matrix-vector
!      multiplies. It should be used when the matrix accese
!      is expensive (disc, recomputed or distributed).
!   CRITE Convergence threshold for eigenvalues.
!      If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted
!      eigenvalues, convergence is signaled.
!   CRITC       Convergence threshold for the coefficients of the last
!               added basis vector(s). If all of those corresponding
!               to unconverged eigenpairs are less than CRITC convergence
!               is signaled.
!   CRITR Convergence threshold for residual vector norms. If
!      all the residual norms ||Ax_i-l_ix_i|| of the targeted
!      x_i are less than CRITR convergence is signaled.
!      If ANY of the criteria are satisfied the algorithm stops
!   MAXITER Upper bound on the number of iterations of the
!      algorithm. When MAXITER is exceeded the algorithm stops.
!      A typical MAXITER can be MAX(200,NUME*40), but it can
!      be increased as needed.
!   WORK Real array of size IWRSZ. Used for both input and output
!      If NIV is in (1.LE.(NIV).LE.(LIM)), on input, WORK
!      must have the NIV initial estimates. These NIV N-element
!      vectors start from WORK(1) and continue one after the
!      other. They must form an orthonormal basis.
!   IWRSZ The size of the real workspace. It must be at least as
!      large as:
!
!      Not this any more:
!           (2*LIM+NOC)*N + LIM*LIM + (NUME+10)*LIM + NUME
!      New workspace:
!     --> this       (2*LIM+NOC)*N + LIM*LIM*2 +      11*LIM + NUME
!
!   IWORK  Integer work array of size IIWSZ. Used as scrath array
!      for indices and for use in the LAPACK routines.
!   IIWSZ The size of the integer workspace. It must be at least
!      as large as:
!                    6*LIM + NUME
!
!      If LIM or NUME needs to be increased, the space should
!               also be increased accordingly. For given IWRSZ and
!      IIWSZ one can calculate how big a problem one can
!      solve (LIM,NUME).
!
!   on exit
!   -------
!   WORK(NOC*N+1)
!      The first NUME*N locations contain the approximations to
!      the NUME extreme eigenvectors. If the lowest eigenpairs
!      are required, (HIEND=false), eigenvectors appear in
!      ascending order, otherwise (HIEND=false), they appear in
!      descending order. If only some are requested, the order
!      is the above one for all the NUME extreme eigenvectors,
!      but convergence has been reached only for the selected
!      ones. The rest are the current approximations to the
!      non-selected eigenvectors.
!   WORK(NOC*N+NUME*N+1)
!      The next NUME locations contain the approximations to
!      the NUME extreme eigenvalues, corresponding to the above
!      NUME eigenvectors. The same ordering and convergence
!      status applies here as well.
!   WORK(NOC*N+NUME*N+NUME+1)
!      The next NUME locations contain the corresponding values
!      of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of
!      the last step of the algorithm.
!   WORK(NOC*N+NUME*N+NUME+NUME+1)
!      The next NUME locations contain the corresponding
!      residual norms of the NUME above eigenvectors, of the
!      last step.
!   NLOOPS      Number of iterations (accesses to matrix)
!
!   HIEND Logical. If .true. on exit the highest eigenpairs are
!      found in descending order. Otherwise, the lowest
!      eigenpairs are arranged in ascending order.
!   IERR An integer denoting the completions status:
!      IERR = 0  denotes normal completion.
!      IERR = -k  denotes error in DSPEVX (k eigenpairs
!         not converged)
!      0<IERR<=2048 denotes some inconsistency as follows:
!      If (INT( MOD(IERR,  2)/1  ) N < LIM
!      If (INT( MOD(IERR,  4)/2  ) LIM < 1
!      If (INT( MOD(IERR,  8)/4  ) ISELEC(1)<1, and no range specified
!      If (INT( MOD(IERR, 16)/8  ) IHIGH > N (in range or ISELEC)
!      If (INT( MOD(IERR, 32)/16 ) IHIGH < ILOW (Invalid range)
!      If (INT( MOD(IERR, 64)/32 ) NEIG >= LIM (Too many wanted)
!      If (INT( MOD(IERR,128)/64 ) Probable duplication in ISELEC
!      If (INT( MOD(IERR,256)/128) NUME >= LIM (max eigen very far)
!      If (INT( MOD(IERR,512)/256) MBLOCK is out of bounds
!      If (INT( MOD(IERR,1024)/512) IWRSZ or IIWSZ is not enough
!      If (INT( MOD(IERR,2048)/1024) Orthogonalization Failed
!      If (INT( MOD(IERR,4096)/2048) NLOOPS > MAXITER
!      If (INT( MOD(IERR,8192)/4096) Not proper number of NIV
!
!      The program will also print an informative message to
!      the standard output when NIV is not proper
!-----------------------------------------------------------------------
!     :    ILOW,IHIGH,ISELEC,NIV,MBLOCK,
!     :    CRITE,CRITC,CRITR,MAXITER,
!     :    WORK,IWRSZ,IWORK,IIWSZ,
!     :    HIEND,NLOOPS,IERR)
!-----------------------------------------------------------------------
!       Reverse Communication
!-----------------------------------------------------------------------
      if (irc==1 .or. irc==3) then
!          ..(=1)preconditioning. Go into DVDRVR
!          ..(=3)or Matrix vector in DVDRVR
         go to 200
      else if (irc == 2) then
!          ..matrix vector for setting up.
         go to 100
      endif
!-----------------------------------------------------------------------
!
! Checking user input errors, and setting up the problem to solve.
!
      ierr = 0
      if (lim > n) ierr = ierr + 1
      if (lim <= 0) ierr = ierr + 2

      hiend = .false.

      if (ilow<=0 .or. ilow>n) then
!          ..Look for user choice of eigenpairs in ISELEC
         if (iselec(1) <= 0) then
!             ..Nothing is given in ISELEC
            ierr = ierr + 4
         else
!             ..Find number of eigenpairs wanted, and their
!           ..min/max indices
            neig = 1
            ilow = iselec(1)
            ihigh = iselec(1)
            do i = 2, lim
               if (iselec(i) <= 0) exit
               ilow = min(ilow,iselec(i))
               ihigh = max(ihigh,iselec(i))
               neig = neig + 1
            end do
!           ..Check if a very large index is asked for
            if (ihigh > n) ierr = ierr + 8
         endif
      else
!          ..Look for a range between ILOW and IHIGH
!          ..Invalid range. IHIGH>N
         if (ihigh > n) ierr = ierr + 8
         neig = ihigh - ilow + 1
!          ..Invalid range. IHIGH<ILOW
         if (neig <= 0) ierr = ierr + 16
         if (neig > lim) then
!           ..Not enough Basis space. Increase LIM or decrease NEIG
            ierr = ierr + 32
         else
!           ..Fill in the ISELEC with the required indices
            do i = 1, neig
               iselec(i) = ilow + i - 1
            end do
         endif
      endif

      if (ierr /= 0) return

      nume = ihigh
!       ..Identify if few of the highest eigenpairs are wanted.
      if (ilow + ihigh - 1 > n) then
         hiend = .true.
         nume = n - ilow + 1
!          ..Change the problem to a minimum eipenpairs one
!          ..by picking the corresponding eigenpairs on the
!          ..opposite side of the spectrum.
         i = 1
         if (neig > 0) then
            iselec(:neig) = n - iselec(:neig) + 1
            i = neig + 1
         endif
      endif
!       ..duplications in ISELEC
      if (neig > nume) ierr = ierr + 64
!       ..Not enough Basis space. Increase LIM or decrease NUME
      if (nume>lim .or. nume==lim .and. nume/=n) ierr = ierr + 128
!     ..Size of Block out of bounds
      if (mblock<1 .or. mblock>neig) ierr = ierr + 256

!       ..Check for enough workspace for Dvdson
      if (iwrsz<lim*(2*n + lim + (nume + 10)) + nume + noc*n .or. iiwsz<6*lim+&
         nume) ierr = ierr + 512

      if (niv > lim) then
!        ..Check number of initial estimates NIV is lower than LIM.
         write (6, *) 'WARNING: Too many initial estimates.?'
         write (6, *) 'The routine needs at most:', lim
         ierr = ierr + 4096
      else if (niv < nume) then
!        ..check if enough initial estimates.
         write (6, *) 'WARNING: Not enough initial estimates'
         write (6, *) nume - niv, ' Lanczos vectors will be added'
      endif

      if (ierr /= 0) return
!
! Assigning space for the real work arrays
!
      iortho  = 1
      ibasis  = iortho  + n*noc
      ieigval = ibasis  + n*lim
      iab     = ieigval + lim
      is      = iab     + n*lim
      itemps  = is      + lim*(lim + 1)/2
      isvec   = itemps  + lim*(lim + 1)/2
!CC   iscra1  = iSvec   + LIM*(NUME+1)
      iscra1  = isvec   + lim*lim
      ioldval = iscra1  + 8*lim
!
! Assigning space for the integer work arrays
!
      iscra2 = 1
      iscra3 = iscra2 + 5*lim
      iicv = iscra3 + lim
!
! Initialize tha basis, the AB, the S.
!

  100 continue
      call initdvd (irc, irev, n, noc, niv, nume + 1, lim, hiend, work(iscra1)&
         , work(iortho), work(ibasis), work(iab), work(is))
!       ----------------------------------------------------------------
!       ..Reverse Communication for possible matrix vector
      if (irc == 2) then
         irev(2) = ibasis + (irev(2)-1)*n
         irev(3) = iab + (irev(3)-1)*n
         return
      endif
!       ----------------------------------------------------------------
!
! Call main driver routine.
!
      nloops = 1

  200 continue
      call dvdrvr (irc, irev, n, hiend, lim, mblock, noc, nume, niv, neig, &
         iselec, crite, critc, critr, maxiter, work(ieigval), work(ibasis), &
         work(iortho), work(iab), work(is), work(itemps), work(isvec), work(&
         iscra1), iwork(iscra2), iwork(iscra3), iwork(iicv), work(ioldval), &
         nloops, ierr)
!       ----------------------------------------------------------------
!       some Reverse Communication
      if (irc == 1) then
!        ..Preconditioning
         irev(2) = iab + (irev(2)-1)*n
         irev(3) = ibasis + (irev(3)-1)*n
         irev(4) = ieigval
         irev(5) = iscra3
         irev(6) = iscra1
         irev(7) = iscra1 + lim
         return
      else if (irc == 3) then
!        ..Matrix-vector
         irev(2) = ibasis + (irev(2)-1)*n
         irev(3) = iab + (irev(3)-1)*n
         return
      endif
!       ----------------------------------------------------------------

      if (hiend) call dscal (nume, -1._real64, work(ieigval), 1)
!
! -Copy the eigenvalues after the eigenvectors
! -Next, copy the difference of eigenvalues between the last two steps
! -Next, copy the residuals for the first NUME estimates
!
      call dcopy (nume, work(ieigval), 1, work(ibasis+n*nume), 1)
      call dcopy (nume, work(ioldval), 1, work(ibasis+(n+1)*nume), 1)
      call dcopy (nume, work(iscra1), 1, work(ibasis+(n+2)*nume), 1)
!
! Set IRC=0 for normal exit with no reverse communication
!
      irc = 0
      return
      end subroutine dvdson
end module dvdson_i
