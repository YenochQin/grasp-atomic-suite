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
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!   Editted by C. F. Fischer                                       5/10/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
      integer  :: irc
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
      integer  :: ierr
      real(real64)  :: crite
      real(real64)  :: critc
      real(real64)  :: critr
      logical  :: hiend
      integer  :: irev(*)
      integer  :: iselec(lim)
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


!=======================================================================
      subroutine adds(n, lim, kpass, nncv, basis, ab, s)
!=======================================================================
!     Called by: DVDSON
!
!     Calculates the new column in the S matrix. S has a
!     new row and column, but being symmetric only the new column is
!     stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i.
!
!     subroutines called:
!     DDOT, DSCAL
!
!-----------------------------------------------------------------------
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
      !USE ddot_I
      use suite_linear_algebra, only: ddot
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer , intent(in) :: lim
      integer , intent(in) :: kpass
      integer , intent(in) :: nncv
      real(real64)  :: basis(n*lim)
      real(real64)  :: ab(n*lim)
      real(real64), intent(inout) :: s(lim*(lim + 1)/2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: idstart, isstart, iv, ibstart, ibv
      real(real64) :: ss
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   -------
!   N  The order of the matrix A
!   kpass The current dimension of the expanding sub-basis
!   NNCV Number of new basis vectors.
!   Basis the basis vectors, including the new NNCV ones.
!   AB  The matrix D=AB. (with new NNCV columns)
!   on exit
!   -------
!   S           The small matrix with NNCV new columns at the last part
!-----------------------------------------------------------------------
!
! The new S is calculated by adding the new last columns
! S(new)=B^T D(new).
!
      idstart = kpass*n + 1
      isstart = kpass*(kpass + 1)/2
      do iv = 1, nncv
         ibstart = 1
         do ibv = 1, kpass + iv
            ss = ddot(n,basis(ibstart),1,ab(idstart),1)
            s(isstart+ibv) = ss
            ibstart = ibstart + n
         end do
         isstart = isstart + kpass + iv
         idstart = idstart + n
      end do

      return
      end subroutine adds


!=======================================================================
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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


!=======================================================================
      subroutine initdvd(irc, irev, n, noc, niv, nume, lim, hiend, scra1, &
         orthobasis, basis, ab, s)
!=======================================================================
!     Initializes the basis and the auxiliary arrays AB and S.
!     If not enough initial estimates exist the basis will be
!     supplemented with Lanczos vectors of the current NIV vectors.
!     eg. (b1,b2) --> (b1,Ab1,AAb1, b2,Ab2,AAb2) for a NUME of 6.
!
!     OrthoBasis is the Basis with the NOC orthogonalization
!     constraint vectors in the begining. This equialence holds:
!           equivalence(Basis(1),OrthoBasis(NOC*N+1))
!-----------------------------------------------------------------------
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
      use dinit_i
      !USE dcopy_I
      use mgs_nrm_i
      !USE dscal_I
      use adds_i
      use suite_linear_algebra, only: dcopy, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: irc
      integer, intent(in) :: n
      integer, intent(in) :: noc
      integer, intent(inout) :: niv
      integer, intent(in) :: nume
      integer, intent(in) :: lim
      logical, intent(in) :: hiend
      integer, intent(out) :: irev(*)
      real(real64)  :: scra1(*)
      real(real64)  :: orthobasis(n*(noc + lim))
      real(real64)  :: basis(n*lim)
      real(real64)  :: ab(n*lim)
      real(real64)  :: s(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ist, iend, new, kpass

      save ist, iend, new, kpass
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!       ..Reverse Communication
      if (irc == 2) go to 100
!-----------------------------------------------------------------------
!
! If no initial estimates pick one
!
      if (niv == 0) then
         call dinit (n, 1._real64/sqrt(dble(n)), basis, 1)
         niv = 1
      endif
!
! Compute AB. Also fill basis with orthonormalized ABs until enough NIVs.
!
      ist = 1
      iend = niv
      new = niv
   10 continue
      irc = 2
      irev(1) = new
      irev(2) = ist
      irev(3) = ist
      return
  100 continue
      if (iend < nume) then
         new = min(nume - iend,iend - ist + 1)
         call dcopy (n*new, ab((ist-1)*n+1), 1, basis(1+iend*n), 1)
!            ..orthonormalize OrthoBasis i.e.,
!            B-1..B-noc,B1,...,Biend,( Biend+1,...Biend+new )
         call mgs_nrm (n, noc + iend, new, scra1, orthobasis)
         ist = iend + 1
         iend = iend + new
         go to 10
      endif
      niv = iend
      !xhh print*, 'niv=',niv, 'nume=', nume
!
! Scale if HIEND for highest eigepairs
!
      if (hiend) call dscal (n*niv, -1._real64, ab, 1)
!
! Also find the small matrix S = B^TAB.
!
      kpass = 0
      call adds (n, lim, kpass, niv, basis, ab, s)

      irc = 0
      return
      end subroutine initdvd



!=======================================================================
      subroutine multbc(n, k, m, c, temp, b)
!=======================================================================
!     called by: DVDRVR
!
!     Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
!     Used for collapsing the expanding basis to current estimates,
!     when basis becomes too large, or for returning the results back
!
!       Subroutines called
!       DINIT, DGEMV, DCOPY
!-----------------------------------------------------------------------
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
      !USE dgemv_I
      !USE dcopy_I
      use suite_linear_algebra, only: dcopy, dgemv
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: k
      integer, intent(in) :: m
      real(real64)  :: c(k*m)
      real(real64)  :: temp(m)
      real(real64)  :: b(n*k)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: irow
!-----------------------------------------------
!-----------------------------------------------------------------------
      do irow = 1, n
         call dgemv ('Transp', k, m, 1._real64, c, k, b(irow), n, 0._real64, temp, 1)
         call dcopy (m, temp, 1, b(irow), n)
      end do

      return
      end subroutine multbc


!=======================================================================
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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


!=======================================================================
      subroutine ovflow(n, nume, kpass, scra1, basis, ab, s, svec, eigval)
!=======================================================================
!     Called by: DVDRVR
!     Called when the upper limit (LIM) has been reached for the basis
!     expansion (by KPASS). The basis is truncated B'=BC and
!     similarly the array D'=DC. The new S is computed as
!     S'(i,j)=l(i)delta(i,j) where l(i) eigenvalues, and delta of
!     Kronecker, i,j=1,NUME. The new eigenvectors of the small matrix
!     are the unit vectors.
!
!     Subroutines called:
!     DCOPY, DINIT, MULTBC
!-----------------------------------------------------------------------
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
      use multbc_i
      use dinit_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: nume
      integer, intent(inout) :: kpass
!CFF   ... add dimensionto SCRA1
      real(real64)  :: scra1(nume)
      real(real64)  :: basis(n*kpass)
      real(real64)  :: ab(n*kpass)
      real(real64), intent(out) :: s((kpass*(kpass + 1))/2)
      real(real64)  :: svec(kpass*nume)
      real(real64) , intent(in) :: eigval(kpass)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ind, icur, i
!-----------------------------------------------
!   on entry
!   -------
!   NUME The largest index of eigenvalues wanted.
!   SVEC  the kpass eigenvectors of the smaller system solved
!   EIGVAL  the eigenvalues of this small system
!   on exit
!   -------
!   Basis
!   AB
!   S   The new small matrix to be solved.
!-----------------------------------------------------------------------
! Truncate  the basis and the AB array.
!
      call multbc (n, kpass, nume, svec, scra1, basis)
      call multbc (n, kpass, nume, svec, scra1, ab)
!
! calculation of the new upper S=diag(l1,...,l_NUME) and
! its matrix Svec of eigenvectors (e1,...,e_NUME)
!
      call dinit ((nume*(nume + 1))/2, 0._real64, s, 1)
      call dinit (nume*nume, 0._real64, svec, 1)
      ind = 0
      icur = 0
      do i = 1, nume
         s(ind+i) = eigval(i)
         svec(icur+i) = 1
         icur = icur + nume
         ind = ind + i
      end do

      kpass = nume

      return
      end subroutine ovflow


!=======================================================================
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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


!=======================================================================
      subroutine mgs_nrm(n, kp, new, scra, b)
!=======================================================================
!       Orthogonalizis the new vectors in B (from kp+1...kp+new)
!     to the previous kp B vectors and to themselves
!     using Modified Gram Schmidt. Then normalizes them.
!       The procedure is repeated twice.
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:12:31   2/12/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!     USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !use dgemv_I
      !use ddot_I
      !use dscal_I
      use suite_linear_algebra, only: ddot, dgemv, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer , intent(in) :: kp
      integer, intent(in) :: new
      real(real64), intent(inout) :: scra(new)
      real(real64), intent(inout) :: b((kp + new)*n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, newstart, kcur, k, jcur, j, mm
      real(real64) :: dnm
!-----------------------------------------------

!
! MODIFIED GRAM-SCHMIDT  (twice)
!
        !if (myid.eq.0) then
      do i = 1, 2

         newstart = kp*n + 1
!
!  First record contribution from the kp vectors to the new ones.
!
         kcur = 1
         do k = 1, kp
            jcur = newstart
            call dgemv ('T', n, new, 1._real64, b(jcur), n, b(kcur), 1, 0._real64, scra, &
               1)
            do j = 1, new
!           call daxpy(N,-scra(j),B(kcur),1,B(jcur),1)
               do mm = 0, n - 1
                  b(jcur+mm) = b(jcur+mm) - scra(j)*b(kcur+mm)
               end do
               jcur = jcur + n
            end do
            kcur = kcur + n
         end do
!
!  Then orthogonalize the new ones among themselves.
!
         do k = 1, new
            jcur = kcur + n
!  The current vector should be normalized
!
            dnm = ddot(n,b(kcur),1,b(kcur),1)
            dnm = sqrt(dnm)
            call dscal (n, 1/dnm, b(kcur), 1)
!
            ! No trailing vectors: avoid forming B one element past its end.
            if (k < new) then
            call dgemv ('T', n, new - k, 1._real64, b(jcur), n, b(kcur), 1, 0._real64, &
               scra, 1)
            end if
            do j = k + 1, new
!              call daxpy(N,-scra(j-k),B(kcur),1,B(jcur),1)
               do mm = 0, n - 1
                  b(jcur+mm) = b(jcur+mm) - scra(j-k)*b(kcur+mm)
               end do
               jcur = jcur + n
            end do
            kcur = kcur + n
         end do
      end do
        !end if

!        call MPI_BCAST(B,(kp+New)*n,MPI_DOUBLE_PRECISION,0,
!     :                 MPI_COMM_WORLD,ierr)

      return
      end subroutine mgs_nrm
