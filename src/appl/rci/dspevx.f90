!***********************************************************************
!                                                                      *

      subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w&
         , z, ldz, work, iwork, ifail, info)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use cons_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_i
      use dlamch_i
      use dlansp_i
      use dcopy_i
      use dopgtr_i
      use dopmtr_i
      use dscal_i
      use dsptrd_i
      use dstebz_i
      use dstein_i
      use dsteqr_i
      use dsterf_i
      use dswap_i
      use xerbla_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n
      integer  :: il
      integer  :: iu
      integer  :: m
      integer  :: ldz
      integer  :: info
      real(double), intent(in) :: vl
      real(double), intent(in) :: vu
      real(double), intent(in) :: abstol
      character  :: jobz
      character  :: range
      character  :: uplo
      integer, dimension(*)  :: iwork
      integer, dimension(*)  :: ifail
      real(double), dimension(*)  :: ap
      real(double), dimension(*)  :: w
      real(double), dimension(ldz,*)  :: z
      real(double), dimension(*)  :: work
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, iinfo, imax, indd, inde, indee, indibl, indisp, indiwo, &
         indtau, indwrk, iscale, itmp1, j, jj, nsplit
      real(double) :: abstll, anrm, bignum, eps, rmax, rmin, safmin, sigma, &
         smlnum, tmp1, vll, vuu
      logical :: alleig, indeig, valeig, wantz
      character :: order
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      intrinsic min, sqrt
!-----------------------------------------------
!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSPEVX computes selected eigenvalues and, optionally, eigenvectors
!  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
!  can be selected by specifying either a range of values or a range of
!  indices for the desired eigenvalues.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  RANGE   (input) CHARACTER*1
!          = 'A': all eigenvalues will be found;
!          = 'V': all eigenvalues in the half-open interval (VL,VU]
!                 will be found;
!          = 'I': the IL-th through IU-th eigenvalues will be found.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!
!          On exit, AP is overwritten by values generated during the
!          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!          and first superdiagonal of the tridiagonal matrix T overwrite
!          the corresponding elements of A, and if UPLO = 'L', the
!          diagonal and first subdiagonal of T overwrite the
!          corresponding elements of A.
!
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues. VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute error tolerance for the eigenvalues.
!          An approximate eigenvalue is accepted as converged
!          when it is determined to lie in an interval [a,b]
!          of width less than or equal to
!
!                  ABSTOL + EPS *   max( |a|,|b| ) ,
!
!          where EPS is the machine precision.  If ABSTOL is less than
!          or equal to zero, then  EPS*|T|  will be used in its place,
!          where |T| is the 1-norm of the tridiagonal matrix obtained
!          by reducing AP to tridiagonal form.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!          If this routine returns with INFO>0, indicating that some
!          eigenvectors did not converge, try setting ABSTOL to
!          2*DLAMCH('S').
!
!          See "Computing Small Singular Values of Bidiagonal Matrices
!          with Guaranteed High Relative Accuracy," by Demmel and
!          Kahan, LAPACK Working Note #3.
!
!  M       (output) INTEGER
!          The total number of eigenvalues found.  0 <= M <= N.
!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the selected eigenvalues in ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
!          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!          contain the orthonormal eigenvectors of the matrix A
!          corresponding to the selected eigenvalues, with the i-th
!          column of Z holding the eigenvector associated with W(i).
!          If an eigenvector fails to converge, then that column of Z
!          contains the latest approximation to the eigenvector, and the
!          index of the eigenvector is returned in IFAIL.
!          If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are
!          supplied in the array Z; if RANGE = 'V', the exact value of M
!          is not known in advance and an upper bound must be used.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
!
!  IWORK   (workspace) INTEGER array, dimension (5*N)
!
!  IFAIL   (output) INTEGER array, dimension (N)
!          If JOBZ = 'V', then if INFO = 0, the first M elements of
!          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge.
!          If JOBZ = 'N', then IFAIL is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, then i eigenvectors failed to converge.
!                Their indices are stored in array IFAIL.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = lsame(jobz,'V')
      alleig = lsame(range,'A')
      valeig = lsame(range,'V')
      indeig = lsame(range,'I')
!
      info = 0
      if (.not.(wantz .or. lsame(jobz,'N'))) then
         info = -1
      else if (.not.(alleig .or. valeig .or. indeig)) then
         info = -2
      else if (.not.(lsame(uplo,'L') .or. lsame(uplo,'U'))) then
         info = -3
      else if (n < 0) then
         info = -4
      else if (valeig .and. n>0 .and. vu<=vl) then
         info = -7
      else if (indeig .and. il<1) then
         info = -8
      else if (indeig .and. (iu<min(n,il) .or. iu>n)) then
         info = -9
      else if (ldz<1 .or. wantz .and. ldz<n) then
         info = -14
      endif
!
      if (info /= 0) then
         call xerbla ('DSPEVX', (-info))
         return
      endif
!
!     Quick return if possible
!
      m = 0
      if (n == 0) return
!
      if (n == 1) then
         if (alleig .or. indeig) then
            m = 1
            w(1) = ap(1)
         else
            if (vl<ap(1) .and. vu>=ap(1)) then
               m = 1
               w(1) = ap(1)
            endif
         endif
         if (wantz) z(1,1) = one
         return
      endif
!
!     Get machine constants.
!
      safmin = dlamch('Safe minimum')
      eps = dlamch('Precision')
      smlnum = safmin/eps
      bignum = one/smlnum
      rmin = sqrt(smlnum)
      rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
!
!     Scale matrix to allowable range, if necessary.
!
      iscale = 0
      abstll = abstol
      if (valeig) then
         vll = vl
         vuu = vu
      endif
      anrm = dlansp('M',uplo,n,ap,work)
      if (anrm>zero .and. anrm<rmin) then
         iscale = 1
         sigma = rmin/anrm
      else if (anrm > rmax) then
         iscale = 1
         sigma = rmax/anrm
      endif
      if (iscale == 1) then
         call dscal ((n*(n + 1))/2, sigma, ap, 1)
         if (abstol > 0) abstll = abstol*sigma
         if (valeig) then
            vll = vl*sigma
            vuu = vu*sigma
         endif
      endif
!
!     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
      indtau = 1
      inde = indtau + n
      indd = inde + n
      indwrk = indd + n
      call dsptrd (uplo, n, ap, work(indd), work(inde), work(indtau), iinfo)
!
!     If all eigenvalues are desired and ABSTOL is less than or equal
!     to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails
!     for some eigenvalue, then try DSTEBZ.
!
      if ((alleig .or. indeig .and. il==1 .and. iu==n) .and. abstol<=zero) then
         call dcopy (n, work(indd), 1, w, 1)
         indee = indwrk + 2*n
         if (.not.wantz) then
            call dcopy (n - 1, work(inde), 1, work(indee), 1)
            call dsterf (n, w, work(indee), info)
         else
            call dopgtr (uplo, n, ap, work(indtau), z, ldz, work(indwrk), iinfo&
               )
            call dcopy (n - 1, work(inde), 1, work(indee), 1)
            call dsteqr (jobz, n, w, work(indee), z, ldz, work(indwrk), info)
            if (info == 0) then
               ifail(:n) = 0
            endif
         endif
         if (info == 0) then
            m = n
            go to 20
         endif
         info = 0
      endif
!
!     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
!
      if (wantz) then
         order = 'B'
      else
         order = 'E'
      endif
      indibl = 1
      indisp = indibl + n
      indiwo = indisp + n
      call dstebz (range, order, n, vll, vuu, il, iu, abstll, work(indd), work(&
         inde), m, nsplit, w, iwork(indibl), iwork(indisp), work(indwrk), iwork&
         (indiwo), info)
!
      if (wantz) then
         call dstein (n, work(indd), work(inde), m, w, iwork(indibl), iwork(&
            indisp), z, ldz, work(indwrk), iwork(indiwo), ifail, info)
!
!        Apply orthogonal matrix used in reduction to tridiagonal
!        form to eigenvectors returned by DSTEIN.
!
         call dopmtr ('L', uplo, 'N', n, m, ap, work(indtau), z, ldz, work(&
            indwrk), info)
      endif
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
   20 continue
      if (iscale == 1) then
         if (info == 0) then
            imax = m
         else
            imax = info - 1
         endif
         call dscal (imax, one/sigma, w, 1)
      endif
!
!     If eigenvalues are not in order, then sort them, along with
!     eigenvectors.
!
      if (wantz) then
         do j = 1, m - 1
            i = 0
            tmp1 = w(j)
            do jj = j + 1, m
               if (w(jj) >= tmp1) cycle
               i = jj
               tmp1 = w(jj)
            end do
!
            if (i == 0) cycle
            itmp1 = iwork(indibl+i-1)
            w(i) = w(j)
            iwork(indibl+i-1) = iwork(indibl+j-1)
            w(j) = tmp1
            iwork(indibl+j-1) = itmp1
            call dswap (n, z(1,i), 1, z(1,j), 1)
            if (info == 0) cycle
            itmp1 = ifail(i)
            ifail(i) = ifail(j)
            ifail(j) = itmp1
         end do
      endif
!
      return
!
!     End of DSPEVX
!
      end subroutine dspevx
