module initdvd_i
    implicit none (type, external)
    private
    public :: initdvd
contains
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
end module initdvd_i
