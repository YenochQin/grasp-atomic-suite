module ovflow_i
    implicit none (type, external)
    private
    public :: ovflow
contains
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use multbc_i, only: multbc
      use dinit_i, only: dinit
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
      call dinit ((nume*(nume + 1))/2, 0.0_real64, s, 1)
      call dinit (nume*nume, 0.0_real64, svec, 1)
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
end module ovflow_i
