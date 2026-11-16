module adds_i
    implicit none (type, external)
    private
    public :: adds
contains
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
end module adds_i
