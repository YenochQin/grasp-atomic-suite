module iniestmpi_i
    implicit none (type, external)
    private
    public :: iniestmpi
contains
!***********************************************************************
!
      subroutine iniestmpi (nmax, ncf, niv, basis, hmx, jcol, irow)
!
!       Routine for providing initial estimates from the diagonal
!     of the matrix. This way was used by Dvdson in atomic structure
!     calculations. It should be used to obtain estimates when nothing
!     else is available.
!
!       nmax is typically 1000

!  Structure of the input sparse matrix hmx:
!    . It's a 1-d array
!    . Length: 1 to jcol(ncf)
!    . Number of non-zero elements for column j is:
!          jcol(j) - jcol(j-1) + 1
!    . Row index for element hmx(i) is irow(i)
!  Xinghong He  98-10-28
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use mpi_c

      use suite_linear_algebra, only: dcopy, dspevx
      use dinit_i, only: dinit
      use mpiu_i, only: gdsummpi
      use suite_linear_algebra, only: dlamch
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nmax
      integer, intent(in) :: ncf
      integer  :: niv
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
      real(real64)  :: basis(*)
      real(real64), intent(in) :: hmx(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ns, joffspar, j, joffnorm, ir, nfound, info

      integer, pointer :: iwork(:),ifail(:)
      real(real64), pointer :: ap(:), eigval(:),vec(:), work(:)
      real(real64) :: abstol

!     real dLAMCH
!     external dLAMCH

!-----------------------------------------------

      ns = min(nmax,ncf)
      call alloc (ap, (ns*(ns + 1))/2, "AP", "INIESTmpi")
      call dinit ((ns*(ns + 1))/2, 0.0_real64, ap, 1)

!  Expand the sparse form to normal form for upper-right sub-matrix

      joffspar = 0                               ! offset for sparse form
!     DO J = 1, NS
      do j = mpi_rank + 1, ns, mpi_process_count
         joffnorm = (j*(j - 1))/2                ! offset for normal form
         do ir = joffspar + 1, jcol(j)
            ap(irow(ir)+joffnorm) = hmx(ir)
         end do
         joffspar = jcol(j)
      end do

!  Merge ap from all nodes and then send to all nodes

      call gdsummpi (ap, (ns*(ns+1))/2)

      call alloc(eigval, ns,     "EIGVAL", "INIESTmpi")
      call alloc (vec,   ns*niv, "VEC",    "INIESTmpi")
      call alloc (work,  8*ns,   "WORK",   "INIESTmpi")
!GG      CALL ALLOC (iwork, 8*ns,   'IWORK',  'INIESTmpi')
      call alloc (iwork, 5*ns,   "IWORK",  "INIESTmpi")
      call alloc (ifail, ns,     "IFAIL",  "INIESTmpi")

      abstol = 2*dlamch("S")

      call dspevx ("Vectors also", "In a range", "Upper triangular", ns, ap, &
      -1.0_real64, -1.0_real64, 1, niv, abstol, nfound, eigval, vec, ns, work, iwork, ifail, &
         info)
      mpi_error_code = -abs(info)

!  Build the Basis.

      call dinit (ncf*niv, 0.0_real64, basis, 1)

!  scatter the vectors

      do j = 1, niv
         call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)), 1, basis(ncf*(j-1)+1), 1)
      end do

      call dcopy (niv, eigval, 1, basis(niv*ncf+1), 1)

       call dalloc (ap,     "AP",     "INIESTmpi")
       call dalloc (eigval, "EIGVAL", "INIESTmpi")
       call dalloc (vec,    "VEC",    "INIESTmpi")
       call dalloc (work,   "WORK",   "INIESTmpi")
       call dalloc (iwork,  "IWORK",  "INIESTmpi")
       call dalloc (ifail,  "IFAIL",  "INIESTmpi")

      return
      end subroutine iniestmpi
end module iniestmpi_i
