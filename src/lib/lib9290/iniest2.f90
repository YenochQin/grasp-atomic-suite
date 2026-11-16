!***********************************************************************
!
      subroutine iniest2(nmax, ncf, niv, basis, hmx, jcol, irow)
!
!  Serial version of iniestmpi.
!  Structure of the input sparse matrix hmx:
!    . It's a 1-d array
!    . Length: 1 to jcol(ncf)
!    . Number of non-zero elements for column j is:
!          jcol(j) - jcol(j-1) + 1
!    . Row index for element hmx(i) is irow(i)
!  Xinghong He  98-10-28
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:36:59   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use dinit_i
      use memory_man
      use suite_linear_algebra, only: dcopy, dspevx
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nmax
      integer, intent(in) :: ncf
      integer  :: niv
      integer, dimension(0:*), intent(in) :: jcol
      integer, dimension(*), intent(in) :: irow
      real(double), dimension(*) :: basis
      real(double), dimension(*), intent(in) :: hmx
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ns, joffspar, j, joffnorm, ir, nfound, info, ierr

      integer, dimension(:), pointer :: iwork,ifail
      real(double), dimension(:), pointer :: ap, eigval,vec, work

!-----------------------------------------------

      ns = min(nmax,ncf)
      call alloc (ap, (ns*(ns + 1))/2, 'AP', 'INIEST2')
      call dinit ((ns*(ns + 1))/2, 0.d0, ap, 1)

!  Expand the sparse form to normal form for upper-right sub-matrix

      joffspar = 0                               ! offset for sparse form
      do j = 1, ns
         joffnorm = (j*(j - 1))/2                ! offset for normal form
         do ir = joffspar + 1, jcol(j)
            ap(irow(ir)+joffnorm) = hmx(ir)
         end do
         joffspar = jcol(j)
      end do

!  Merge ap from all nodes and then send to all nodes

      call alloc(eigval, ns, 'EIGVAL', 'INIEST2')
      call alloc (vec, ns*niv, 'VEC', 'INIEST2')
      call alloc (work, 8*ns, 'WORK', 'INIEST2')
      call alloc (iwork, 8*ns, 'IWORK', 'INIEST2')
      call alloc (ifail, ns, 'IFAIL', 'INIEST2')

      call dspevx ('Vectors also', 'In a range', 'Upper triangular', ns, ap, &
         -1.d0, -1.d0, 1, niv, 0.d0, nfound, eigval, vec, ns, work, iwork, ifail, &
         info)
      ierr = -abs(info)

!  Build the Basis.

      call dinit (ncf*niv, 0.d0, basis, 1)

!  scatter the vectors

      do j = 1, niv
         call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)), 1, basis(ncf*(j-1)+1), 1)
      end do

      call dcopy (niv, eigval, 1, basis(niv*ncf+1), 1)

      !deallocate(ap)
      !deallocate(eigval)
      !deallocate(vec)
      !deallocate(work)
      !deallocate(iwork)
      !deallocate(ifail)
       call dalloc (ap, 'AP', 'INIEST2')
       call dalloc (eigval, 'EIGVAL', 'INIEST2')
       call dalloc (vec, 'VEC', 'INIEST2')
       call dalloc (work, 'WORK', 'INIEST2')
       call dalloc (iwork, 'IWORK', 'INIEST2')
       call dalloc (ifail, 'IFAIL', 'INIEST2')

      return
      end subroutine iniest2
