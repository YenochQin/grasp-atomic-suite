module iniest_i
    implicit none (type, external)
    private
    public :: iniest
contains
      subroutine iniest(n, nb, niv, hmx, jcol, irow, basis, iblock, jblock)
!-----------------------------------------------------------------------
!     Routine for providing initial estimates from the diagonal
!     of the matrix. This way was used by Dvdson in atomic structure
!       calculations. It should be used to obtain estimates when nothing
!       else is available.
!
!-----------------------------------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dinit_i
      !USE dspevx_I
      !USE vec_I
      !USE dcopy_I
      use suite_linear_algebra, only: dcopy, dspevx
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: nb
      integer  :: niv
      integer, intent(in) :: jblock
      integer, intent(in) :: jcol(0:*)
      integer, intent(in) :: irow(*)
      integer, intent(in) :: iblock(*)
      real(real64), intent(in) :: hmx(*)
      real(real64)  :: basis(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer :: ns, icount, &
         j, istart, ish, ii, nfound, info, ierr, ic
      !REAL(DOUBLE) :: EIGVAL, WORK

      integer, pointer :: iqiwork(:),iqif(:),iqjsh(:)
      integer, allocatable, target :: iwork(:), ifail(:), jsh(:)
      real(real64), pointer :: iqap(:),iqeig(:),iqvec(:),iqwork(:)
      real(real64), allocatable, target :: ap(:), eigval(:), vec(:), work(:)

!-----------------------------------------------
!     !pointer (iqap,ap(1)),(iqeig,eigval(1)),(iqvec,vec(1))
!     !pointer (iqwork,work(1)),(iqiwork,iwork(1)),(iqif,IFAIL(1))
!     !pointer (iqjsh,jsh(1))

!******************************************************************

!***** alloc space for 100*100 lower triangular
      ns = min(1000,nb)
!once here        NS = min(800, NB)
        !print*, 'NS=',NS,' NIV=',NIV
      !CALL ALLOC (IQAP, NS*(NS + 1)/2, 8)
      allocate (ap(ns*(ns + 1)/2),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns*(ns + 1)/2,'iniest::ap',ierr);
      iqap=>ap

      call dinit (ns*(ns + 1)/2, 0._real64, ap, 1)

      !CALL ALLOC (IQEIG, NS, 8)
      allocate (eigval(ns),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns,'iniest::eigval',ierr);
      iqeig=>eigval

      !CALL ALLOC (IQVEC, NS*NIV, 8)
      allocate (vec(ns*niv),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns*niv,'iniest::vec',ierr);
      iqvec=>vec

      !CALL ALLOC (IQWORK, 8*NS, 8)
      allocate (work(8*ns),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,8*ns,'iniest::work',ierr);
      iqwork=>work

      !CALL ALLOC (IQIWORK, 5*NS, 4)
      allocate (iwork(5*ns),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,5*ns,'iniest::iwork',ierr);
      iqiwork=>iwork

      !CALL ALLOC (IQIF, NS, 4)
      allocate (ifail(ns),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns,'iniest::IFAIL',ierr);
      iqif=>ifail

      !CALL ALLOC (IQJSH, NS, 4)
      allocate (jsh(ns),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns,'iniest::jsh',ierr);
      iqjsh=>jsh

      icount = 0

!**** separate upper left block of size NS*NS
      do j = 1, n
         if (icount >= ns) exit
         if (iblock(j) /= jblock) cycle
         icount = icount + 1
         jsh(icount) = j
         istart = jcol(j-1) + 1
         ish = j - icount
  100    continue
         ap(irow(istart)-ish+(icount-1)*(2*ns-icount)/2) = hmx(istart)
!           print*, 'ic = ',icount, ' i=', irow(istart)-ish
         istart = istart + 1
! check the block structure for zero elements
         if (istart > jcol(j)) go to 102
         ish = ish + count(iblock(irow(istart-1):irow(istart))/=jblock)
         if (istart > jcol(j)) go to 102
         if (irow(istart) - ish <= ns) go to 100
!           goto 100
  102    continue
      end do
      call dspevx ('Vectors also', 'In a range', 'Lower triangular', ns, ap, &
         -1._real64, -1._real64, 1, niv, 0._real64, nfound, eigval, vec, ns, work, iwork, ifail, &
         info)
      ierr = -abs(info)
      if (ierr /= 0) write (6, *) 'iniest ierr =', ierr
!           print '(D14.7,X,I2)', (eigval(i),i, i=1,NIV)

!******************************************************************


!
!       ..Build the Basis.
!
      call dinit (n*niv, 0._real64, basis, 1)
      istart = 0
      do j = 1, niv
!scatter the vectors
         do ic = 1, ns
            basis(istart+jsh(ic)) = vec((j - 1)*ns + ic)
         end do
!              call dcopy(NS,vec((J-1)*NS+1),1,BASIS(istart+ish(j)),1)
         istart = istart + n
      end do
      call dcopy (niv, eigval, 1, basis(niv*n+1), 1)
      deallocate(ap) !CALL DALLOC (IQAP)
      deallocate(eigval) !CALL DALLOC (IQEIG)
      deallocate(vec) !CALL DALLOC (IQVEC)
      deallocate(work) !CALL DALLOC (IQWORK)
      deallocate(iwork) !CALL DALLOC (IQIWORK)
      deallocate(ifail) !CALL DALLOC (IQIF)
      deallocate(jsh) !CALL DALLOC (IQJSH)

      return
      contains
      subroutine mem_fail(unit, count, variable, status)
!     The legacy external MEM_FAIL had no implementation in this suite.
      integer, intent(in) :: unit, count, status
      character(*), intent(in) :: variable
      write(unit, '(A,A,A,I0,A,I0)') 'Unable to allocate ', variable, &
         ': elements=', count, ', status=', status
      error stop 'INIEST allocation failed'
      end subroutine mem_fail
      end subroutine iniest
end module iniest_i
