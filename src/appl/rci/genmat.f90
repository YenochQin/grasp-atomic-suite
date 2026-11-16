!***********************************************************************
!                                                                      *
      subroutine genmat(atwinv, jblock, myid, nprocs, elsto, irestart, slf_en)
!
!        Generate Hamiltonian matrix for all blocks
!        This routine calls setham to do the computation. It makes
!        sure that the hamiltonian matrix is complete.
!  Only node-0 has the correct, non-zero elsto. And in any case, elsto
!  is not added to the matrix elements in files *.res
!  This routine has been re-written to work in both single- and
!  multi- processors.
!
! Xinghong He 1998-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use eigv_c
      use iccu_c
      use orb_c
      use where_c
      use hblock_c
      use hmat_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use posfile_i
      use setham_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: jblock
      integer  :: myid
      integer  :: nprocs
      integer, intent(out) :: irestart
      real(double)  :: atwinv
      real(double)  :: elsto
      real(double), dimension(*) :: slf_en
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iread, ios, ncfdum, iccutdum, myiddum, nprocsdum, i, &
         ios2, nelc, ir, irowdum, ipos, nrows, j, icstrt
      real(double) :: stoel, dum, eav0
!-----------------------------------------------
!
      nelmnt = 0                      ! Counting continues in setham
      eav    = 0.d0
      elsto  = 0.d0

! See how much had been done (Hamiltonian matrix)
! irestart is set;
! iread accumulated;
! nelmnt, eav, elsto obtained (to be further modified in setham)

      iread = 0                       ! # of rows read, initialization necessary

      read (imcdf, iostat=ios) ncfdum, iccutdum, myiddum, nprocsdum
      irestart = ios

      if (ios == 0) then

         if (ncf/=ncfdum .or. iccut(1)/=iccutdum .or. myid/=myiddum .or. &
            nprocs/= nprocsdum) then
            write (istde, *) ncf, ncfdum, iccut(1), iccutdum, myid, myiddum, &
               nprocs, nprocsdum, 'check'
            stop 'genmat:1'
         endif

         do i = myid + 1, ncf, nprocs
            read (imcdf, iostat=ios2) nelc, stoel, (dum,ir=2,nelc), eav0, (&
               irowdum,ir=1,nelc)
                            ! Lower triangle row-mode, diagonal last
            if (ios2 == 0) then
               iread = iread + 1
               nelmnt = nelmnt + nelc
               eav = eav + eav0
               elsto = stoel
            else
               exit
            endif
         end do
         ipos = 7 + nw + nw + iread + 1
      else
         ipos = 7 + nw + nw
      endif

! Find the maximum number of rows

      nrows = (ncf - myid - 1 + nprocs)/nprocs
      if (ncf < nprocs) nrows = ncf/(myid + 1)

! Report the number of rows read.
! A more suitable report on all nodes can be done here, but this will
! set a synchronization point.

      write (istde, *) iread, ' (total ', nrows, ') rows read from .res'
      if (myid == 0) write (24, *) iread, ' (total ', nrows, &
         ') rows read from .res'

! Position the file for the next record from setham

      do i = 1, jblock - 1
         j = (ncfblk(i) - myid - 1 + nprocs)/nprocs
         if (ncfblk(i) < nprocs) j = ncfblk(i)/(myid + 1)
         ipos = ipos + j + 1
      end do
      call posfile (0, imcdf, ipos)

      if (ios /= 0) write (imcdf) ncf, iccut(1), myid, nprocs

      if (iread < nrows) then
         icstrt = iread*nprocs + myid + 1
!     ...Generate the rest of the Hamiltonian matrix
         call setham (myid, nprocs, jblock, elsto, icstrt, nelmnt, atwinv, &
            slf_en)
      else
         nelmnttmp = nelmnt
         ncftmp = ncf
      endif

      return
      end subroutine genmat
