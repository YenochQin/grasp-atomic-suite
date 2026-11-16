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
      use, intrinsic :: iso_fortran_env, only: real64
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
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jblock
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      integer, intent(out) :: irestart
      real(real64)  :: atwinv
      real(real64)  :: elsto
      real(real64) :: slf_en(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iread, ios, ncfdum, iccutdum, myiddum, nprocsdum, i, &
         ios2, nelc, ir, irowdum, ipos, nrows, j, icstrt
      real(real64) :: stoel, dum, eav0
!-----------------------------------------------
!
      stored_matrix_element_count = 0                      ! Counting continues in setham
      hamiltonian_energy_shift    = 0._real64
      elsto  = 0._real64

! See how much had been done (Hamiltonian matrix)
! irestart is set;
! iread accumulated;
! nelmnt, eav, elsto obtained (to be further modified in setham)

      iread = 0                       ! # of rows read, initialization necessary

      read (imcdf, iostat=ios) ncfdum, iccutdum, myiddum, nprocsdum
      irestart = ios

      if (ios == 0) then

         if (csf_count/=ncfdum .or. iccut(1)/=iccutdum .or. myid/=myiddum .or. &
            nprocs/= nprocsdum) then
            write (stderr_unit, *) csf_count, ncfdum, iccut(1), iccutdum, myid, myiddum, &
               nprocs, nprocsdum, 'check'
            stop 'genmat:1'
         endif

         do i = myid + 1, csf_count, nprocs
            read (imcdf, iostat=ios2) nelc, stoel, (dum,ir=2,nelc), eav0, (&
               irowdum,ir=1,nelc)
                            ! Lower triangle row-mode, diagonal last
            if (ios2 == 0) then
               iread = iread + 1
               stored_matrix_element_count = stored_matrix_element_count + nelc
               hamiltonian_energy_shift = hamiltonian_energy_shift + eav0
               elsto = stoel
            else
               exit
            endif
         end do
         ipos = 7 + orbital_count + orbital_count + iread + 1
      else
         ipos = 7 + orbital_count + orbital_count
      endif

! Find the maximum number of rows

      nrows = (csf_count - myid - 1 + nprocs)/nprocs
      if (csf_count < nprocs) nrows = csf_count/(myid + 1)

! Report the number of rows read.
! A more suitable report on all nodes can be done here, but this will
! set a synchronization point.

      write (stderr_unit, *) iread, ' (total ', nrows, ') rows read from .res'
      if (myid == 0) write (24, *) iread, ' (total ', nrows, &
         ') rows read from .res'

! Position the file for the next record from setham

      do i = 1, jblock - 1
         j = (csf_counts_by_block(i) - myid - 1 + nprocs)/nprocs
         if (csf_counts_by_block(i) < nprocs) j = csf_counts_by_block(i)/(myid + 1)
         ipos = ipos + j + 1
      end do
      call posfile (0, imcdf, ipos)

      if (ios /= 0) write (imcdf) csf_count, iccut(1), myid, nprocs

      if (iread < nrows) then
         icstrt = iread*nprocs + myid + 1
!     ...Generate the rest of the Hamiltonian matrix
         call setham (myid, nprocs, jblock, elsto, icstrt, stored_matrix_element_count, atwinv, &
            slf_en)
      else
         nelmnttmp = stored_matrix_element_count
         ncftmp = csf_count
      endif

      return
      end subroutine genmat
