!***********************************************************************
!
      subroutine lodstate(idblk)
!
!   Print block info and ask ASF serial numbers for each block
!
!   Input:
!      nblock, idblk(), ncfblk()
!
!   Output:
!      nevblk(), ncmaxblk()
!      ncmin, iccmin(1:ncmin) -- via items, in common block
!
!    Memories allocated here but not de-allocated are
!        pccmin (via items, deallocated here only if ncmin=0)
!
!  Written by Xinghong He                                  Jul 17 1997
!  Updated by Xinghong He                                  Jun 10 1998
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use def_c
      use hblock_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use items_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = 8) :: idblk(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, jblock, ncf, ncminold, ierr, ntmp
      character :: str*256
!-----------------------------------------------
!

      write (stderr_unit, *) 'There are ', symmetry_block_count, ' blocks  ', &
         '(block   J/Parity   NCF):'
      write (stderr_unit, '( 4(I3, 1X, A5, I8, 5X) )') (j,idblk(j)(1:5),csf_counts_by_block(j),j=1&
         ,symmetry_block_count)

      eigenvector_counts_by_block = 0
      ncmaxblk = 0

      write (stderr_unit, *)
      write (stderr_unit, *) 'Enter ASF serial numbers for each block'

      ncmin = 0
  123 continue
      do jblock = 1, symmetry_block_count
         ncf = csf_counts_by_block(jblock)
  234    continue
         write (stderr_unit, *) 'Block ', jblock, '   ncf = ', ncf, ' id = ', idblk(&
            jblock)(1:5)

!        ...Read and parse the list of levels

         read (*, '(A)') str
         write(734,'(a)') trim(str) ! write to rscf.log file see, rscf

!        ...ICCMIN is allocated and accumulated in items
!        ...ncmin is both input and output parameters to items
         ncminold = ncmin
         call items (ncmin, ncf, str, ierr)
         if (ncmin == 0) call dalloc (iccmin, 'ICCMIN', 'LODSTATE' )
         if (ierr < 0) go to 234
         eigenvector_counts_by_block(jblock) = ncmin - ncminold

!        ...Determine ncmaxblk
         ntmp = 0
         do i = ncminold + 1, ncmin
            ntmp = max(ntmp,iccmin(i))
         end do
         ncmaxblk(jblock) = ntmp
      end do

      if (ncmin == 0) then
         write (stderr_unit, *) 'At least one state should be selected'
         go to 123
      endif

      return
      end subroutine lodstate
