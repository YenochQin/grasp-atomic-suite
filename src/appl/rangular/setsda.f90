!***********************************************************************
!                                                                      *
      subroutine setsda(outsdadummy,nnonz,lprint,nb,myid,nprocs,fhead)
!                                                                      *
!   This routine examines lists                                        *
!                               (IC,IR,npos)                           *
!   to set up the array  IENDC  required by the Davidson eigensolver   *
!   of Stathopoulos and Fischer.                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
!   Modified for block interactions by C. F. Fischer        May 1997   *
!   Modified for multi-processors   by Xinghong He       03 Jul 1998   *
!
!  Currently shared by mcpblk, mcpmpi
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:52:18   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use orb_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use outsda_i, only: outsda
      implicit none
      procedure(outsda) :: outsdadummy
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nnonz
      integer, intent(in) :: nb
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      logical  :: lprint
      character, intent(in) :: fhead*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mb, iend, iclast, i, ic, npos, ierr
      integer, dimension(:), pointer :: iendc
      integer, dimension(:), pointer :: irow
      character :: mcplab*3
!-----------------------------------------------

      if (myid == 0) write (6, *) &
         'Analysing sparse matrix array definition file ...', 30

      read (30) mcplab, mb
      if (nb /= mb) then
         write (istde, *) 'setsda: nb = ', nb, '.NE. mb (=', mb, ')'
         stop
      endif
!
!   Allocate storage for IENDC(0:NCF)
!
      call alloc (iendc, 0,  ncf,'IENDC','SETSDA'  )
      call alloc (irow, nnonz, 'IROW', 'SETSDA' )
!
!   Analyse data on file 30; set up IENDC and IROW
! In multiprocessor environment, iendc of each node will have the
! same length (ncf+1); but will have its own part filled. irow is
! local, and its length is determined by the local parameter nnonz.

      iend = 0
      iclast = 0
      do i = 1, nnonz
         read (30) ic, irow(i), npos
         if (ic /= iclast) then
            iendc(iclast) = iend
            iclast = ic
         endif
         iend = npos
      end do
!xhh - changed to suits MPI environment as well
!      IENDC(NCF) = IEND
      iendc(ic) = iend
!
!   Sorting complete; rewrite to mcpXXX.30 file
!
      open(29, file=fhead//'.30', status='OLD', form='UNFORMATTED', iostat=ierr&
         , position='APPEND')
      if (ierr /= 0) then
         write (istde, *) ' Error when opening the file mcp.30'
         stop
      endif

      write (29) 'MCP', nb, ncf
      write (29) nnonz
      write (29) (iendc(i),i=myid + 1,ncf,nprocs), (irow(i),i=1,nnonz)
      close(29)
!cjb
!
!   Deallocate storage
!
      call dalloc (iendc, 'IENDC', 'SETSDA')
      call dalloc (irow, 'IROW', 'SETSDA')
!
      call outsdadummy (lprint, nnonz, ncf)
      return
      end subroutine setsda
