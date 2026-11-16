!***********************************************************************
!
      subroutine cslh(name, ncore, nblkin, idblk)
!  A container which calls setcsll to open, read <name> file to get
!     nblock, ncfblk(), idblk(), ncftot.
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Called by rcivu, mcpvu
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:12   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use hblock_c
      use def_c
      use orb_c, ncftot=>ncf
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcsll_i
      use lodcsh_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      integer  :: nblkin
      character  :: name*(*)
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum
      integer :: ierr
!-----------------------------------------------
!

! node-0 does exactly the same as the serial code does

      call alloc (ncfblk, nblkin + 1, 'NCFBLK', 'CSLH')

      call setcsll (21, name, nblkin, nblock, ncfblk, ncftot, idblk)
      call ralloc (ncfblk, nblock + 1, 'NCFBLK', 'CSLH')
      rewind (21)
      read (21, *)
      !..Load header of <name> file
      call lodcsh (21, ncore)

      return
      end subroutine cslh
