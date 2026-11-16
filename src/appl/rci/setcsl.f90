!***********************************************************************
!
      subroutine setcsl(name, ncore, nblkin, idblk)
!
!  A container which calls setcsll to open, read <name>.c file to get
!     nblock, ncfblk(), idblk(), ncf (it is ncftot here).
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use hblock_c
      use def_c
      use orb_c, only: ncftot=>csf_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE setcsll_I
      use lodcsh_i
      use setcsll_i, only: setcsll
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncore
      integer  :: nblkin
      character  :: name*(*)
      character(len=8), dimension(*) :: idblk
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum
!-----------------------------------------------
!
!     POINTER (ncfblk, ncfblk(0:*))
!-----------------------------------------------------------------------

      call alloc (csf_counts_by_block, 0,  nblkin, 'NCFBLK', 'SETCSL' )
      call setcsll (21, name, nblkin, symmetry_block_count, csf_counts_by_block(1:), ncftot, idblk)
      call ralloc (csf_counts_by_block, 0, symmetry_block_count, 'NCFBLK', 'SETCSL' )

      rewind (21)
      read (21, *)

      !..Load header of <name>.c file
      call lodcsh (21, ncore)

      return
      end subroutine setcsl
