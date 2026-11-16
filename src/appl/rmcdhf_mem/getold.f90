!***********************************************************************
!                                                                      *
      subroutine getold(idblk)
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:25:01   1/ 6/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw
      use memory_man
      use blkidx_c
      use damp_c, only: cdamp, odamp
      use def_c, rwtdum=>wt
      use default_c
      use fixd_c, only: lfix, nfix
      use hblock_c
      use iounit_c
      use invt_c, only: noinvt
      use orthct_c
      use orb_c
      use orba_c, only: iorder
      use corre_c, only: lcorre
      use scf_c, only: scnsty,method
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodstate_i
      use getoldwt_i
      use prtrsl_i
      use getrsl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum
      integer , dimension(nnnw) :: indx
      integer :: i, nsubs, norder, loc, noffset, jblock, j
      logical :: getyn, yes
      character :: record*256
!-----------------------------------------------

      !WRITE (istde,*) 'EOL type calculation;'

! lodstate fills
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)

      call alloc (ncmaxblk, nblock, 'NCMAXBLK', 'GETOLD')
      call alloc (nevblk, nblock, 'NEVBLK', 'GETOLD' )

!cjb  ncmaxblk & nevblk initialised in lodstate
!     NEVBLK = 0
!     NCMAXBLK = 0
!cjb

!cjb LODSTATE(NBLOCK, NCFBLK, IDBLK, NEVBLK, NCMAXBLK) -> (IDBLK)
!     CALL LODSTATE (NBLOCK, NCFBLK(1), IDBLK, NEVBLK, NCMAXBLK)
      call lodstate (idblk)
!cjb
!
!   Allocate the storage for and set the weights
!
      call alloc (weight, ncmin, 'WEIGHT', 'GETOLD')

      call getoldwt (ndef, ncmin, weight)
!
!   Eigenvector damping
!
      call alloc (cdamp, ncmin, 'CDAMP', 'GETOLD')
!
      cdamp(:ncmin) = 0.d0
!
!   Print the list of all subshells
!
      write (istde, *) 'Radial functions'
      call prtrsl
!
!   Determine which orbitals are to be varied, which are fixed.
!   Quantities determined here:
!     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
!   Instead of broadcasting these quantities, we broadcast
!   the intermediate result from GETRSL (see below)
!
      lfix(:nw) = .true.

      write (istde, *) 'Enter orbitals to be varied (Updating order)'
      call getrsl (indx, nsubs)

      lfix(indx(:nsubs)) = .false.
!XHH      give a big value, rather than zero to scnsty()
      scnsty(indx(:nsubs)) = 1.d20
      nfix = nw - nsubs
      if (nfix == nw) write (istde, *) &
         'All subshell radial wavefunctions are fixed;', &
         ' performing CI calculation.'

!   Determine orbital updating order

      norder = 0
      do i = 1, nw
         iorder(i) = i
         if (lfix(i)) cycle
         norder = norder + 1
         iorder(i) = indx(norder)
      end do
!
!XHH added a array to store the index of the correlation functions
!
      lcorre(:nw) = .true.

      write (istde, *) 'Which of these are spectroscopic orbitals?'
      call getrsl (indx, nsubs)
      if (nsubs > 0) then
         do i = 1, nsubs
            loc = indx(i)
            if (lfix(loc)) cycle
            method(loc) = 1
            noinvt(loc) = .false.
            odamp(loc) = 0.d0
            lcorre(loc) = .false.
         end do
      endif

! Set NSIC. It will be non-zero if all orbitals to be varied are
! spectroscopic orbitals

      nsic = (nw - nfix)/4
      do i = 1, nw
         if (.not.(.not.lfix(i) .and. lcorre(i))) cycle
         nsic = 0
         exit
      end do
!
      nscf = 24
      nsolv = 3
      orthst = .true.
!
!   Make the allocation for the auxiliary vector required
!   by SUBROUTINE NEWCO
!
      call alloc (rwtdum, ncmin, 'RWTDUM', 'GETOLD')
!
!   Place the block numbers of the all ncmin eigenstate(wanted)
!   in array idxblk
!
      call alloc (idxblk, ncmin, 'IDXBLK', 'GETOLD')
      noffset = 0
      do jblock = 1, nblock
         do j = 1, nevblk(jblock)
            idxblk(j+noffset) = jblock
         end do
         noffset = noffset + nevblk(jblock)
      end do
      if (noffset /= ncmin) stop 'getold: ncmin trouble'

      return
      end subroutine getold
