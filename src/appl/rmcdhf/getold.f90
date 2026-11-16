!***********************************************************************
!                                                                      *
      subroutine getold(idblk)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_orbitals
      use memory_man
      use blkidx_c, only: double, idxblk
      use damp_c, only: mixing_damping_factors, orbital_damping_factors
      use def_c, rwtdum=>wt
      use default_c, only: double, inpci, ndef, ndump
      use fixd_c, only: is_orbital_fixed, fixed_orbital_count
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use invt_c, only: noinvt
      use orthct_c, only: orthst
      use orb_c
      use orba_c, only: iorder
      use corre_c, only: lcorre
      use scf_c, only: scnsty,method
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodstate_i, only: lodstate
      use getoldwt_i, only: getoldwt
      use prtrsl_i, only: prtrsl
      use getrsl_i, only: getrsl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqadum
      integer  :: indx(max_orbitals)
      integer :: i, nsubs, norder, loc, noffset, jblock, j
      logical :: getyn, yes
      character :: record*256
!-----------------------------------------------

      !WRITE (istde,*) 'EOL type calculation;'

! lodstate fills
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)

      call alloc (ncmaxblk, symmetry_block_count, "NCMAXBLK", "GETOLD")
      call alloc (eigenvector_counts_by_block, symmetry_block_count, "NEVBLK", "GETOLD" )

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
      call alloc (weight, ncmin, "WEIGHT", "GETOLD")

      call getoldwt (ndef, ncmin, weight)
!
!   Eigenvector damping
!
      call alloc (mixing_damping_factors, ncmin, "CDAMP", "GETOLD")
!
      mixing_damping_factors(:ncmin) = 0.0_real64
!
!   Print the list of all subshells
!
      write (stderr_unit, *) "Radial functions"
      call prtrsl
!
!   Determine which orbitals are to be varied, which are fixed.
!   Quantities determined here:
!     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
!   Instead of broadcasting these quantities, we broadcast
!   the intermediate result from GETRSL (see below)
!
      is_orbital_fixed(:orbital_count) = .true.

      write (stderr_unit, *) "Enter orbitals to be varied (Updating order)"
      call getrsl (indx, nsubs)

      is_orbital_fixed(indx(:nsubs)) = .false.
!XHH      give a big value, rather than zero to scnsty()
      scnsty(indx(:nsubs)) = 1.0e20_real64
      fixed_orbital_count = orbital_count - nsubs
      if (fixed_orbital_count == orbital_count) write (stderr_unit, *) &
         "All subshell radial wavefunctions are fixed;", &
         " performing CI calculation."

!   Determine orbital updating order

      norder = 0
      do i = 1, orbital_count
         iorder(i) = i
         if (is_orbital_fixed(i)) cycle
         norder = norder + 1
         iorder(i) = indx(norder)
      end do
!
!XHH added a array to store the index of the correlation functions
!
      lcorre(:orbital_count) = .true.

      write (stderr_unit, *) "Which of these are spectroscopic orbitals?"
      call getrsl (indx, nsubs)
      if (nsubs > 0) then
         do i = 1, nsubs
            loc = indx(i)
            if (is_orbital_fixed(loc)) cycle
            method(loc) = 1
            noinvt(loc) = .false.
            orbital_damping_factors(loc) = 0.0_real64
            lcorre(loc) = .false.
         end do
      end if

! Set NSIC. It will be non-zero if all orbitals to be varied are
! spectroscopic orbitals

      nsic = (orbital_count - fixed_orbital_count)/4
      do i = 1, orbital_count
         if (.not.(.not.is_orbital_fixed(i) .and. lcorre(i))) cycle
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
      call alloc (rwtdum, ncmin, "RWTDUM", "GETOLD")
!
!   Place the block numbers of the all ncmin eigenstate(wanted)
!   in array idxblk
!
      call alloc (idxblk, ncmin, "IDXBLK", "GETOLD")
      noffset = 0
      do jblock = 1, symmetry_block_count
         do j = 1, eigenvector_counts_by_block(jblock)
            idxblk(j+noffset) = jblock
         end do
         noffset = noffset + eigenvector_counts_by_block(jblock)
      end do
      if (noffset /= ncmin) error stop "getold: ncmin trouble"

      return
      end subroutine getold
