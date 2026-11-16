!***********************************************************************
!                                                                      *
      subroutine getoldmpi(idblk)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_orbitals
      use memory_man
      use blkidx_c
      use damp_c, only: mixing_damping_factors, orbital_damping_factors
      use def_c, rwtdum=>wt
      use default_c
      use fixd_c, only: is_orbital_fixed, fixed_orbital_count
      use hblock_c
      use iounit_c
      use invt_c, only: noinvt
      use orthct_c
      use orb_c
      use orba_c, only: iorder
      use corre_c, only: lcorre
      use scf_c, only: scnsty,method
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodstate_i
      use getoldwt_i
      use prtrsl_i
      use getrsl_i
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

      call alloc (ncmaxblk, symmetry_block_count, 'NCMAXBLK', 'GETOLDmpi')
      call alloc (eigenvector_counts_by_block, symmetry_block_count, 'NEVBLK', 'GETOLDmpi' )

!cjb  ncmaxblk & nevblk initialised in lodstate
!     NEVBLK = 0
!     NCMAXBLK = 0
!cjb

!cjb LODSTATE(NBLOCK, NCFBLK, IDBLK, NEVBLK, NCMAXBLK) -> (IDBLK)
!cjb   CALL LODSTATE (NBLOCK, NCFBLK(1), IDBLK, NEVBLK, NCMAXBLK)
      if (mpi_rank == 0) then
         call lodstate (idblk)
      end if
!cjb

      call mpi_bcast (ncmin, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast (eigenvector_counts_by_block(1), symmetry_block_count, mpi_integer, 0,            &
                      mpi_comm_world, mpi_error_code)
      call mpi_bcast (ncmaxblk(1), symmetry_block_count, mpi_integer, 0,          &
                      mpi_comm_world, mpi_error_code)

      ! actually, ncmin will always be positive
      if (ncmin /= 0) then
         if (mpi_rank /=  0) then   ! lodstate allocated it from node-0
!GG            CALL alloc (pccmin, ncmin, 'PCCMIN', 'GETOLD')
            call alloc (iccmin, ncmin, 'ICCMIN', 'GETOLDmpi')
         endif

         call mpi_bcast (iccmin(1), ncmin, mpi_integer, 0,         &
                         mpi_comm_world, mpi_error_code)
      endif

!
!   Allocate the storage for and set the weights
!
      call alloc (weight, ncmin, 'WEIGHT', 'GETOLDmpi')

      if (mpi_rank == 0) then
         call getoldwt (ndef, ncmin, weight)
      end if
      call mpi_bcast (weight(1), ncmin, mpi_double_precision, 0,  &
                      mpi_comm_world, mpi_error_code)
!
!   Eigenvector damping
!
      call alloc (mixing_damping_factors, ncmin, 'CDAMP', 'GETOLDmpi')
!
      mixing_damping_factors(:ncmin) = 0._real64
!
!   Print the list of all subshells
!
      if (mpi_rank == 0) then
         write (stderr_unit, *) 'Radial functions'
         call prtrsl
      end if
!
!   Determine which orbitals are to be varied, which are fixed.
!   Quantities determined here:
!     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
!   Instead of broadcasting these quantities, we broadcast
!   the intermediate result from GETRSL (see below)
!
      is_orbital_fixed(:orbital_count) = .true.

      if (mpi_rank == 0) then
         write (stderr_unit, *) 'Enter orbitals to be varied (Updating order)'
         call getrsl (indx, nsubs)
      end if
      call mpi_bcast (nsubs, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      if (nsubs > 0) call mpi_bcast (indx(1), nsubs, mpi_integer, 0,  &
                                     mpi_comm_world, mpi_error_code)

      is_orbital_fixed(indx(:nsubs)) = .false.
!XHH      give a big value, rather than zero to scnsty()
      scnsty(indx(:nsubs)) = 1.e20_real64
      fixed_orbital_count = orbital_count - nsubs
      if (fixed_orbital_count == orbital_count) then
         if(mpi_rank == 0)                                                 &
         write (stderr_unit,*)'All subshell radial wavefunctions are fixed;',&
                        ' performing CI calculations with RCI.'
      end if

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

      if (mpi_rank == 0) then
         write (stderr_unit, *) 'Which of these are spectroscopic orbitals?'
         call getrsl (indx, nsubs)
      end if
      call mpi_bcast (nsubs, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)
      if (nsubs .gt. 0) then
         call mpi_bcast (indx(1), nsubs, mpi_integer, 0,            &
                         mpi_comm_world, mpi_error_code)
!GG      IF (NSUBS > 0) THEN
         do i = 1, nsubs
            loc = indx(i)
            if (is_orbital_fixed(loc)) cycle
            method(loc) = 1
            noinvt(loc) = .false.
            orbital_damping_factors(loc) = 0._real64
            lcorre(loc) = .false.
         end do
      endif

! Set NSIC. It will be non-zero if all orbitals to be varied are
! spectroscopic orbitals

      nsic = 4 + (orbital_count - fixed_orbital_count)/4
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
      call alloc (rwtdum, ncmin, 'RWTDUM', 'GETOLDmpi')
!
!   Place the block numbers of the all ncmin eigenstate(wanted)
!   in array idxblk
!
      call alloc (idxblk, ncmin, 'IDXBLK', 'GETOLDmpi')
      noffset = 0
      do jblock = 1, symmetry_block_count
         do j = 1, eigenvector_counts_by_block(jblock)
            idxblk(j+noffset) = jblock
         end do
         noffset = noffset + eigenvector_counts_by_block(jblock)
      end do
      if (noffset /= ncmin) error stop 'getold: ncmin trouble'

      return
      end subroutine getoldmpi
