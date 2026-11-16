!***********************************************************************
!                                                                      *
      subroutine scf(eol, rwffile2)
      use, intrinsic :: iso_fortran_env, only: real32
!                                                                      *
!   This  subroutine  performs  the SCF iterations. The procedure is   *
!   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
!   (1986) 290.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!               [RSCF92]: improv, matrix, MAXARR, newco,               *
!                         ORBOUT, ORTHSC, setlag.                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   Block version by Xinghong He          Last revision: 05 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use blkidx_c
      use default_c
      use def_c
      use debug_c
      use eigv_c
      use fixd_c
      use hblock_c
      use iounit_c
      use lagr_c
      use mcpa_c
      use mpi_s
      use pos_c
      use peav_c
      use orb_c
      use orba_c
      use scf_c
      use orthct_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use matrix_i
      use newco_i
      use setlag_i
      use improv_i
      use maxarr_i
      use prwf_i
      use orthsc_i
      use orbout_i
      use endsum_i
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) :: eol
      character  :: rwffile2*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, nit, jseq, kount, k
      real(real64) :: wtaev, wtaev0, dampmx
      logical :: convg, lsort, dvdfirst
!-----------------------------------------------


      total_csf_count = csf_count
      !IF (myid .EQ. 0) PRINT *, '===SCF==='

!=======================================================================
!   Determine Orthonomalization order --- lsort
!=======================================================================

      if (ndef == 0) then
         lsort = .false.
      else
  123    continue
         write (stderr_unit, *) 'Orthonomalization order? '
         write (stderr_unit, *) '     1 -- Update order'
         write (stderr_unit, *) '     2 -- Self consistency connected'
         read (stdin_unit, *) j
         if (j == 1) then
            lsort = .false.
         else if (j == 2) then
            lsort = .true.
         else
            write (stderr_unit, *) 'Input is wrong, redo...'
            go to 123
         endif
      endif

!=======================================================================
!   Deallocate storage that will no longer be used
!=======================================================================

!GG      CALL DALLOC (JQSA, 'JQSA', 'SCF')

!=======================================================================
!   Allocate and fill in auxiliary arrays
!=======================================================================

      call alloc (csf_offsets_by_block, symmetry_block_count, 'NCFPAST', 'SCF')
      call alloc (eigenvalue_offsets_by_block, symmetry_block_count, 'NCMINPAST', 'SCF')
      call alloc (mixing_offsets_by_block, symmetry_block_count, 'NEVECPAST', 'SCF')
      call alloc (eavblk, symmetry_block_count, 'EAVBLK', 'SCF')

      csf_offsets_by_block(1) = 0
      eigenvalue_offsets_by_block(1) = 0
      mixing_offsets_by_block(1) = 0
      do i = 2, symmetry_block_count
         csf_offsets_by_block(i) = csf_offsets_by_block(i-1) + csf_counts_by_block(i - 1)
         eigenvalue_offsets_by_block(i) = eigenvalue_offsets_by_block(i-1) + eigenvector_counts_by_block(i - 1)
         mixing_offsets_by_block(i) = mixing_offsets_by_block(i-1) + eigenvector_counts_by_block(i - 1)*csf_counts_by_block(i - 1)
      end do

      !*** Size of the eigenvector array for all blocks
      mixing_coefficient_storage_size = mixing_offsets_by_block(symmetry_block_count) + &
         & eigenvector_counts_by_block(symmetry_block_count)*csf_counts_by_block(symmetry_block_count)

      if (eol) then
         call alloc (shifted_eigenvalues, ncmin, 'EVAL', 'SCF')
         call alloc (csf_mixing_coefficients, mixing_coefficient_storage_size, 'EVEC', 'SCF')
!GG         CALL ALLOC (IATJPO, NCMIN, 'IATJPO', 'SCF')
!GG         CALL ALLOC (IASPAR, NCMIN, 'IASPAR', 'SCF')
      endif

!=======================================================================
!
!=======================================================================
      nddim = 0
      nxdim = 0
      nydim = 0

!     This call should only be made AFTER the call to newco
!     CALL setlag (EOL)

!   For (E)OL calculations, determine the level energies and
!   mixing coefficients

!CFF   .. set the logical variable dvdfirst
      dvdfirst = .true.
      if (eol) then
         call matrix (dvdfirst)
         call newco (wtaev)
      endif
      wtaev0 = 0.0
      dvdfirst = .false.
      do nit = 1, nscf
         if (mpi_rank == 0) write (*, 301) nit

!   For all pairs constrained through a Lagrange multiplier, compute
!   the Lagrange multiplier

         call setlag (eol)

!   Improve all orbitals in turn

         dampmx = 0.0
         if (mpi_rank == 0) write (*, 302)
         do j = 1, orbital_count
            jseq = iorder(j)
            if (is_orbital_fixed(jseq)) cycle
            call improv (eol, jseq, lsort, dampmx)
         end do
!
!   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
!   improve it
!
!         write(istde,*) 'nsic=',nsic
         do kount = 1, nsic
            call maxarr (k)
            if (k == 0) then
               convg = .true.
               go to 3
            else
               if (scnsty(k) <= accy) then
                  convg = .true.
                  go to 3
               endif
            endif
            call improv (eol, k, lsort, dampmx)
         end do

         call maxarr (k)

         if (k == 0) then
            convg = .true.
         else
            if (scnsty(k) <= accy) then
               convg = .true.
            else
               convg = .false.
            endif
         endif

    3    continue
         if (ldbpr(24) .and. mpi_rank==0) call prwf (0)

!   Perform Gram-Schmidt process
!   For OL calculation, orthst is true and orbitals are orthonormalized
!   in subroutine improv. For AL calculation, orthst is false.
         if (.not.orthst) call orthsc

!   Write the subshell radial wavefunctions to the .rwf file

         if (mpi_rank == 0) call orbout (rwffile2)

         if (eol) then
            call matrix(dvdfirst)
            call newco (wtaev)
         endif
!        Make this a relative convergence test
!        IF(ABS(WTAEV-WTAEV0).LT.1.0D-9.and.
!    &   DAMPMX.LT.1.0D-4) CONVG=.true.
!        PRINT *, 'WTAEV, WTAEV0', WTAEV, WTAEV0,
!    &             ABS((WTAEV-WTAEV0)/WTAEV)
!cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
!cjb     IF(DABS(WTAEV-WTAEV0).LT.1.0D-8.and.              &
!cjb                   DAMPMX.LT.1.0D-2) CONVG=.true.
         if (abs((wtaev - wtaev0)/wtaev) < 0.001_real32*accy) convg = .true.
         wtaev0 = wtaev
         if (.not.convg) cycle
         if (ldbpr(25) .and. .not.ldbpr(24) .and. mpi_rank==0) call prwf (0)
            !IF (EOL) CALL matrix (dvdfirst)
         go to 5

      end do

      if (mpi_rank == 0) write (stderr_unit, *) ' Maximum iterations in SCF Exceeded.'

    5 continue
      do i = 31, 32 + max_mcp_file_multipole_rank
         close(i)                                ! The MCP coefficient files
      end do

      if (mpi_rank == 0) then
         !CLOSE (23)     ! The .rwf file
         close(25)                               ! The .mix file
!
!   Complete the summary - moved from rscf92 for easier alloc/dalloc
!
         call endsum
      endif
!
!   Deallocate storage
!
      call dalloc (wt, 'WT', 'SCF')                       !Either getold or getald

      if (nec > 0) then
         call dalloc (iecc, 'IECC', 'SCF')
         call dalloc (ecv, 'ECV', 'SCF')
!GG         CALL DALLOC (IQAR, 'IQAR', 'SCF')
         call dalloc (csf_subshell_occupations, 'IQA', 'SCF')
      endif

      if (nddim > 0) then
         call dalloc (da, 'DA', 'SCF')
         call dalloc (nda, 'NDA', 'SCF')
      endif

      if (nxdim > 0) then
         call dalloc (xa, 'XA', 'SCF')
         call dalloc (nxa, 'NXA', 'SCF')
      endif

      if (nydim > 0) then
         call dalloc (ya, 'YA', 'SCF')
         call dalloc (nya, 'NYA', 'SCF')
      endif

      if (eol) then
         call dalloc (shifted_eigenvalues, 'EVAL', 'SCF')
         call dalloc (csf_mixing_coefficients, 'EvEC', 'SCF')
!GG         CALL DALLOC (IATJPO, 'IATJPO', 'SCF')
!GG         CALL DALLOC (IASPAR, 'IASPAR', 'SCF')
         call dalloc (ncmaxblk, 'NCMAXBLK', 'SCF') ! getold.f
         call dalloc (eavblk, 'EAVBLK', 'SCF')     ! getold.f
         call dalloc (idxblk, 'IDXBLK', 'SCF')     ! Allocated in getold.f
         call dalloc (iccmin, 'ICCMIN', 'SCF')       ! Allocated in items.f<-getold.f
      endif
!
      call dalloc (csf_offsets_by_block, 'NCFPAST', 'SCF')
      call dalloc (eigenvalue_offsets_by_block, 'NCMINPAST', 'SCF')
      call dalloc (mixing_offsets_by_block, 'NEVECPAST', 'SCF')

  301 format(/,' Iteration number ',1i3,/,' --------------------')
  302 format(41x,'Self-            Damping'/,&
         'Subshell    Energy    Method   P0    ',&
         'consistency  Norm-1  factor  JP',' MTP INV NNP'/)

      return
      end subroutine scf
