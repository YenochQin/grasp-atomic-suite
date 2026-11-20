!***********************************************************************
!                                                                      *
      subroutine scfmpi(extended_optimal_level, output_wavefunction_path)
!                                                                      *
!   This  subroutine  performs  the SCF iterations. The procedure is   *
!   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
!   (1986) 290.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!               [RSCF92]: improvmpi, matrixmpi, MAXARR, newcompi,      *
!                         ORBOUT, ORTHSC, setlagMpi.                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
!   ifort -i8 version by Alexander Kramida (AK) Last rev. 22 Mar 2016  *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use array_memory_management, only: alloc, dalloc
      use transition_block_indices, only: idxblk
      use calculation_input_defaults, only: inpci, ndef, ndump
      use atomic_calculation_state, only: accy, iccmin, ncmin, nscf, nsic, wt
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use orbital_lagrange_multipliers, only: ecv, iecc, nec
      use mcp_file_metadata, only: max_mcp_file_multipole_rank
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use block_mean_energies, only: eavblk
      use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count
      use scf_orbital_order, only: iorder, max_orbitals
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use scf_orthogonality_control, only: orthst
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use matrixmpi_i, only: matrixmpi
      use newcompi_i, only: newcompi
      use setlagmpi_i, only: setlagmpi
      use improvmpi_i, only: improvmpi
      use maxarr_i, only: maxarr
      use prwf_i, only: prwf
      use radial_orbital_orthogonalization, only: orthsc
      use orbout_i, only: orbout
      use endsum_i, only: endsum
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) :: extended_optimal_level
      character  :: output_wavefunction_path*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: orbital_index, i, nit, jseq, kount, k, l_mpi
      real(real64) :: wtaev, wtaev0, dampmx
      logical :: convg, lsort, dvdfirst
!-----------------------------------------------
!CFF   .. set the logical variable dvdfirst
      dvdfirst = .true.
!-----------------------------------------------------------------------
!AK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
!      ISIZE = sizeof(NCF)
      l_mpi = mpi_logical
!      if (ISIZE.EQ.8) L_MPI = MPI_LOGICAL8

      total_csf_count = csf_count
      !IF (myid .EQ. 0) PRINT *, '===SCF==='

!=======================================================================
!   Determine Orthonomalization order --- lsort
!=======================================================================

      if (ndef == 0) then
         lsort = .false.
      else
         if (mpi_rank == 0) then
  123       continue
            write (stderr_unit, *) "Orthonomalization order? "
            write (stderr_unit, *) "     1 -- Update order"
            write (stderr_unit, *) "     2 -- Self consistency connected"
            read (stdin_unit, *) orbital_index
            if (orbital_index == 1) then
               lsort = .false.
            else if (orbital_index == 2) then
               lsort = .true.
            else
               write (stderr_unit, *) "Input is wrong, redo..."
               goto 123
            end if
         end if
         call mpi_bcast (lsort, 1, l_mpi, 0, mpi_comm_world, mpi_error_code)
      end if

!=======================================================================
!   Deallocate storage that will no longer be used
!=======================================================================

!GG      CALL DALLOC (JQSA, 'JQSA', 'SCFmpi')

!=======================================================================
!   Allocate and fill in auxiliary arrays
!=======================================================================

      call alloc (csf_offsets_by_block, symmetry_block_count, "NCFPAST", "SCFmpi")
      call alloc (eigenvalue_offsets_by_block, symmetry_block_count, "NCMINPAST", "SCFmpi")
      call alloc (mixing_offsets_by_block, symmetry_block_count, "NEVECPAST", "SCFmpi")
      call alloc (eavblk, symmetry_block_count, "EAVBLK", "SCFmpi")

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

      if (extended_optimal_level) then
         call alloc (shifted_eigenvalues, ncmin, "EVAL", "SCFmpi")
         call alloc (csf_mixing_coefficients, mixing_coefficient_storage_size, "EVEC", "SCFmpi")
         call alloc (iatjpo, ncmin, "IATJPO", "SCFmpi")
         call alloc (iaspar, ncmin, "IASPAR", "SCFmpi")
      end if

!=======================================================================
!
!=======================================================================
      nddim = 64
      nxdim = 64
      nydim = 64
!      NDDIM = 0
!      NYDIM = 0
!      NYDIM = 0
!AK    Must initialize DA,XA,YA,... here since MPI routines rely on their existence in all processes
      call alloc (nda, nddim, "NDA","SCFmpi")
      call alloc (da, nddim, "DA","SCFmpi")
      call alloc (nxa, nxdim, "NXA","SCFmpi")
      call alloc (xa, nxdim, "XA","SCFmpi")
      call alloc (nya, nydim, "NYA","SCFmpi")
      call alloc (ya, nydim, "YA","SCFmpi")

!     This call should only be made AFTER the call to newco
!     CALL setlag (EOL)

!   For (E)OL calculations, determine the level energies and
!   mixing coefficients

!CFF   .. set the logical variable dvdfirst
      if (extended_optimal_level) then
         call matrixmpi (dvdfirst)
         call newcompi (wtaev)
      end if
      wtaev0 = 0.0
      dvdfirst = .false.
      do nit = 1, nscf
         if (mpi_rank == 0) write (*, '(/," Iteration number ",1i3,/," --------------------")') nit

!   For all pairs constrained through a Lagrange multiplier, compute
!   the Lagrange multiplier

         call setlagmpi (extended_optimal_level)

!   Improve all orbitals in turn

         dampmx = 0.0
         if (mpi_rank == 0) write (*, '(41x,"Self-            Damping"/,"Subshell    Energy    Method   P0    ",&
            &"consistency  Norm-1  factor  JP"," MTP INV NNP"/)')
         do orbital_index = 1, orbital_count
            jseq = iorder(orbital_index)
            if (is_orbital_fixed(jseq)) cycle
            call improvmpi (extended_optimal_level, jseq, lsort, dampmx)
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
               goto 3
            else
               if (scnsty(k) <= accy) then
                  convg = .true.
                  goto 3
               end if
            end if
            call improvmpi (extended_optimal_level, k, lsort, dampmx)
         end do

         call maxarr (k)

         if (k == 0) then
            convg = .true.
         else
            if (scnsty(k) <= accy) then
               convg = .true.
            else
               convg = .false.
            end if
         end if

    3    continue
         if (ldbpr(24) .and. mpi_rank==0) call prwf (0)

!   Perform Gram-Schmidt process
!   For OL calculation, orthst is true and orbitals are orthonormalized
!   in subroutine improv. For AL calculation, orthst is false.
         if (.not.orthst) call orthsc

!   Write the subshell radial wavefunctions to the .rwf file

         if (mpi_rank == 0) call orbout (output_wavefunction_path)

         if (extended_optimal_level) then
            call matrixmpi(dvdfirst)
            call newcompi (wtaev)
         end if
!        Make this a relative convergence test
!        IF(ABS(WTAEV-WTAEV0).LT.1.0D-9.and.
!    &   DAMPMX.LT.1.0D-4) CONVG=.true.
!        PRINT *, 'WTAEV, WTAEV0', WTAEV, WTAEV0,
!    &             ABS((WTAEV-WTAEV0)/WTAEV)
!GG         IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.001*ACCY) CONVG = .TRUE.
!cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
!cjb     IF(DABS(WTAEV-WTAEV0).LT.1.0D-8.and.              &
!cjb                   DAMPMX.LT.1.0D-2) CONVG=.true.
         if (abs((wtaev - wtaev0)/wtaev) < 0.001_real64*accy) convg = .true.
         wtaev0 = wtaev
         if (.not.convg) cycle
         if (ldbpr(25) .and. .not.ldbpr(24) .and. mpi_rank==0) call prwf (0)
            !IF (EOL) CALL matrixmpi (dvdfirst)
         goto 5

      end do

      if (mpi_rank==0) write (stderr_unit,*)" Maximum iterations in SCF Exceeded."

    5 continue
      do i = 31, 32 + max_mcp_file_multipole_rank
         close(i)                                ! The MCP coefficient files
      end do

      if (mpi_rank == 0) then
         !CLOSE (23)     ! The .rwf file
         close(25)                               ! The .mix file
      end if
!
!   Complete the summary - moved from rscf92 for easier alloc/dalloc
!
      if (mpi_rank == 0) call endsum
!
!   Deallocate storage
!
      call dalloc (wt, "WT", "SCFmpi")                       !Either getold or getald

      if (nec > 0) then
         call dalloc (iecc, "IECC", "SCFmpi")
         call dalloc (ecv, "ECV", "SCFmpi")
         call dalloc (csf_subshell_occupations, "IQA", "SCFmpi")
      end if

      if (nddim > 0) then
         call dalloc (da, "DA", "SCFmpi")
         call dalloc (nda, "NDA", "SCFmpi")
         nddim = 0
      end if

      if (nxdim > 0) then
         call dalloc (xa, "XA", "SCFmpi")
         call dalloc (nxa, "NXA", "SCFmpi")
         nxdim = 0
      end if

      if (nydim > 0) then
         call dalloc (ya, "YA", "SCFmpi")
         call dalloc (nya, "NYA", "SCFmpi")
         nydim = 0
      end if

      if (extended_optimal_level) then
         call dalloc (shifted_eigenvalues, "EVAL", "SCFmpi")
         call dalloc (csf_mixing_coefficients, "EvEC", "SCFmpi")
         call dalloc (iatjpo, "IATJPO", "SCFmpi")
         call dalloc (iaspar, "IASPAR", "SCFmpi")
         call dalloc (ncmaxblk, "NCMAXBLK", "SCFmpi") ! getold.f
         call dalloc (eavblk, "EAVBLK", "SCFmpi")     ! getold.f
         call dalloc (idxblk, "IDXBLK", "SCFmpi")     ! Allocated in getold.f
         call dalloc (iccmin, "ICCMIN", "SCFmpi")       ! Allocated in items.f<-getold.f
      end if
!
      call dalloc (csf_offsets_by_block, "NCFPAST", "SCFmpi")
      call dalloc (eigenvalue_offsets_by_block, "NCMINPAST", "SCFmpi")
      call dalloc (mixing_offsets_by_block, "NEVECPAST", "SCFmpi")


      return
      end subroutine scfmpi
