!***********************************************************************
!                                                                      *
      subroutine matrixmpi(dvdfirst)
      use mcp_file_metadata, only: mcp_input_units
      use, intrinsic :: iso_fortran_env, only: output_unit
      use calculation_file_units, only: mixing_output_unit
!                                                                      *
!   Calls routines to form the Hamiltonian matrix and to diagonalise   *
!   it. The total angular momenta and parity of each  ASF  is found;   *
!   the eigenvectors are normalised so  that the sign of the largest   *
!   element is positive.                                               *
!                                                                      *
!   Call(s) to: [LIB92] ALLOC, DALLOC.                                 *
!               [RSCF92]: SETHAM, HMOUT.                               *
!               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL                       *
!                                                                      *
!                                         Last revision: 24 Dec 1992   *
! Block version by Xinghong He                           07 Aug 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use array_memory_management, only: alloc, dalloc
      use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
      use atomic_calculation_state, only: iccmin, ncmin, ncmax
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mcp_file_metadata, only: max_mcp_file_multipole_rank
      use mpi_state, only: mpi_allreduce, mpi_barrier, mpi_comm_world, mpi_double_precision, mpi_error_code, &
         & mpi_process_count, mpi_rank, mpi_sum
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use block_mean_energies, only: eavblk
      use level_symmetry_metadata,          only: jpgg
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setham_i, only: setham
      use maneigmpi_i, only: maneigmpi
      use blas_lapack_interfaces, only: dcopy, ddot, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) ::  dvdfirst
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfile, nelec, ntmp, symmetry_block_index, jblockt, i, ncfpat, ncminpat, &
         nevecpat, ncft, ncoeff, lab, ncontr, itmp, csf_index, idiag, level_index, iofset, &
         jother, ia, jofset, iattmp, iastmp, nelmntgg
      real(real64), pointer :: cmvl(:)
      real(real64) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
      logical :: first
      character :: mcp_channel_label*3

      save first
!
      data first/ .true./
!
!     POINTER (cmvl(1))
!
!-----------------------------------------------------------------------

      if (mpi_rank == 0) write (output_unit, *)

!   Allocate memory for CMVL once (the maximum size)
!   Save previous estimate of eigenvectors

      if (.not.first) then
         call alloc (cmvl, mixing_coefficient_storage_size, "CMVL", "MATRIXmpi")
         call dcopy (mixing_coefficient_storage_size, csf_mixing_coefficients, 1, cmvl, 1)
      end if

!=======================================================================
!   Position the files - MCP files (unit NFILE) for reading
!   and mixing coefficients file (unit 25) for writing
!=======================================================================

      do nfile = 30, 32 + max_mcp_file_multipole_rank
         rewind (mcp_input_units(nfile))
         if (nfile == 30) then
            read (mcp_input_units(nfile))
            read (mcp_input_units(nfile))
            read (mcp_input_units(nfile))
         end if
         read (mcp_input_units(nfile))
         read (mcp_input_units(nfile))
         read (mcp_input_units(nfile))
      end do

      ! To put in ncmin and nvecsiz. Values read here are the same
      ! as those from elsewhere (shch as common blocks)
      if (mpi_rank == 0) then
         rewind (mixing_output_unit)
         read (mixing_output_unit)                               ! 'G92MIX'
         read (mixing_output_unit) nelec, total_csf_count, orbital_count, ntmp, ntmp, symmetry_block_count
         backspace (mixing_output_unit)
         write (mixing_output_unit) nelec, total_csf_count, orbital_count, ncmin, mixing_coefficient_storage_size, &
            symmetry_block_count
      end if

!=======================================================================
!   Do the job block by block
!=======================================================================

               !------------------------------------------------
      do symmetry_block_index = 1, symmetry_block_count                      ! block do-loop
               !------------------------------------------------

!=======================================================================
!   Read indeces of non-zero elements from mcp.30 file. Note the
!   format has been changed to lower-triangle-by-rows.
! Length of iendc can be reduced
!=======================================================================

         read (mcp_input_units(30)) mcp_channel_label, jblockt, csf_count
         if (jblockt/=symmetry_block_index .or. csf_count/=csf_counts_by_block(symmetry_block_index)) error stop &
            "MATRIXMPI: jblockt .NE. jblock .OR. ncf1 .NE. ncf2"
         read (mcp_input_units(30)) nelmntgg
         stored_matrix_element_count = int(nelmntgg,kind=int64)
         call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "MATRIXmpi")
         call alloc (hamiltonian_values, stored_matrix_element_count, "EMT", "MATRIXmpi")
        call alloc (stored_segment_end_offsets, 0, csf_count, "IENDC", "MATRIXMPI" )

!                                                ! may not be necessary if iendc is ALWAYS used
!GG         IENDC(:NCF) = 0                        ! the way it is assigned here.
         stored_segment_end_offsets(0:csf_count) = 0

      !...EMT will be accumulated in setham
         hamiltonian_values(:stored_matrix_element_count) = 0.0_real64

         read (mcp_input_units(30)) (stored_segment_end_offsets(i),i=mpi_rank+1,csf_count,mpi_process_count), &
            (matrix_partner_indices(i) &
            & ,i=1,stored_matrix_element_count)

         ncfpat = csf_offsets_by_block(symmetry_block_index)
         ncminpat = eigenvalue_offsets_by_block(symmetry_block_index)
         nevecpat = mixing_offsets_by_block(symmetry_block_index)

!=======================================================================
!   Skip current block if no eigenlaue is required
!=======================================================================

         if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
            do nfile = 31, 32 + max_mcp_file_multipole_rank
               read (mcp_input_units(nfile)) mcp_channel_label, jblockt, ncft, ncoeff
               if (jblockt /= symmetry_block_index)                  &
                          error stop "MATRIXmpi: jblockt .NE. jblock"
               if (ncft /= csf_count) error stop "MATRIXmpi: ncft .NE. ncf"

               read (mcp_input_units(nfile)) lab, ncontr
               do while(lab/=0 .or. ncontr/=0)
                  read (mcp_input_units(nfile)) (itmp,itmp,tmp,i=1,ncontr)
                  read (mcp_input_units(nfile)) lab, ncontr
               end do
            end do

            call dalloc (stored_segment_end_offsets, "IENDC", "MATRIXmpi")
            call dalloc (matrix_partner_indices, "IROW", "MATRIXmpi")

            cycle

         end if

!=======================================================================
!   Generate the Hamiltonian matrix - average energy is removed here
!=======================================================================

         call setham (symmetry_block_index, mpi_rank, mpi_process_count)
!
!   Determine average energy
!
         hamiltonian_energy_shift = 0.0_real64
         do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(stored_segment_end_offsets(csf_index))
         end do
         call mpi_barrier(mpi_comm_world,mpi_error_code)
         call mpi_allreduce (hamiltonian_energy_shift,tmp,1,mpi_double_precision,mpi_sum,  &
                             mpi_comm_world,mpi_error_code)
         hamiltonian_energy_shift = tmp / csf_count
!GG         EAV = EAV/NCF
         eavblk(symmetry_block_index) = hamiltonian_energy_shift

      ! Print Hamiltonian matrix and average energy
      ! hmout is not general
      !call hmout (0, 1, ncf)

         if (mpi_rank == 0) write (*, '(" Average energy = ",1p,d18.10," Hartrees")') hamiltonian_energy_shift

      ! Subtract the average energy from the diagonal elements
      ! to reduce the condition number of the matrix
!         DO I = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
!            EMT(IENDC(NPROCS*(I-1)+MYID+1)) = EMT(IENDC(NPROCS*(I-1)+MYID+1))&
!                - EAV
         do i = mpi_rank + 1, csf_count, mpi_process_count
            idiag = stored_segment_end_offsets(i)     ! new mode: each row ends in diagonal
            hamiltonian_values(idiag) = hamiltonian_values(idiag) - hamiltonian_energy_shift
         end do

!=======================================================================
!   Compute and store eigenpairs
!=======================================================================
         call maneigmpi (dvdfirst, ldbpg(3),                      &
                   symmetry_block_index, ncfpat, ncminpat, nevecpat, total_csf_count)

!=======================================================================
!   Damp and Schmidt orthogonalise eigenvectors for OL calculations
!=======================================================================

         if (.not.first) then

            do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)

               iofset = (level_index - 1)*csf_count + nevecpat
               jother = level_index

            ! cdamp has the original non-block feature
               cdampj = mixing_damping_factors(level_index + ncminpat)
               if (cdampj == 0.0_real64) cycle         ! So SURE ???

               omcdaj = 1.0_real64 - cdampj

           !...Damp eigenvector and determine the new dominant component
  123          continue
               amax = 0.0_real64
               do i = 1, csf_count
                  evecij = omcdaj*csf_mixing_coefficients(i+iofset) + cdampj*cmvl(i + iofset)
                  csf_mixing_coefficients(i+iofset) = evecij
                  wa = abs(evecij)
                  if (wa <= amax) cycle
                  amax = wa
                  ia = i
               end do

            !...compute the normalization factor
               sum = 0.0_real64
               do i = 1, csf_count
                  sum = sum + csf_mixing_coefficients(i+iofset)**2
               end do
               dnfac = 1.0_real64/sqrt(sum)

            !...Renormalize and invert as necessary
               if (csf_mixing_coefficients(ia+iofset) < 0.0_real64) dnfac = -dnfac
               call dscal (csf_count, dnfac, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)

            !...Schmidt orthogonalise
  234          continue
               jother = jother - 1
               if (jother < 1) cycle
               jofset = (jother - 1)*csf_count + nevecpat
               ovrlap = ddot(csf_count,csf_mixing_coefficients(iofset+1:iofset+(csf_count)),1, &
                  & csf_mixing_coefficients(jofset+1:jofset+(csf_count)),1)
               if (ovrlap /= 0.0_real64) then          ! So SURE ???
                  omcdaj = 1.0_real64
                  cdampj = -ovrlap
                  call dcopy (csf_count, csf_mixing_coefficients(jofset+1:jofset+(csf_count)), 1, cmvl(iofset + &
                     & 1:iofset+(csf_count)), 1)
                  goto 123
               else
                  goto 234
               end if
            end do
         end if

!   Write out the eigenpair information: ASF symmetries, eigenvalues,
!   and eigenvectors to GRASP92 mixing coefficients File

         if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
            iattmp = 999
            iastmp = 999
         else
!GGGG
         iattmp = iabs(jpgg(symmetry_block_index))
         if(jpgg(symmetry_block_index) >= 0) then
            iastmp = 1
         else
            iastmp = -1
         end if
!GG            IATTMP = IATJPO(NCMINPAT + 1)
!GG            IASTMP = IASPAR(NCMINPAT + 1)
         end if

         if (mpi_rank == 0) then
            write (mixing_output_unit) symmetry_block_index, csf_count, eigenvector_counts_by_block(symmetry_block_index), &
               & iattmp, iastmp
            write (mixing_output_unit) (iccmin(i + ncminpat),i=1,eigenvector_counts_by_block(symmetry_block_index))
            write (mixing_output_unit) hamiltonian_energy_shift, (shifted_eigenvalues(i + ncminpat),i=1, &
               eigenvector_counts_by_block(symmetry_block_index))
            write (mixing_output_unit) ((csf_mixing_coefficients(i+(level_index-1)*csf_count+nevecpat),i=1,csf_count), &
               & level_index=1, &
               eigenvector_counts_by_block(symmetry_block_index)&
               )
         end if

         call dalloc (hamiltonian_values, "EMT", "MATRIXmpi")
         call dalloc (stored_segment_end_offsets, "IENDC", "MATRIXmpi")
         call dalloc (matrix_partner_indices, "IROW", "MATRIXmpi")

!----------------------
      end do
!----------------------
!
!   Deallocate the temporary storage
!
      if (.not.first) then
         call dalloc (cmvl, "CMVL", "MATRIXmpi")
      else
         first = .false.
      end if


      return
      end subroutine matrixmpi
