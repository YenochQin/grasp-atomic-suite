!***********************************************************************
!                                                                      *
      subroutine matrixmpi(dvdfirst)
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
!   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
!      Spin-angular coefficiants have been put on memory               *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use memory_man
      use damp_c, only: double, max_orbitals, mixing_damping_factors, orbital_damping_factors
      use def_c, only: iccmin, ncmin, ncmax
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use eigv_c
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use hmat_c
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use mcpa_c, only: double, max_mcp_file_multipole_rank
      use mpi_c
      use orb_c
      use pos_c
      use peav_c, only: double, eavblk
      use syma_c,          only: jpgg
      use rmcdhf_mem_c,    only: imin1_30, imax1_30,imin2_30, imax2_30
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setham_i, only: setham
      use maneigmpi_i, only: maneigmpi
      use read1_mem_i, only: read1_mem
      use read03_mem_i, only: read03_mem
      use suite_linear_algebra, only: dcopy, ddot, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) ::  dvdfirst
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfile, nelec, ntmp, jblock, jblockt, i, ncfpat, ncminpat, &
         nevecpat, ncft, ncoeff, lab, ncontr, itmp, ir, idiag, j, iofset, &
         jother, ia, jofset, iattmp, iastmp, nelmntgg
      real(real64), pointer :: cmvl(:)
      real(real64) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
      logical :: first
      character :: mcplab*3
      integer(int64) :: ncontr_tot(20), ncontr_tot2(20)

      save first
!
      data first/ .true./
!
!     POINTER (cmvl(1))
!
!-----------------------------------------------------------------------

      if (mpi_rank == 0) write (6, *)

!   Allocate memory for CMVL once (the maximum size)
!   Save previous estimate of eigenvectors

      if (.not.first) then
         call alloc (cmvl, mixing_coefficient_storage_size, "CMVL", "MATRIXmpi")
         call dcopy (mixing_coefficient_storage_size, csf_mixing_coefficients, 1, cmvl, 1)
      end if

      ! To put in ncmin and nvecsiz. Values read here are the same
      ! as those from elsewhere (shch as common blocks)
      if (mpi_rank == 0) then
         rewind (25)
         read (25)                               ! 'G92MIX'
         read (25) nelec, total_csf_count, orbital_count, ntmp, ntmp, symmetry_block_count
         backspace (25)
         write (25) nelec, total_csf_count, orbital_count, ncmin, mixing_coefficient_storage_size, symmetry_block_count
      end if

!=======================================================================
!   Do the job block by block
!=======================================================================

      ncontr_tot = 0
      ncontr_tot2 = 0
               !------------------------------------------------
      do jblock = 1, symmetry_block_count                      ! block do-loop
               !------------------------------------------------

!=======================================================================
!   Read indeces of non-zero elements from mcp.30 file. Note the
!   format has been changed to lower-triangle-by-rows.
! Length of iendc can be reduced
!=======================================================================

         nelmntgg = imax2_30(jblock) - imin2_30(jblock)
         stored_matrix_element_count = int(nelmntgg,kind=int64)
         csf_count = imax1_30(jblock) - imin1_30(jblock)
         call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "MATRIXmpi")
         call alloc (hamiltonian_values, stored_matrix_element_count, "EMT", "MATRIXmpi")
         call alloc (stored_segment_end_offsets, 0, csf_count, "IENDC", "MATRIXMPI" )

!                                                ! may not be necessary if iendc is ALWAYS used
!GG         IENDC(:NCF) = 0                        ! the way it is assigned here.
         stored_segment_end_offsets(0:csf_count) = 0

      !...EMT will be accumulated in setham
         hamiltonian_values(:stored_matrix_element_count) = 0.0_real64
         call read03_mem(jblock)

         ncfpat = csf_offsets_by_block(jblock)
         ncminpat = eigenvalue_offsets_by_block(jblock)
         nevecpat = mixing_offsets_by_block(jblock)

!=======================================================================
!   Skip current block if no eigenlaue is required
!=======================================================================

         if (eigenvector_counts_by_block(jblock) == 0) then
            do nfile = 31, 32 + max_mcp_file_multipole_rank
               ncontr_tot2(nfile-30) = ncontr_tot2(nfile-30) + 1
               call read1_mem(nfile, ncontr_tot2(nfile-30), lab, ncontr)
               do while(lab/=0 .or. ncontr/=0)
                  ncontr_tot(nfile-30) = ncontr_tot(nfile-30) + ncontr
                  ncontr_tot2(nfile-30) = ncontr_tot2(nfile-30) + 1
                  call read1_mem(nfile, ncontr_tot2(nfile-30), lab, ncontr)
               end do
            end do

            call dalloc (stored_segment_end_offsets, "IENDC", "MATRIXmpi")
            call dalloc (matrix_partner_indices, "IROW", "MATRIXmpi")

            cycle

         end if

!=======================================================================
!   Generate the Hamiltonian matrix - average energy is removed here
!=======================================================================

         call setham (jblock, mpi_rank, mpi_process_count, ncontr_tot, ncontr_tot2)
!
!   Determine average energy
!
         hamiltonian_energy_shift = 0.0_real64
         do ir = mpi_rank + 1, csf_count, mpi_process_count
            hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(stored_segment_end_offsets(ir))
         end do
         call mpi_barrier(mpi_comm_world,mpi_error_code)
         call mpi_allreduce (hamiltonian_energy_shift,tmp,1,mpi_double_precision,mpi_sum,  &
                             mpi_comm_world,mpi_error_code)
         hamiltonian_energy_shift = tmp / csf_count
!GG         EAV = EAV/NCF
         eavblk(jblock) = hamiltonian_energy_shift

      ! Print Hamiltonian matrix and average energy
      ! hmout is not general
      !call hmout (0, 1, ncf)

         if (mpi_rank == 0) write (*, 302) hamiltonian_energy_shift

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
                   jblock, ncfpat, ncminpat, nevecpat, total_csf_count)

!=======================================================================
!   Damp and Schmidt orthogonalise eigenvectors for OL calculations
!=======================================================================

         if (.not.first) then

            do j = 1, eigenvector_counts_by_block(jblock)

               iofset = (j - 1)*csf_count + nevecpat
               jother = j

            ! cdamp has the original non-block feature
               cdampj = mixing_damping_factors(j + ncminpat)
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

         if (eigenvector_counts_by_block(jblock) == 0) then
            iattmp = 999
            iastmp = 999
         else
!GGGG
         iattmp = iabs(jpgg(jblock))
         if(jpgg(jblock) >= 0) then
            iastmp = 1
         else
            iastmp = -1
         end if
!GG            IATTMP = IATJPO(NCMINPAT + 1)
!GG            IASTMP = IASPAR(NCMINPAT + 1)
         end if

         if (mpi_rank == 0) then
            write (25) jblock, csf_count, eigenvector_counts_by_block(jblock), iattmp, iastmp
            write (25) (iccmin(i + ncminpat),i=1,eigenvector_counts_by_block(jblock))
            write (25) hamiltonian_energy_shift, (shifted_eigenvalues(i + ncminpat),i=1,eigenvector_counts_by_block(jblock))
            write (25) ((csf_mixing_coefficients(i+(j-1)*csf_count+nevecpat),i=1,csf_count),j=1,eigenvector_counts_by_block(jblock)&
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

  302 format(" Average energy = ",1p,d18.10," Hartrees")

      return
      end subroutine matrixmpi
