!***********************************************************************
!                                                                      *
      subroutine matrix(dvdfirst)
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
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:52:04   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double, long
      use memory_man
      use damp_c
      use def_c,           only: iccmin, ncmin, ncmax
      use debug_c
      use eigv_c
      use hblock_c
      use hmat_c
      use iounit_c
      use mcpa_c
      use mpi_s
      use orb_c
      use pos_c
      use peav_c
      use syma_c,          only: jpgg
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setham_i
      use maneig_i
      use suite_linear_algebra, only: dcopy, ddot, dscal
      implicit none
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
      real(double), dimension(:), pointer :: cmvl
      real(double) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
      logical :: first
      character :: mcplab*3

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
         call alloc (cmvl, mixing_coefficient_storage_size, 'CMVL', 'MATRIX')
         call dcopy (mixing_coefficient_storage_size, csf_mixing_coefficients, 1, cmvl, 1)
      endif

!=======================================================================
!   Position the files - MCP files (unit NFILE) for reading
!   and mixing coefficients file (unit 25) for writing
!=======================================================================

      do nfile = 30, 32 + max_mcp_file_multipole_rank
         rewind (nfile)
         if (nfile == 30) then
            read (nfile)
            read (nfile)
            read (nfile)
         endif
         read (nfile)
         read (nfile)
         read (nfile)
      end do

      ! To put in ncmin and nvecsiz. Values read here are the same
      ! as those from elsewhere (shch as common blocks)
      if (mpi_rank == 0) then
         rewind (25)
         read (25)                               ! 'G92MIX'
         read (25) nelec, total_csf_count, orbital_count, ntmp, ntmp, symmetry_block_count
         backspace (25)
         write (25) nelec, total_csf_count, orbital_count, ncmin, mixing_coefficient_storage_size, symmetry_block_count
      endif

!=======================================================================
!   Do the job block by block
!=======================================================================

               !------------------------------------------------
      do jblock = 1, symmetry_block_count                      ! block do-loop
               !------------------------------------------------

!=======================================================================
!   Read indeces of non-zero elements from mcp.30 file. Note the
!   format has been changed to lower-triangle-by-rows.
! Length of iendc can be reduced
!=======================================================================

         read (30) mcplab, jblockt, csf_count
         if (jblockt/=jblock .or. csf_count/=csf_counts_by_block(jblock)) error stop &
            'matrx: jblockt .NE. jblock .OR. ncf1 .NE. ncf2'
         read (30) nelmntgg
         stored_matrix_element_count = int(nelmntgg,kind=long)
         call alloc (matrix_partner_indices, stored_matrix_element_count, 'IROW', 'MATRIX')
         call alloc (hamiltonian_values, stored_matrix_element_count, 'EMT', 'MATRIX')
!cjb          ALLOC (IENDC, NCF + 1,...)  ->  ALLOC (IENDC, 0, NCF,...)
!cjb     CALL ALLOC (IENDC, NCF + 1, 'IENDC', 'MATRIX' )
         call alloc (stored_segment_end_offsets, 0, csf_count, 'IENDC', 'MATRIX' )

!        ! may not be necessary if iendc is ALWAYS used
         stored_segment_end_offsets(0:csf_count) = 0                        ! the way it is assigned here.

      !...EMT will be accumulated in setham
         hamiltonian_values(:stored_matrix_element_count) = 0.d0

         read (30) (stored_segment_end_offsets(i),i=mpi_rank + 1,csf_count,mpi_process_count), &
            & (matrix_partner_indices(i),i=1,stored_matrix_element_count)

         ncfpat = csf_offsets_by_block(jblock)
         ncminpat = eigenvalue_offsets_by_block(jblock)
         nevecpat = mixing_offsets_by_block(jblock)

!=======================================================================
!   Skip current block if no eigenlaue is required
!=======================================================================

         if (eigenvector_counts_by_block(jblock) == 0) then
            do nfile = 31, 32 + max_mcp_file_multipole_rank
               read (nfile) mcplab, jblockt, ncft, ncoeff
               if (jblockt /= jblock) error stop 'matrx: jblockt .NE. jblock'
               if (ncft /= csf_count) error stop 'matrx: ncft .NE. ncf'

               read (nfile) lab, ncontr
               do while(lab/=0 .or. ncontr/=0)
                  read (nfile) (itmp,itmp,tmp,i=1,ncontr)
                  read (nfile) lab, ncontr
               end do
            end do

            call dalloc (stored_segment_end_offsets, 'IENDC', 'MATRIX')
            call dalloc (matrix_partner_indices, 'IROW', 'MATRIX')

            cycle

         endif

!=======================================================================
!   Generate the Hamiltonian matrix - average energy is removed here
!=======================================================================

         call setham (jblock, mpi_rank, mpi_process_count)
!
!   Determine average energy
!
         hamiltonian_energy_shift = 0.d0
         do ir = mpi_rank + 1, csf_count, mpi_process_count
            hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(stored_segment_end_offsets(ir))
         end do
!         DO IR = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
!            EAV = EAV + EMT(IENDC(NPROCS*(IR-1)+MYID+1))
!         END DO

         hamiltonian_energy_shift = hamiltonian_energy_shift/csf_count
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

         call maneig (dvdfirst, ldbpg(3),                      &
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
               if (cdampj == 0.d0) cycle         ! So SURE ???

               omcdaj = 1.d0 - cdampj

           !...Damp eigenvector and determine the new dominant component
  123          continue
               amax = 0.d0
               do i = 1, csf_count
                  evecij = omcdaj*csf_mixing_coefficients(i+iofset) + cdampj*cmvl(i + iofset)
                  csf_mixing_coefficients(i+iofset) = evecij
                  wa = abs(evecij)
                  if (wa <= amax) cycle
                  amax = wa
                  ia = i
               end do

            !...compute the normalization factor
               sum = 0.d0
               do i = 1, csf_count
                  sum = sum + csf_mixing_coefficients(i+iofset)**2
               end do
               dnfac = 1.d0/sqrt(sum)

            !...Renormalize and invert as necessary
               if (csf_mixing_coefficients(ia+iofset) < 0.d0) dnfac = -dnfac
               call dscal (csf_count, dnfac, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)

            !...Schmidt orthogonalise
  234          continue
               jother = jother - 1
               if (jother < 1) cycle
               jofset = (jother - 1)*csf_count + nevecpat
               ovrlap = ddot(csf_count,csf_mixing_coefficients(iofset+1:iofset+(csf_count)),1, &
                  & csf_mixing_coefficients(jofset+1:jofset+(csf_count)),1)
               if (ovrlap /= 0.d0) then          ! So SURE ???
                  omcdaj = 1.d0
                  cdampj = -ovrlap
                  call dcopy (csf_count, csf_mixing_coefficients(jofset+1:jofset+(csf_count)), 1, cmvl(iofset + &
                     & 1:iofset+(csf_count)), 1)
                  go to 123
               else
                  go to 234
               endif
            end do
         endif

!   Write out the eigenpair information: ASF symmetries, eigenvalues,
!   and eigenvectors to GRASP92 mixing coefficients File

         if (eigenvector_counts_by_block(jblock) == 0) then
            iattmp = 999
            iastmp = 999
         else
!GGGG
            iattmp = iabs(jpgg(jblock))
            if(jpgg(jblock) >= 0) then
!cjb iastmp error corrected by Pawel Syty
!PS            iastmp = 2
               iastmp = 1
            else
!cjb iastmp error corrected by Pawel Syty
!PS            iastmp = 1
               iastmp = -1
            end if
!GG            IATTMP = IATJPO(NCMINPAT + 1)
!GG            IASTMP = IASPAR(NCMINPAT + 1)
         endif

         if (mpi_rank == 0) then
            write (25) jblock, csf_count, eigenvector_counts_by_block(jblock), iattmp, iastmp
            write (25) (iccmin(i + ncminpat),i=1,eigenvector_counts_by_block(jblock))
            write (25) hamiltonian_energy_shift, (shifted_eigenvalues(i + ncminpat),i=1,eigenvector_counts_by_block(jblock))
            write (25) ((csf_mixing_coefficients(i+(j-1)*csf_count+nevecpat),i=1,csf_count),j=1,eigenvector_counts_by_block(jblock)&
               )
         endif

         call dalloc (hamiltonian_values, 'EMT', 'MATRIX')
         call dalloc (stored_segment_end_offsets, 'IENDC', 'MATRIX')
         call dalloc (matrix_partner_indices, 'IROW', 'MATRIX')

!----------------------
      end do
!----------------------
!
!   Deallocate the temporary storage
!
      if (.not.first) then
         call dalloc (cmvl, 'CMVL', 'MATRIX')
      else
         first = .false.
      endif

  302 format(' Average energy = ',1p,d18.10,' Hartrees')

      return
      end subroutine matrix
