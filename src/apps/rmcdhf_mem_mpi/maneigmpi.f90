!***********************************************************************
!                                                                      *
      subroutine maneigmpi(dvdfirst, print_details, symmetry_block_index, &
                        ncfpat, ncminpat, nevecpat, total_csf_count)
!                                                                      *
!   This module  manages the  operation of the  eigensolvers and the   *
!   storage of the eigenpairs.  There are two principal branches:      *
!                                                                      *
!      (1) Matrix of order 1: the trivial case                         *
!      (2) Matrix of order greater than 1: eigenpairs are found        *
!          using DVDSON; this involves up to three steps:              *
!                    (a) The matrix is analysed to determine its       *
!                        block structure (only irreducibe matrices     *
!                        are correctly treated by DVDSON)              *
!                    (b) Eigenpairs are extracted for each block       *
!                    (c) The appropriate eigenpairs are selected and   *
!                        stored                                        *
!                                                                      *
!   We  assume that  the sparse representation  of the matrix  is in   *
!   core.                                                              *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, RALLOC.          *
!               [RSCF92]: POSNFL, SPICMV2.
!               [DVDSON]: DVDSON                                       *
!               [AUXBLAS]: DINIT/SINIT                                 *
!               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL, DSWAP/SSWAP          *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
!   Modified by Xinghong He               Last revision: 17 Aug 1998   *
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
      use atomic_calculation_state, only: iccmin
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, &
         & mpi_process_count, mpi_rank, mpi_sum
      use orbital_csf_basis, only: csf_count
      use level_symmetry_metadata,          only: jpgg
      use active_symmetry_block, only: active_block_index => jblock
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mpi_sparse_matrix_product, only: spicmvmpi
      use mpi_initial_eigenvectors, only: iniestmpi
      use davidson_eigensolver, only: dvd_operator, gdvd
      use itjpo_i, only: itjpo
      use ispar_i, only: ispar
      use blas_lapack_interfaces, only: dcopy, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) ::  dvdfirst
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: ncfpat
      integer, intent(in) :: ncminpat
      integer, intent(in) :: nevecpat
      integer             :: total_csf_count
      logical             :: print_details
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvect, subspace_limit, real_workspace_size, nvex, initial_vector_count, iteration_limit, &
         & submatrix_order_limit, expansion_block_size, &
         lowest_pair_index, highest_pair_index, integer_workspace_size, ic, iteration_count, matrix_vector_product_count, nend, &
            & level_index, jstate, &
         iofset, csf_index, ia
      integer, pointer :: integer_workspace(:), selected_pair_indices(:)
      real(real64) :: pnwork, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, orthogonality_tolerance, amax, &
         & wa, dnfac
      real(real64), pointer :: real_workspace(:), matrix_diagonal(:), atmp(:)
      logical :: solve_high_end
!-----------------------------------------------

      !PRINT *, 'maneig ...'

!      ...spicmvmpi needs this COMMON /WCHBLK/JBLOCKK
      active_block_index = symmetry_block_index
!
!=======================================================================
!   Trivial case
!=======================================================================
      if (csf_count == 1) then
         shifted_eigenvalues(ncminpat+1) = 0.0_real64
         csf_mixing_coefficients(nevecpat+1) = 1.0_real64
         goto 123                               ! Don't like big ELSE
      end if
!
!=======================================================================
!   Non-trivial case - Use Davidson eigensolver
!=======================================================================
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECT
!   plus a smaller number if storage is severely constrained
!
      nvect = ncmaxblk(symmetry_block_index)
!CFF   ... make this more like rscf
      subspace_limit   = min (csf_count, 2*nvect + 80)
!GG      lim   = MIN (ncf, 2*nvect + 40)
!      LIM = MIN(NCF,2*NVECT + 20)
      real_workspace_size = 2*csf_count*subspace_limit + subspace_limit*subspace_limit*2 + 11*subspace_limit + nvect
      call alloc (real_workspace, real_workspace_size, "WORK", "MANEIGmpi")

      !...At most 14 ? restriction removed xhh 98-05-19
      !nvex = MIN (nvect,ncfblk(jblock),14)
      nvex = min(nvect,csf_counts_by_block(symmetry_block_index))
      initial_vector_count = nvex
      iteration_limit = min(nvect*200,csf_count)
!      N1000 = 2000
      submatrix_order_limit = 4000
!
!   Initial estimates for eigenvectors
!
!CFF
      if (dvdfirst .or. (csf_count <= submatrix_order_limit) ) then
!GG         CALL INIESTmpi (N1000, NCF, NIV, WORK, EMT, IENDC, IROW)
         call iniestmpi (submatrix_order_limit, csf_count, initial_vector_count, real_workspace, hamiltonian_values, &
            & stored_segment_end_offsets(0:csf_count), &
            & matrix_partner_indices)
      else
!CFF        .. use current estimates
      nend = csf_count * nvex
         do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
            real_workspace( nend + iccmin(level_index+ncminpat) ) = shifted_eigenvalues(ncminpat+level_index)
            call dcopy (csf_count, csf_mixing_coefficients( nevecpat + csf_count*(level_index-1) + 1:nevecpat + csf_count*( &
               & level_index-1)+ &
               & (csf_count)),    1,    &
                          real_workspace( csf_count*(iccmin(level_index+ncminpat)-1) + 1 :( csf_count*(iccmin(level_index+ &
                             & ncminpat)-1) + 1 )+((csf_count)-1)), 1)
         end do
      end if

! iniest looks for eigenvectors of n1000*n1000 matrix so there
! is no need to call dvdson if block size <= n1000

      if (csf_count > submatrix_order_limit) then
         if (mpi_rank == 0) write (*,*) "Calling dvdson!!!", iteration_limit,nvect

         ! Call Davidson eigensolver

         expansion_block_size = 1
         lowest_pair_index = 1
         highest_pair_index = nvex
         integer_workspace_size = 6*subspace_limit + nvect
         eigenvalue_tolerance = 1.0e-17_real64
!         CRITC = 1.0D-08
!         CRITR = 1.0D-08
!         ORTHO = MAX(1D-8,CRITR)
         coefficient_tolerance = 1.0e-9_real64
         residual_tolerance = 1.0e-9_real64
         orthogonality_tolerance = max (1e-9_real64, residual_tolerance)
!
!   Store the diagonals in a separate array and make it global
!
         call alloc (matrix_diagonal, csf_count, "DIAG", "MANEIGmpi")
         call alloc (atmp, csf_count, "ATMP", "MANEIGmpi")
         do csf_index = 1, csf_count
            atmp(csf_index) = 0.0_real64
            matrix_diagonal(csf_index) = 0.0_real64    ! this one may not be necessary
         end do

         do ic = mpi_rank + 1, csf_count, mpi_process_count
!GG            DIAG(IC) = EMT(IENDC(IC))
            atmp(ic) = hamiltonian_values(stored_segment_end_offsets(ic))
         end do
         call mpi_allreduce (atmp, matrix_diagonal, csf_count, mpi_double_precision,  &
                             mpi_sum, mpi_comm_world, mpi_error_code)
         call dalloc(atmp,"ATMP", "MANEIGmpi")

         call alloc (integer_workspace, integer_workspace_size, "IWORK", "MANEIGmpi")
         call alloc (selected_pair_indices, subspace_limit, "JWORK", "MANEIGmpi")
         if (csf_count>1000) then
            call gdvd (spicmvmpi,csf_count,subspace_limit,matrix_diagonal,lowest_pair_index,highest_pair_index, &
               & selected_pair_indices,initial_vector_count,expansion_block_size, &
            eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, orthogonality_tolerance, iteration_limit, &
               & real_workspace, real_workspace_size, integer_workspace, integer_workspace_size, &
            solve_high_end, iteration_count, matrix_vector_product_count, mpi_error_code)
         end if
         call dalloc (matrix_diagonal, "DIAG", "MANEIGmpi")
         call dalloc (integer_workspace, "IWORK", "MANEIGmpi")
         call dalloc (selected_pair_indices, "JWORK", "MANEIGmpi")

         if (mpi_rank == 0) then
            write (*,'("DVDSON: ",1i3," loops; ",1i3," matrix-vector multiplies.")') iteration_count, matrix_vector_product_count
            if (mpi_error_code /= 0) then
               write (*,'(" Returned from DVDSON with IERR = ",1i4)') mpi_error_code
            end if
         end if
      end if
!
!   Pick up the eigen pairs and store in EVAL and EVEC
!
      nend = csf_count*nvex
      do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
         shifted_eigenvalues(ncminpat+level_index) = real_workspace(nend + iccmin(level_index + ncminpat))
         call dcopy (csf_count, real_workspace(csf_count*(iccmin(level_index + ncminpat) - 1) + 1:csf_count*(iccmin(level_index &
            & + ncminpat) - 1)+ &
            & (csf_count)), 1, &
      & csf_mixing_coefficients( nevecpat+csf_count*(level_index-1)+1:nevecpat+csf_count*(level_index-1)+(csf_count)), 1)
      end do
!     print *, ncminpat,(eval(ncminpat+j),j=1,nevblk(jblock)),
!    1 'zou,from maneig'
!
!   Deallocate storage
!
      call dalloc (real_workspace, "WORK", "MANEIGmpi")

  123 continue
      do jstate = 1, eigenvector_counts_by_block(symmetry_block_index)
!
!   Find the dominant component of each eigenvector
!
         iofset = nevecpat + csf_count*(jstate - 1)

         amax = 0.0_real64
         do csf_index = 1, csf_count
            wa = abs(csf_mixing_coefficients(csf_index+iofset))
            if (wa <= amax) cycle
            amax = wa
            ia = csf_index
         end do
!
!   Find the angular momentum and parity of the dominant component
!
!GG         IATJPO(JSTATE+NCMINPAT) = ITJPO(IA + NCFPAT)
!GG         IASPAR(JSTATE+NCMINPAT) = ISPAR(IA + NCFPAT)
!
!   Redefine eigenvectors so that the dominant component
!   is positive
!
         if (csf_mixing_coefficients(ia+iofset) >= 0.0_real64) cycle
         dnfac = -1.0_real64
         call dscal (csf_count, dnfac, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)
!===============================================================

      end do


      return
      end subroutine maneigmpi
