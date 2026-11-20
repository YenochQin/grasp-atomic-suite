!***********************************************************************
!                                                                      *
      subroutine maneig(iatjpo, iaspar)
!                                                                      *
!   This module  manages the  operation of the  eigensolvers and the   *
!   storage of the eigenpairs.  There are two principal branches:      *
!                                                                      *
!      (1) Matrix of order 1: the trivial case                         *
!      (2) Matrix of order greater than 1: there are two branches      *
!             (i) Matrices of order less than or equal to IOLPCK:      *
!                 eigenpairs are found using LAPACK SUBROUTINEs        *
!            (ii) Matrices of order greater than IOLPCK: eigenpairs    *
!                 are found using DVDSON; this involves up to three    *
!                 steps:                                               *
!                    (a) The matrix is analysed to determine its       *
!                        block structure (only irreducibe matrices     *
!                        are correctly treated by DVDSON)              *
!                    (b) Eigenpairs are extracted for each block       *
!                    (c) The appropriate eigenpairs are selected and   *
!                        stored                                        *
!                 Different methods of storage and different           *
!                 versions of the matrix-vector multiply are used      *
!                 depending upon the order and density of the matrix   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, posfile,         *
!                        RALLOC.                                       *
!               [RCI92]: DNICMV, SPICMV2, SPODMV.                      *
!               [DVDSON]: DVDSON.                                      *
!               [AUXBLAS]: DINIT/SINIT.                                *
!               [BLAS]: DCOPY/SCOPY, DSWAP/SSWAP.                      *
!               [LAPACK]: DSPEVX/SSPEVX.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
!   Modified Misha Saparov                                  Feb 1997   *
!            Charlotte F. Fischer                           May 1997   *
!                 Except for the disk version, all matrices have       *
!                 diagonals  shifted by EAV                            *
!  All arrays allocated here are de-allocated except pneval and pnevec *
!  which will be de-allocated in matrix.                               *
!  Block Version  By Xinghong He          Last revision: 18 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use array_memory_management, only: alloc, dalloc, ralloc
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use rci_matrix_file_position, only: nposition
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use orbital_csf_basis,           only: csf_count, orbital_count, csf_subshell_occupations
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use where_c, only: imcdf
      use active_symmetry_block, only: active_block_index => jblock
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dnicmv_i, only: dnicmv
      use sparse_symmetric_matvec, only: spicmv2
      use spodmv_i, only: spodmv
      use file_positioning, only: posfile
      use real_vector_initialization, only: dinit
      use dspevx_i, only: dspevx
      use iniestsd_i, only: iniestsd
      use davidson_eigensolver, only: dvd_operator, gdvd
      use sparse_eigensolver_initial_basis, only: iniest2
      use iniestdm_i, only: iniestdm
      use csf_quantum_queries, only: itjpo, ispar
      use blas_lapack_interfaces, only: dcopy, dspevx, dswap
      use blas_lapack_interfaces, only: dlamch
      implicit none (type, external)
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: iatjpo
      integer, intent(out) :: iaspar
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!cjb  INTEGER, PARAMETER      :: IOLPCK = 1000
      integer, parameter      :: iolpck = 2000
! GG      REAL(DOUBLE), PARAMETER :: ABSTOL = 1.0D-10
!cjb  NINCOR
!cjb  INTEGER, PARAMETER      :: NINCOR = 1         ! To enforce DISK
      integer, parameter      :: nincor = 268435456 ! = 2 GB or  memory
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer(int64) :: nstore
      integer :: mpi_rank, mpi_process_count, row_count, i, ndense, file_csf_count, iccutdum, file_mpi_rank, &
         file_mpi_process_count, iofset, nelc, ir, nvecmn, nvex, m, info, loc, nbrkev,   &
         imv, ndense_l, subspace_limit, real_workspace_size, integer_workspace_size, iteration_limit, expansion_block_size, nend, &
         lowest_pair_index, highest_pair_index, initial_vector_count, iteration_count, matrix_vector_product_count, ierr, j, ia
      real(real64) :: elsto,  dummy, &
          diatmp, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, orthogonality_tolerance, dmungo, amax, wa, abstol
! GG          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA
      logical :: solve_high_end, ldisc, sparse
      character(len=8) :: number_text
      real(real64), pointer :: w(:), z(:), real_workspace(:), matrix_diagonal(:)
      integer, pointer :: integer_workspace(:), ifail(:), selected_pair_indices(:)
!-----------------------------------------------------------------------
      abstol = 2*dlamch("S")
      mpi_rank = 0
      mpi_process_count = 1
      !IF (MYID == 0) WRITE (6, *) 'Calling maneig...'

!  (nrows+1) is the number of records of the present block's .res file

      row_count = (csf_count - mpi_rank - 1 + mpi_process_count)/mpi_process_count
      if (csf_count < mpi_process_count) row_count = csf_count/(mpi_rank + 1)
      !CALL posfile (1, imcdf, nrows+1)
      call posfile (0, imcdf, nposition)

      if (csf_count == 1) then
!-----------------------------------------------------------------------
!
! (1) - Trivial    ncf = 1
!
!-------------------------------------------------------
         write (24, *) "Trivial eigenvalue problem."

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         call alloc (shifted_eigenvalues, 1,"EVAL", "MANEIG" )
         call alloc (csf_mixing_coefficients, 1, "EVECO", "MANEIG")
         shifted_eigenvalues(1) = 0.0_real64
         csf_mixing_coefficients(1) = 1.0_real64

         ! Still read through the .res file
!GG
!GG Gediminas NIST 2005.11.03
!GG         READ (imcdf)
         do i = 1, row_count + 1
            read (imcdf)
         end do

      else                                       !if-2
!-----------------------------------------------------------------------
!
! (2) - Non trivial
!
!-------------------------------------------------------
!
!   Matrix of order greater than 1; how many elements in a triangle?
!
         ndense = (csf_count*(csf_count + 1))/2

         if (csf_count <= iolpck) then
!-----------------------------------------------------------------------
!
! (2.1) - LAPACK    Dense, Memory,
!
!-------------------------------------------------------
            if (mpi_rank == 0) then
               write (6, *) &
                  "LAPACK routine DSPEVX selected for eigenvalue problem."
               write (24, *) &
                  "LAPACK routine DSPEVX selected for eigenvalue problem."
            end if

!   Allocate storage for the dense representation of the matrix
!   and initialize emt

            call alloc (hamiltonian_values, ndense, "EMT", "MANEIG")
            call dinit (ndense, 0.0_real64, hamiltonian_values, 1)

!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            call alloc (real_workspace, csf_count,"WORK", "MANEIG" )
            call alloc (matrix_partner_indices, csf_count, "IROW", "MANEIG")
            read (imcdf) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
            if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) stop &
               "maneig:1"

            do i = mpi_rank + 1, csf_count, mpi_process_count
               iofset = (i*(i - 1))/2
               read (imcdf) nelc, elsto, (real_workspace(ir),ir=1,nelc), (matrix_partner_indices(ir),ir=1,&
                  nelc)
               ! In the row-mode of the lower triangle,
               ! diagonal is the last one
               do ir = 1, nelc - 1
                  hamiltonian_values(iofset+matrix_partner_indices(ir)) = real_workspace(ir)
               end do
               hamiltonian_values(iofset+matrix_partner_indices(nelc)) = real_workspace(nelc) - hamiltonian_energy_shift

            end do

            call dalloc (real_workspace, "WORK", "MANEIG")
            call dalloc (matrix_partner_indices, "IROW", "MANEIG")

!   Find the eigenpairs
!
!    ivec() - serial numbers of eigenstates of the current block
!    iccmin() - serial numbers of eigenstates of all blocks.
!    nvecmn - minimum serial number of the eigenstates of the block
!    nvecmx - maximum .............
!    nvex - clear from def: NVECMX-NVECMN+1

            nvecmn = csf_count
            do i = 1, nvec
               nvecmn = min(nvecmn,ivec(i))
            end do
            nvex = nvecmx - nvecmn + 1
            call alloc (w, nvex, "W", "MANEIG")
            call alloc (z, csf_count*nvex,"Z", "MANEIG" )
            call alloc (real_workspace, csf_count*8,"WORK", "MANEIG" )
            call alloc (integer_workspace, csf_count*5,"IWORK", "MANEIG" )
! GG            CALL ALLOC (IFAIL, NVEX, 'IFAIL', 'MANEIG')
            call alloc (ifail, csf_count, "IFAIL", "MANEIG")
            call dspevx ("V", "I", "U", csf_count, hamiltonian_values, dummy, dummy, nvecmn, nvecmx&
               , abstol, m, w, z, csf_count, real_workspace, integer_workspace, ifail, info)
            if (info /= 0) stop "maneig: Failure in DSPEVX [LAPACK]"
            call dalloc (real_workspace, "WORK", "MANEIG")
            call dalloc (integer_workspace, "IWORK", "MANEIG")
            call dalloc (ifail, "IFAIL", "MANEIG")
            call dalloc (hamiltonian_values, "EMT", "MANEIG")

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            call alloc (shifted_eigenvalues, nvec,"EVAL", "MANEIG" )
            call alloc (csf_mixing_coefficients, csf_count*nvec, "EVEC", "MANEIG")

            do i = 1, nvec
               loc = ivec(i)
               shifted_eigenvalues(i) = w(loc - nvecmn + 1)
               iofset = csf_count*(i - 1)
               loc = csf_count*(loc - nvecmn)
               call dcopy (csf_count, z(loc + 1:loc+(csf_count)), 1, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)
            end do
            call dalloc (w, "W", "MANEIG")
            call dalloc (z, "Z", "MANEIG")

         else
!-----------------------------------------------------------------------

! (2.2) - DVDSON --- preparation work
!
!-------------------------------------------------------
           write (24,*)"DVDSON routine selected for eigenvalue problem;"

!   Sparse or dense matrix multiply? On disc or in core?

!GG            NBRKEV = (NCF + 1)*(NCF + 1)/3       ! Normal
            !NBRKEV =  1                      ! To enforce DENSE
            !NBRKEV = (NCF*(NCF+1)) / 2 + 1   ! To enforde SPARSE
!--------------------------------------------------------------
! Uncomment out the following to force sparse storage
!
!
!--------------------------------------------------------------
               sparse = .true.
               nstore = stored_matrix_element_count + stored_matrix_element_count/2 + (csf_count + 1)/2

            call alloc (matrix_diagonal, csf_count, "DIAG", "MANEIG")

            if (nstore > nincor) then
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal
!
!-------------------------------------------------------
               write (24, *) " matrix stored on disc;"

!   Disk storage; necessarily sparse; one column of the matrix in
!   memory

               ldisc = .true.
               sparse = .true.
               imv = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block

               read (imcdf) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
               if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) stop &
                  "maneig:2"

               do i = mpi_rank + 1, csf_count, mpi_process_count
                  read (imcdf) nelc, elsto, (dummy,ir=2,nelc), diatmp
                  matrix_diagonal(i) = diatmp - hamiltonian_energy_shift
               end do

            else
!-----------------------------------------------------------------------
!
! (2.2.2) - DVDSON --- Memory, load all
!
!-------------------------------------------------------
!
!   Core storage; load matrix into memory
!
               ldisc = .false.
               if (sparse) then
!-----------------------------------------------------------------------

! (2.2.2.1) - DVDSON --- Memory, load all, sparse
!
!-------------------------------------------------------
                  if (mpi_rank == 0) write (24, *) &
                     " matrix stored in sparse representation in core;"

                  imv = 2
                  write (6, *) "nelmnt = ", stored_matrix_element_count
                  call alloc (hamiltonian_values, stored_matrix_element_count, "EMT", "MANEIG")
                  call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "MANEIG")
                  call alloc (stored_segment_end_offsets, 0, csf_count , "IENDC", "MANEIG")
                  iofset = 0
                  stored_segment_end_offsets(0) = 0

                  read (imcdf) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                  if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                     stop "maneig:3"
                  do i = mpi_rank + 1, csf_count, mpi_process_count
                     read (imcdf) nelc, elsto, (hamiltonian_values(ir+iofset),ir=1,nelc), (&
                        matrix_partner_indices(ir + iofset),ir=1,nelc)
                     hamiltonian_values(nelc+iofset) = hamiltonian_values(nelc+iofset) - hamiltonian_energy_shift
                     matrix_diagonal(i) = hamiltonian_values(nelc+iofset)
                     iofset = iofset + nelc
                     stored_segment_end_offsets(i) = iofset
!                     WRITE (31 + MYID, *) I, IENDC(I), DIAG(I)

                  end do
               else
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  write (24, *) &
                     " matrix stored in full representation in core;"

                  imv = 3

! Find NDENSE_L, the number of elements on the node (dense form)

                  ndense_l = 0
                  do i = mpi_rank + 1, csf_count, mpi_process_count
                     ndense_l = ndense_l + i
                  end do

                  call alloc (hamiltonian_values, ndense_l, "EMT", "MANEIG")
                  call dinit (ndense_l, 0.0_real64, hamiltonian_values, 1)
                  call alloc (real_workspace, csf_count, "WORK", "MANEIG")
                  call alloc (matrix_partner_indices, csf_count, "IROW", "MANEIG")

                  read (imcdf) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                  if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                     stop "maneig:4"

                  iofset = 0
                  do i = mpi_rank + 1, csf_count, mpi_process_count
                     read (imcdf) nelc, elsto, (real_workspace(ir),ir=1,nelc), (matrix_partner_indices(ir),&
                        ir=1,nelc)
                     real_workspace(nelc) = real_workspace(nelc) - hamiltonian_energy_shift
                     matrix_diagonal(i) = real_workspace(nelc)
                     do ir = 1, nelc
                        hamiltonian_values(iofset+matrix_partner_indices(ir)) = real_workspace(ir)
                     end do
                     iofset = iofset + i
                  end do
                  call dalloc (real_workspace, "WORK", "MANEIG")
                  call dalloc (matrix_partner_indices, "IROW", "MANEIG")

               end if
!               ...Memory mode - sparse or dense
!-----------------------------------------------------------------------
!  (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
!-----------------------------------------------------------------------
            end if
!            ...Disk or Memory
!-----------------------------------------------------------------------
!  (2.2.3e)      *** E n d   o f   D V D S O N
!-----------------------------------------------------------------------
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained
!
            subspace_limit = min(csf_count,2*nvecmx + 60)
!           lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
            real_workspace_size = 2*csf_count*subspace_limit + subspace_limit*subspace_limit*2 + 11*subspace_limit + nvecmx
            call alloc (real_workspace, real_workspace_size, "WORK", "MANEIG")
            integer_workspace_size = 6*subspace_limit + nvecmx
            call alloc (integer_workspace, integer_workspace_size, "IWORK", "MANEIG")
!*changed by Misha 02/12/97
            eigenvalue_tolerance = 1.0e-17_real64
            coefficient_tolerance = 1.0e-9_real64
            residual_tolerance = 1.0e-9_real64
            orthogonality_tolerance = max(1e-9_real64,residual_tolerance)
! end of changes

!            maxitr = MAX (nvecmx*100, ncf/10)
            iteration_limit = max(nvecmx*200,csf_count/10)
            !maxitr = MIN (nvect*100, ncf)  ! FROM RSCFVU !!!
            call alloc (selected_pair_indices, subspace_limit,"JWORK", "MANEIG" )

            call alloc (shifted_eigenvalues, nvecmx, "EVAL", "MANEIG")
            call alloc (csf_mixing_coefficients, csf_count*nvecmx, "EVEC", "MANEIG")

            dmungo = 10.0e99_real64
            call dinit (nvecmx, dmungo, shifted_eigenvalues, 1)

!   Compute the eigenpairs in each block

            nvex = nvecmx
            if (ldisc) then
               expansion_block_size = nvex
            else
               expansion_block_size = 1
            end if
            nend = csf_count*nvex

            lowest_pair_index = 1
            highest_pair_index = nvex
            initial_vector_count = nvex
!************************************************************************
!
!   Call Davidson eigensolver
!
            select case (imv)
            case (1)
!******************** sparse and matrix on disk **********************
               write (6, *) " Sparse - Disk,  iniestsd"
               call posfile (0, imcdf, nposition)! was within iniestsd before
               call iniestsd (2000, csf_count, mpi_rank, mpi_process_count, initial_vector_count, real_workspace, imcdf, &
                  & hamiltonian_energy_shift)

               !NIV = 0   ! Why equal 0 ???
               !WRITE (6, *) ' Calling gdvd(spodmv,...'
               call gdvd (spodmv, csf_count, subspace_limit, matrix_diagonal, lowest_pair_index, highest_pair_index, &
                  & selected_pair_indices, initial_vector_count, &
                  expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                     & orthogonality_tolerance, iteration_limit, real_workspace, real_workspace_size, &
                  integer_workspace, integer_workspace_size, solve_high_end, iteration_count, matrix_vector_product_count, ierr)

            case (2)
!******************** sparse and matrix in memory ********************
               write (6, *) " Sparse - Memory,  iniest2"
!            CALL INIEST2 (1000, NCF,NIV,WORK,EMT,IENDC,IROW)
               call iniest2 (2000, csf_count, initial_vector_count, real_workspace, hamiltonian_values, &
                  & stored_segment_end_offsets, matrix_partner_indices)
               !WRITE (*, *) NCF, NIV, (WORK(I),I=NCF*NIV + 1,NCF*NIV + NIV)
               !WRITE (*, *) LIM, ILOW, IHIGH, MBLOCK, MAXITR, LWORK, LIWORK
               !WRITE (*, *) IERR
               call gdvd (spicmv2, csf_count, subspace_limit, matrix_diagonal, lowest_pair_index, highest_pair_index, &
                  & selected_pair_indices, initial_vector_count, &
                  expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                     & orthogonality_tolerance, iteration_limit, real_workspace, real_workspace_size, &
                  integer_workspace, integer_workspace_size, solve_high_end, iteration_count, matrix_vector_product_count, ierr)
               !WRITE (*, *) 'after gdvd...'
               !WRITE (*, *) NCF, NIV, (WORK(I),I=NCF*NIV + 1,NCF*NIV + NIV)
               !WRITE (*, *) LIM, ILOW, IHIGH, MBLOCK, MAXITR, LWORK, LIWORK
               !WRITE (*, *) HIEND, NLOOPS, NMV, IERR

               call dalloc (hamiltonian_values, "EMT", "MANEIG")
               call dalloc (matrix_partner_indices, "IROW", "MANEIG")
               call dalloc (stored_segment_end_offsets, "IENDC", "MANEIG")

            case (3)
!*************************** dense and in memory **********************
               write (6, *) " Dense - Memory,  iniestdm"
!            CALL INIESTDM (1000,NCF,NIV,WORK,EMT)
               call iniestdm (2000, csf_count, initial_vector_count, real_workspace, hamiltonian_values)
               call gdvd (dnicmv, csf_count, subspace_limit, matrix_diagonal, lowest_pair_index, highest_pair_index, &
                  & selected_pair_indices, initial_vector_count, &
                  expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                     & orthogonality_tolerance, iteration_limit, real_workspace, real_workspace_size, &
                  integer_workspace, integer_workspace_size, solve_high_end, iteration_count, matrix_vector_product_count, ierr)
               call dalloc (hamiltonian_values, "EMT", "MANEIG")
            end select
!************************************************************************
            call dalloc (matrix_diagonal, "DIAG", "MANEIG")
            call dalloc (integer_workspace, "IWORK", "MANEIG")
            call dalloc (selected_pair_indices, "JWORK", "MANEIG")

            write (24, *) " ", iteration_count, " iterations;"
            write (24, *) " ", matrix_vector_product_count, " matrix-vector multiplies."

            if (ierr /= 0) then
               write (stderr_unit, *) "MANEIG: Returned from DVDSON with"
               write (stderr_unit, *) " IERR = ", ierr, "."
               stop "maneig: DVDSON wrong"
            end if

!   Put the eigenpairs in order, overwriting as necessary

            call dcopy (nvex, real_workspace(nend+1:nend+(nvex)), 1, shifted_eigenvalues, 1)
            call dcopy (csf_count*nvex, real_workspace(1:csf_count*nvex), 1, csf_mixing_coefficients, 1)
            call dalloc (real_workspace, "WORK", "MANEIG")

!   Rearrange and reallocate storage for the eigenpairs
!   as necessary

            if (nvec < nvecmx) then
               call alloc (integer_workspace, nvecmx, "IWORK", "MANEIG")
               do i = 1, nvecmx
                  integer_workspace(i) = i
               end do
               do i = 1, nvec
                  iofset = ivec(i)
                  loc = integer_workspace(i)
                  if (iofset == loc) cycle
                  call dswap (1, shifted_eigenvalues(iofset:(iofset)+((1)-1)), 1, shifted_eigenvalues(i:(i)+((1)-1)), 1)
                  integer_workspace(i) = integer_workspace(iofset)
                  integer_workspace(iofset) = loc
                  iofset = csf_count*(iofset - 1)
                  loc = csf_count*(i - 1)
                  call dswap (csf_count, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1, &
                     & csf_mixing_coefficients(loc+1:loc+(csf_count)), 1)
               end do
               call dalloc (integer_workspace, "IWORK", "MANEIG")
               call ralloc (shifted_eigenvalues,  nvec, "EVAL", "MANEIG")
               call ralloc (csf_mixing_coefficients, csf_count*nvec, "EVEC", "MANEIG" )

            end if

         end if
! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

      end if
! (3e)               *** E N D   O F    A L L

!--------------------------------------------------------------------
! Only the following quantities are needed after this routine is
! finished:
!   eval(), evec(), iatjpo, iaspar
!--------------------------------------------------------------------
!
!   Clean up eigenvectors; determine their J/P values
!
      do j = 1, nvec

!         Find the dominant component of each eigenvector

         iofset = (j - 1)*csf_count

         amax = 0.0_real64
         do i = 1, csf_count
            wa = abs(csf_mixing_coefficients(i+iofset))
            if (wa <= amax) cycle
            amax = wa
            ia = i
         end do

!          Find the angular momentum and parity of the dominant component

         iatjpo = itjpo(ia)
         iaspar = ispar(ia)

!          Change sign of eigenvactor if dominant component is negative

         if (csf_mixing_coefficients(ia+iofset) >= 0.0_real64) cycle
         csf_mixing_coefficients(1+iofset:csf_count+iofset) = -csf_mixing_coefficients(1+iofset:csf_count+iofset)
      end do

      return
      end subroutine maneig
