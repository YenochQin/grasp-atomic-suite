module diagonalize_hamiltonian_block_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: diagonalize_hamiltonian_block
contains

    !***********************************************************************
    !                                                                      *
    subroutine diagonalize_hamiltonian_block(iatjpo, iaspar, nelmnt_a)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: summary_output_unit
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
        !
        !cjb comments [re-introduced] from GRASP-2013
        !   iatjpo - Output. (2j+1) of the dominant coefficients
        !   iaspar - Output, Parity (1 or -1) of the dominant coefficients
        !   nelmnt - Number of non-zero matrix elements on the current process *
        !            stored in the common block                                *
        !   nelmnt_a - The maximum number of non-zero matrix elements for all  *
        !             nodes                                                    *
        !                                                                      *
        !     value  meaning  of  method used IV                               *
        !       1    LAPACK                  =>  ncf < IOLPCK                  *
        !       2    Davidson, Dense-Memory  =>  dense requires less memory    *
        !       3    Davidson, Sparse-Memory =>  sparse requires less memory   *
        !       4    Davidson, Sparse-Disk   =>  memory requirement too large  *
        !                                        for at least one  process     *
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
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use array_memory_management, only: alloc, dalloc, ralloc
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use rci_matrix_file_position, only: nposition
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
            matrix_partner_indices
        use orbital_csf_basis, only: csf_count
        use selected_level_indices, only: active_level_indices, active_level_count, highest_active_level_index
        use ci_restart_units, only: restart_matrix_unit
        use standard_io_units, only: stderr_unit
        use mpi_state, only: mpi_error_code, mpi_process_count, mpi_rank

        ! Procedure interfaces
        use multiply_sparse_hamiltonian_interface, only: multiply_sparse_hamiltonian

        ! Module dependencies
        use sparse_symmetric_matvec, only: multiply_sparse_symmetric_matrix

        ! Procedure interfaces
        use multiply_restart_hamiltonian_interface, only: multiply_restart_hamiltonian

        ! Module dependencies
        use file_positioning, only: position_calculation_file
        use real_vector_initialization, only: initialize_real_vector

        ! Procedure interfaces
        use initialize_sparse_eigenvectors_interface, only: initialize_sparse_eigenvectors

        ! Module dependencies
        use davidson_eigensolver, only: dvd_operator, solve_davidson_eigenpairs
        use sparse_eigensolver_initial_basis, only: estimate_sparse_initial_basis

        ! Procedure interfaces
        use initialize_dense_eigenvectors_interface, only: initialize_dense_eigenvectors

        ! Module dependencies
        use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
        use mpi_sparse_matrix_product, only: multiply_mpi_sparse_matrix

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dcopy, dspevx, dswap

        ! Module dependencies
        use mpi_initial_eigenvectors, only: estimate_mpi_initial_eigenvectors
        use mpi_utilities, only: allreduce_mpi_real_sum
        use mpi_utilities, only: finalize_mpi_runtime

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dlamch
        implicit none (type, external)

        ! Arguments
        !-----------------------------------------------
        !  E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        integer, intent(out) :: iatjpo
        integer, intent(out) :: iaspar
        integer(int64)        :: nelmnt_a

        ! Local constants
        integer, parameter      :: dense_initial_csf_limit = 2000
        ! GG      REAL(DOUBLE), PARAMETER :: ABSTOL = 1.0D-10
        !cjb  NINCOR
        !cjb  INTEGER, PARAMETER      :: NINCOR = 1         ! To enforce DISK
        integer, parameter      :: in_memory_matrix_element_limit = 268435456 ! = 2 GB or  memory

        ! Local declarations
        !cff   ... variable for deciding between sparse or dense
        !     Memory need:
        !     sparse : nstore_s = nelmnt + (isize*(nelmnt + ncf+1))/8
        !     dense  : nstore_d = (ncf*ncf+1))/(2*nprocs)    (on average)
        !     Use Davidson (disk)  if nstore_s > NINCOR
        !     Use Davidson (memory, sparse) if nstore_s < nstore_d
        !     Use Davidson (memory, dense) if nstore_d < nstore_s
        !-----------------------------------------------
        integer(int64) :: nstore_s, nstore_d
        integer :: row_count, i, ndense, file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count, iofset, nelc, ir, &
            nvecmn, nvex, m, info, loc, imv, ndense_l, subspace_limit, real_workspace_size, integer_workspace_size, &
            iteration_limit, expansion_block_size, nend, lowest_pair_index, highest_pair_index, initial_vector_count, &
            iteration_count, matrix_vector_product_count, j, ia, idummy
        real(real64) :: elsto,  dummy, &
            diatmp, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, dmungo, amax, &
            & wa, abstol

        logical :: solve_high_end, ldisc, sparse
        real(real64), pointer :: w(:), z(:), real_workspace(:), matrix_diagonal(:)
        integer, pointer :: integer_workspace(:), ifail(:), selected_pair_indices(:)
        !-----------------------------------------------
        !
        !-----------------------------------------------------------------------
        abstol = 2*dlamch("S")
        if (mpi_rank == 0) write (output_unit, *) "Calling maneig..."

        !  (nrows+1) is the number of records of the present block's .res file

        row_count = (csf_count - mpi_rank - 1 + mpi_process_count)/mpi_process_count
        if (csf_count < mpi_process_count) row_count = csf_count/(mpi_rank + 1)
        !CALL posfile (1, restart_matrix_unit, nrows+1)
        call position_calculation_file (0, restart_matrix_unit, nposition)

        if (csf_count == 1) then
            !-----------------------------------------------------------------------
            !
            ! (1) - Trivial    ncf = 1
            !
            !-------------------------------------------------------
            if (mpi_rank == 0) write (summary_output_unit, *) "Trivial eigenvalue problem."

            !   Matrix of order 1: the trivial case; we assume that the value
            !   of EAV is available

            call alloc (shifted_eigenvalues, 1,"EVAL", "MANEIG" )
            call alloc (csf_mixing_coefficients, 1, "EVECO", "MANEIG")
            shifted_eigenvalues(1) = 0.0_real64
            csf_mixing_coefficients(1) = 1.0_real64

            ! Still read through the .res file
            !GG
            !GG Gediminas NIST 2005.11.03
            !GG         READ (restart_matrix_unit)
            do i = 1, row_count + 1
                read (restart_matrix_unit)
            end do

        else                                       !if-2
            !-----------------------------------------------------------------------
            !
            ! (2) - Non trivial
            !
            !-------------------------------------------------------
            if(mpi_rank==0) write(*,*) "ncf=",csf_count, " iolpck=", dense_initial_csf_limit
            if (csf_count <= dense_initial_csf_limit) then
                !-----------------------------------------------------------------------
                !
                ! (2.1) - LAPACK    Dense, Memory,
                !
                !-------------------------------------------------------
                if (mpi_rank == 0) then
                    write (output_unit, *) &
                        "LAPACK routine DSPEVX selected for eigenvalue problem."
                    write (summary_output_unit, *) &
                        "LAPACK routine DSPEVX selected for eigenvalue problem."
                end if

                !   Allocate storage for the dense representation of the matrix
                !   and initialize emt

                ndense = (csf_count*(csf_count + 1))/2
                call alloc (hamiltonian_values, ndense, "EMT", "MANEIG")
                call initialize_real_vector (ndense, 0.0_real64, hamiltonian_values, 1)

                !   Read the matrix into position from the disc file; it's already
                !   been properly positioned.

                call alloc (real_workspace, csf_count,"WORK", "MANEIG" )
                call alloc (matrix_partner_indices, csf_count, "IROW", "MANEIG")
                read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) stop &
                    "maneig:1"

                do i = mpi_rank + 1, csf_count, mpi_process_count
                    iofset = (i*(i - 1))/2
                    read (restart_matrix_unit) nelc, elsto, (real_workspace(ir),ir=1,nelc), (matrix_partner_indices(ir),ir=1,&
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

                ! Let each node have a complete copy of EMT

                call allreduce_mpi_real_sum (hamiltonian_values, ndense)

                !   Find the eigenpairs
                !
                !    ivec() - serial numbers of eigenstates of the current block
                !    iccmin() - serial numbers of eigenstates of all blocks.
                !    nvecmn - minimum serial number of the eigenstates of the block
                !    nvecmx - maximum .............
                !    nvex - clear from def: NVECMX-NVECMN+1

                nvecmn = csf_count
                do i = 1, active_level_count
                    nvecmn = min(nvecmn,active_level_indices(i))
                end do
                nvex = highest_active_level_index - nvecmn + 1
                call alloc (w, nvex, "W", "MANEIG")
                call alloc (z, csf_count*nvex,"Z", "MANEIG" )
                call alloc (real_workspace, csf_count*8,"WORK", "MANEIG" )
                call alloc (integer_workspace, csf_count*5,"IWORK", "MANEIG" )
                ! GG            CALL ALLOC (IFAIL, NVEX, 'IFAIL', 'MANEIG')
                call alloc (ifail, csf_count, "IFAIL", "MANEIG")
                call dspevx ("V", "I", "U", csf_count, hamiltonian_values, dummy, dummy, nvecmn, highest_active_level_index&
                    , abstol, m, w, z, csf_count, real_workspace, integer_workspace, ifail, info)
                if (info /= 0)                                             &
                    call finalize_mpi_runtime("maneig: Failure in DSPEVX [LAPACK]",mpi_rank)
                call dalloc (real_workspace, "WORK", "MANEIG")
                call dalloc (integer_workspace, "IWORK", "MANEIG")
                call dalloc (ifail, "IFAIL", "MANEIG")
                call dalloc (hamiltonian_values, "EMT", "MANEIG")

                !   Store the eigenpairs in their proper positions EVAL() and EVEC()

                call alloc (shifted_eigenvalues, active_level_count,"EVAL", "MANEIG")
                call alloc (csf_mixing_coefficients, csf_count*active_level_count, "EVEC", "MANEIG")

                do i = 1, active_level_count
                    loc = active_level_indices(i)
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
                if (mpi_rank == 0) write (summary_output_unit,*)"DVDSON routine selected" &
                    //" for eigenvalue problem;"
                !--------------------------------------------------------------
                if (mpi_rank == 0) print *,"zou: nelmnt",stored_matrix_element_count
                !            print *, 'myid, nprocs,nincor', myid, nprocs, nincor

                !           If on any process more memory is needed than NINCOR, use disc.
                !           If not, should the matrix be sparse or dense (distributed)?
                !           sparse requires = nelmnt + (isize*(nelmnt+ncf+1))/8
                !           dense  storage = (ncf*(ncf+1)/2/nproc - on average
                !           (because of integer arithmetic nproc-1 is added to numerator

                !           This is a comparion of matrix elements only.
                !
                nstore_s = nelmnt_a
                nstore_d = ((csf_count+mpi_process_count-1)/mpi_process_count)
                nstore_d = nstore_d*((csf_count+1)/2)
                if (mpi_rank == 0) write(*,*) "nstore_s=", nstore_s,   &
                    "nstore_d=", nstore_d, "nincor=", in_memory_matrix_element_limit
                ldisc = .false.
                if ( nstore_s > in_memory_matrix_element_limit) then
                    !              .. we need to use disk
                    ldisc =.true.
                    sparse=.true.
                else

                    if (nstore_s < nstore_d ) then
                        sparse = .true.
                    else
                        sparse = .false.
                    end if
                end if

                call alloc (matrix_diagonal,csf_count,"DIAG", "MANEIG")
                call initialize_real_vector (csf_count, 0.0_real64, matrix_diagonal, 1)
                if (ldisc) then
                    !-----------------------------------------------------------------------
                    !
                    ! (2.2.1) - DVDSON --- Disk, load diagonal
                    !
                    !-------------------------------------------------------
                    if (mpi_rank == 0)  then
                        write (*, *) " matrix stored on disc;"
                        write (summary_output_unit, *) " matrix stored on disc;"
                    end if

                    !   Disk storage; necessarily sparse; one column of the matrix in
                    !   memory

                    imv = 1

                    !   Load diagonal - Each node will have the same, complete copy
                    !   after this if block
                    print *, mpi_rank, "Reading diagonal, unit imcdf =", restart_matrix_unit
                    read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                    if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                        & stop &
                        "maneig:2"

                    do i = mpi_rank + 1, csf_count, mpi_process_count
                        read (restart_matrix_unit) nelc, elsto, (dummy,ir=2,nelc), diatmp, &
                            (idummy, ir=1,nelc)
                        matrix_diagonal(i) = diatmp - hamiltonian_energy_shift
                    end do
                    print *, " finished"
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
                        if (mpi_rank == 0) then
                            write (summary_output_unit, *) " sparse matrix stored in memory"
                            write ( *, *) " sparse matrix stored in memory"
                        end if

                        imv = 2
                        write (output_unit, *) "nelmnt = ", stored_matrix_element_count
                        call alloc (hamiltonian_values, stored_matrix_element_count, "EMT", "MANEIG")
                        call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "MANEIG")
                        call alloc (stored_segment_end_offsets, 0, csf_count , "IENDC", "MANEIG")
                        iofset = 0
                        stored_segment_end_offsets(0) = 0

                        read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                        if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                            stop "maneig:3"
                        do i = mpi_rank + 1, csf_count, mpi_process_count
                            read (restart_matrix_unit) nelc, elsto, (hamiltonian_values(ir+iofset),ir=1,nelc), (&
                                matrix_partner_indices(ir + iofset),ir=1,nelc)
                            hamiltonian_values(nelc+iofset) = hamiltonian_values(nelc+iofset) - hamiltonian_energy_shift
                            matrix_diagonal(i) = hamiltonian_values(nelc+iofset)
                            iofset = iofset + nelc
                            stored_segment_end_offsets(i) = iofset
                        end do
                    else
                        !-----------------------------------------------------------------------
                        !
                        ! (2.2.2.2) - DVDSON --- Memory, load all, dense
                        !
                        !-------------------------------------------------------
                        if (mpi_rank == 0) then
                            write (summary_output_unit,*) " full matrix stored in in memory"
                            write ( *,*) " full matrix stored in in memory"
                        end if

                        imv = 3

                        ! Find NDENSE_L, the number of elements on the node (dense form)

                        ndense_l = 0
                        do i = mpi_rank + 1, csf_count, mpi_process_count
                            ndense_l = ndense_l + i
                        end do

                        call alloc (hamiltonian_values, ndense_l, "EMT", "MANEIG")
                        call initialize_real_vector (ndense_l, 0.0_real64, hamiltonian_values, 1)
                        call alloc (real_workspace, csf_count, "WORK", "MANEIG")
                        call alloc (matrix_partner_indices, csf_count, "IROW", "MANEIG")

                        read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                        if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                            stop "maneig:4"

                        iofset = 0
                        do i = mpi_rank + 1, csf_count, mpi_process_count
                            read (restart_matrix_unit) nelc, elsto, (real_workspace(ir),ir=1,nelc), (matrix_partner_indices(ir),&
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

                ! Make diagonals global, no matter it is disk or memory mode

                call allreduce_mpi_real_sum (matrix_diagonal, csf_count)
                !
                !   Allocate storage for workspace; see the header of DVDSON for
                !   the expression below; the value of LIM can be reduced to NVECMX
                !   plus a smaller number if storage is severely constrained
                !
                subspace_limit = min(csf_count,2*highest_active_level_index + 60)
                !           lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
                real_workspace_size = 2*csf_count*subspace_limit + subspace_limit*subspace_limit*2 + 11*subspace_limit + &
                    highest_active_level_index
                call alloc (real_workspace, real_workspace_size, "WORK", "MANEIG")
                real_workspace(1:real_workspace_size) = 0.0_real64
                integer_workspace_size = 6*subspace_limit + highest_active_level_index
                call alloc (integer_workspace, integer_workspace_size, "IWORK", "MANEIG")
                !*changed by Misha 02/12/97
                eigenvalue_tolerance = 1.0e-17_real64
                !GG            CRITC = 1.0D-08
                !GG            CRITR = 1.0D-08

                coefficient_tolerance = 1.0e-9_real64
                residual_tolerance = 1.0e-9_real64
                ! end of changes

                !            maxitr = MAX (nvecmx*100, ncf/10)
                iteration_limit = max(highest_active_level_index*100,csf_count/10)
                !maxitr = MIN (nvect*100, ncf)  ! FROM RSCFVU !!!
                call alloc (selected_pair_indices, subspace_limit,"JWORK", "MANEIG" )

                call alloc (shifted_eigenvalues, highest_active_level_index, "EVAL", "MANEIG")
                call alloc (csf_mixing_coefficients, csf_count*highest_active_level_index, "EVEC", "MANEIG")

                dmungo = 10.0e99_real64
                call initialize_real_vector (highest_active_level_index, dmungo, shifted_eigenvalues, 1)

                !   Compute the eigenpairs in each block

                nvex = highest_active_level_index
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
                    if (mpi_rank == 0) print *, " Sparse - Disk, iniestsd"
                    call position_calculation_file (0, restart_matrix_unit, nposition)! was within iniestsd before
                    call initialize_sparse_eigenvectors (dense_initial_csf_limit, csf_count, mpi_rank, mpi_process_count, &
                        & initial_vector_count, &
                        & real_workspace, &
                        & restart_matrix_unit, &
                        hamiltonian_energy_shift)
                    print *, "Returned from iniestsd "
                    if (csf_count> dense_initial_csf_limit) then
                        print *, "Calling solve_davidson_eigenpairs"
                        call solve_davidson_eigenpairs (multiply_restart_hamiltonian,csf_count,subspace_limit,matrix_diagonal, &
                            & lowest_pair_index, &
                            & highest_pair_index,            &
                            selected_pair_indices,initial_vector_count,expansion_block_size,eigenvalue_tolerance, &
                            & coefficient_tolerance, &
                            & residual_tolerance,iteration_limit,   &
                            real_workspace,real_workspace_size,integer_workspace,integer_workspace_size,solve_high_end, &
                            & iteration_count, &
                            matrix_vector_product_count,mpi_error_code)
                    end if

                case (2)
                    !******************** sparse and matrix in memory ********************
                    if (mpi_rank == 0) print *, " Sparse - Memory, iniestmpi"
                    call estimate_mpi_initial_eigenvectors (dense_initial_csf_limit, &
                        & csf_count,initial_vector_count,real_workspace, &
                        & hamiltonian_values, &
                        & stored_segment_end_offsets,matrix_partner_indices)
                    if(csf_count> dense_initial_csf_limit) then
                        call solve_davidson_eigenpairs (multiply_mpi_sparse_matrix,csf_count,subspace_limit,matrix_diagonal, &
                            & lowest_pair_index,highest_pair_index,        &
                            selected_pair_indices,initial_vector_count,expansion_block_size,eigenvalue_tolerance, &
                            & coefficient_tolerance, &
                            & residual_tolerance,iteration_limit,   &
                            real_workspace,real_workspace_size,integer_workspace,integer_workspace_size,solve_high_end, &
                            & iteration_count, &
                            matrix_vector_product_count,mpi_error_code)
                    end if

                    call dalloc (hamiltonian_values, "EMT", "MANEIG")
                    call dalloc (matrix_partner_indices, "IROW", "MANEIG")
                    call dalloc (stored_segment_end_offsets, "IENDC", "MANEIG")

                case (3)
                    !*************************** dense and in memory **********************
                    if (mpi_rank == 0) print *, " Dense - Memory, iniestdm"
                    call initialize_dense_eigenvectors &
                        & (dense_initial_csf_limit,csf_count,initial_vector_count,real_workspace,hamiltonian_values)
                    if (csf_count> dense_initial_csf_limit) then
                        call solve_davidson_eigenpairs (multiply_sparse_hamiltonian,csf_count,subspace_limit,matrix_diagonal, &
                            & lowest_pair_index, &
                            & highest_pair_index,            &
                            selected_pair_indices,initial_vector_count,expansion_block_size,eigenvalue_tolerance, &
                            & coefficient_tolerance, &
                            & residual_tolerance,iteration_limit, &
                            real_workspace,real_workspace_size,integer_workspace,integer_workspace_size,solve_high_end, &
                            & iteration_count, &
                            matrix_vector_product_count,mpi_error_code)
                    end if
                    call dalloc (hamiltonian_values, "EMT", "MANEIG")
                end select
                !************************************************************************
                call dalloc (matrix_diagonal, "DIAG", "MANEIG")
                call dalloc (integer_workspace, "IWORK", "MANEIG")
                call dalloc (selected_pair_indices, "JWORK", "MANEIG")
                if (mpi_rank == 0) then
                    write (summary_output_unit, *) " ", iteration_count, " iterations;"
                    write (summary_output_unit, *) " ", matrix_vector_product_count, " matrix-vector multiplies."
                end if
                if (mpi_error_code /= 0) then
                    write (stderr_unit, *) "MANEIG: Returned from DVDSON with"
                    write (stderr_unit, *) " IERR = ", mpi_error_code, "."
                    stop "maneig: DVDSON wrong"
                end if

                !   Put the eigenpairs in order, overwriting as necessary

                call dcopy (nvex, real_workspace(nend+1:nend+(nvex)), 1, shifted_eigenvalues, 1)
                call dcopy (csf_count*nvex, real_workspace(1:csf_count*nvex), 1, csf_mixing_coefficients, 1)
                call dalloc (real_workspace, "WORK", "MANEIG")

                !   Rearrange and reallocate storage for the eigenpairs
                !   as necessary

                if (active_level_count < highest_active_level_index) then
                    call alloc (integer_workspace, highest_active_level_index, "IWORK", "MANEIG")
                    do i = 1, highest_active_level_index
                        integer_workspace(i) = i
                    end do
                    do i = 1, active_level_count
                        iofset = active_level_indices(i)
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
                    call ralloc (shifted_eigenvalues,  active_level_count, "EVAL", "MANEIG")
                    call ralloc (csf_mixing_coefficients, csf_count*active_level_count, "EVEC", "MANEIG" )

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
        do j = 1, active_level_count

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

            iatjpo = csf_total_j_code(ia)
            iaspar = csf_parity_code(ia)

            !          Change sign of eigenvactor if dominant component is negative

            if (csf_mixing_coefficients(ia+iofset) >= 0.0_real64) cycle
            csf_mixing_coefficients(1+iofset:csf_count+iofset) = -csf_mixing_coefficients(1+iofset:csf_count+iofset)
        end do

        return
    end subroutine diagonalize_hamiltonian_block
end module diagonalize_hamiltonian_block_interface
