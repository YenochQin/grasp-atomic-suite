! Solve for the requested eigenpairs of a Hamiltonian symmetry block.
module diagonalize_hamiltonian_block_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: diagonalize_hamiltonian_block
contains

    !  This module  manages the  operation of the  eigensolvers and the
    !  storage of the eigenpairs.  There are two principal branches:
    !
    !     (1) Matrix of order 1: the trivial case
    !     (2) Matrix of order greater than 1: there are two branches
    !            (i) Matrices of order less than or equal to IOLPCK:
    !                eigenpairs are found using LAPACK SUBROUTINEs
    !           (ii) Matrices of order greater than IOLPCK: eigenpairs
    !                are found using DVDSON; this involves up to three
    !                steps:
    !                   (a) The matrix is analysed to determine its
    !                       block structure (only irreducibe matrices
    !                       are correctly treated by DVDSON)
    !                   (b) Eigenpairs are extracted for each block
    !                   (c) The appropriate eigenpairs are selected and
    !                       stored
    !                Different methods of storage and different
    !                versions of the matrix-vector multiply are used
    !                depending upon the order and density of the matrix
    !
    ! Written by Farid A. Parpia. Last revision: 27 Sep 1993
    !  Modified Misha Saparov                                  Feb 1997
    !           Charlotte F. Fischer                           May 1997
    !                Except for the disk version, all matrices have
    !                diagonals  shifted by EAV
    ! All arrays allocated here are de-allocated except pneval and pnevec *
    ! which will be de-allocated in matrix.
    ! Block Version  By Xinghong He. Last revision: 18 Jun 1998
    subroutine diagonalize_hamiltonian_block(iatjpo, iaspar)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use array_memory_management, only: alloc, dalloc
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use rci_matrix_file_position, only: nposition
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
            matrix_partner_indices
        use orbital_csf_basis, only: csf_count
        use selected_level_indices, only: active_level_indices, active_level_count, highest_active_level_index
        use ci_restart_units, only: restart_matrix_unit
        use standard_io_units, only: stderr_unit

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
        use davidson_eigensolver, only: davidson_matrix_vector_operator, solve_davidson_eigenpairs
        use sparse_eigensolver_initial_basis, only: estimate_sparse_initial_basis

        ! Procedure interfaces
        use initialize_dense_eigenvectors_interface, only: initialize_dense_eigenvectors

        ! Module dependencies
        use csf_quantum_queries, only: csf_total_j_code, csf_parity_code

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dcopy, dspevx, dswap
        use blas_lapack_interfaces, only: dlamch
        implicit none (type, external)

        ! Arguments
        !
        !  E x t e r n a l   F u n c t i o n s
        integer, intent(out) :: iatjpo
        integer, intent(out) :: iaspar

        !  Local constants
        integer, parameter      :: dense_initial_csf_limit = 2000
        ! NINCOR
        integer, parameter      :: in_memory_matrix_element_limit = 268435456 ! = 2 GB or  memory

        ! Local declarations
        integer(int64) :: matrix_offset, ndense_l
        integer(int64) :: nstore
        integer :: mpi_rank, mpi_process_count, row_count, i, ndense, file_csf_count, iccutdum, file_mpi_rank, &
            file_mpi_process_count, iofset, nelc, ir, nvecmn, nvex, m, info, loc, imv, subspace_limit, &
            real_workspace_size, integer_workspace_size, iteration_limit, expansion_block_size, nend, lowest_pair_index, &
            highest_pair_index, initial_vector_count, iteration_count, matrix_vector_product_count, ierr, j, ia
        real(real64) :: elsto,  dummy, &
            diatmp, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, dmungo, amax, &
            & wa, abstol

        logical :: solve_high_end, ldisc, sparse
        real(real64), pointer :: w(:), z(:), real_workspace(:), matrix_diagonal(:)
        integer, pointer :: integer_workspace(:), ifail(:), selected_pair_indices(:)
        abstol = 2*dlamch("S")
        mpi_rank = 0
        mpi_process_count = 1

        ! (nrows+1) is the number of records of the present block's .res file

        row_count = (csf_count - mpi_rank - 1 + mpi_process_count)/mpi_process_count
        if (csf_count < mpi_process_count) row_count = csf_count/(mpi_rank + 1)
        call position_calculation_file (0, restart_matrix_unit, nposition)

        if (csf_count == 1) then
            ! (1) - Trivial    ncf = 1
            write (summary_output_unit, *) "Trivial eigenvalue problem."

            ! Matrix of order 1: the trivial case; we assume that the value
            ! of EAV is available

            allocate(shifted_eigenvalues(1))
            allocate(csf_mixing_coefficients(1))
            shifted_eigenvalues(1) = 0.0_real64
            csf_mixing_coefficients(1) = 1.0_real64

            !  Still read through the .res file
            ! GG Gediminas NIST 2005.11.03
            do i = 1, row_count + 1
                read (restart_matrix_unit)
            end do

        else                                       !if-2
            ! (2) - Non trivial
            !
            !   Matrix of order greater than 1; how many elements in a triangle?
            if (csf_count <= dense_initial_csf_limit) then
                ndense = (csf_count*(csf_count + 1))/2
                ! (2.1) - LAPACK    Dense, Memory,
                if (mpi_rank == 0) then
                    write (output_unit, *) &
                        "LAPACK routine DSPEVX selected for eigenvalue problem."
                    write (summary_output_unit, *) &
                        "LAPACK routine DSPEVX selected for eigenvalue problem."
                end if

                ! Allocate storage for the dense representation of the matrix
                ! and initialize emt

                allocate(hamiltonian_values(ndense))
                call initialize_real_vector (ndense, 0.0_real64, hamiltonian_values, 1)

                ! Read the matrix into position from the disc file; it's already
                ! been properly positioned.

                call alloc (real_workspace, csf_count,"real_workspace", "diagonalize_hamiltonian_block" )
                allocate(matrix_partner_indices(csf_count))
                read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) stop &
                    "diagonalize_hamiltonian_block:1"

                do i = mpi_rank + 1, csf_count, mpi_process_count
                    matrix_offset = int(i,int64)*(i - 1)/2
                    read (restart_matrix_unit) nelc, elsto, (real_workspace(ir),ir=1,nelc), (matrix_partner_indices(ir),ir=1,&
                        nelc)
                    ! In the row-mode of the lower triangle,
                    ! diagonal is the last one
                    do ir = 1, nelc - 1
                        hamiltonian_values(matrix_offset+matrix_partner_indices(ir)) = real_workspace(ir)
                    end do
                    hamiltonian_values(matrix_offset+matrix_partner_indices(nelc)) = &
                        real_workspace(nelc) - hamiltonian_energy_shift

                end do

                call dalloc (real_workspace, "real_workspace", "diagonalize_hamiltonian_block")
                deallocate(matrix_partner_indices)

                ! Find the eigenpairs
                !
                !  ivec() - serial numbers of eigenstates of the current block
                !  iccmin() - serial numbers of eigenstates of all blocks.
                !  nvecmn - minimum serial number of the eigenstates of the block
                !  nvecmx - maximum .............
                !  nvex - clear from def: NVECMX-NVECMN+1

                nvecmn = csf_count
                do i = 1, active_level_count
                    nvecmn = min(nvecmn,active_level_indices(i))
                end do
                nvex = highest_active_level_index - nvecmn + 1
                call alloc (w, nvex, "W", "diagonalize_hamiltonian_block")
                call alloc (z, csf_count*nvex,"Z", "diagonalize_hamiltonian_block" )
                call alloc (real_workspace, csf_count*8,"real_workspace", "diagonalize_hamiltonian_block" )
                call alloc (integer_workspace, csf_count*5,"integer_workspace", "diagonalize_hamiltonian_block" )
                call alloc (ifail, csf_count, "IFAIL", "diagonalize_hamiltonian_block")
                call dspevx ("V", "I", "U", csf_count, hamiltonian_values, dummy, dummy, nvecmn, highest_active_level_index&
                    , abstol, m, w, z, csf_count, real_workspace, integer_workspace, ifail, info)
                if (info /= 0) stop "diagonalize_hamiltonian_block: Failure in DSPEVX [LAPACK]"
                call dalloc (real_workspace, "real_workspace", "diagonalize_hamiltonian_block")
                call dalloc (integer_workspace, "integer_workspace", "diagonalize_hamiltonian_block")
                call dalloc (ifail, "IFAIL", "diagonalize_hamiltonian_block")
                deallocate(hamiltonian_values)

                ! Store the eigenpairs in their proper positions EVAL() and EVEC()

                allocate(shifted_eigenvalues(active_level_count))
                allocate(csf_mixing_coefficients(int(csf_count,int64)*active_level_count))

                do i = 1, active_level_count
                    loc = active_level_indices(i)
                    shifted_eigenvalues(i) = w(loc - nvecmn + 1)
                    iofset = csf_count*(i - 1)
                    loc = csf_count*(loc - nvecmn)
                    call dcopy (csf_count, z(loc + 1:loc+(csf_count)), 1, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)
                end do
                call dalloc (w, "W", "diagonalize_hamiltonian_block")
                call dalloc (z, "Z", "diagonalize_hamiltonian_block")

            else

                ! (2.2) - DVDSON --- preparation work
                write (summary_output_unit,*)"solve_davidson_eigenpairs routine selected for eigenvalue problem;"

                ! Sparse or dense matrix multiply? On disc or in core?

                ! NBRKEV = (NCF + 1)*(NCF + 1)/3       ! Normal
                ! NBRKEV =  1                      ! To enforce DENSE
                ! NBRKEV = (NCF*(NCF+1)) / 2 + 1   ! To enforde SPARSE
                !
                !  Uncomment out the following to force sparse storage
                sparse = .true.
                nstore = stored_matrix_element_count + (stored_matrix_element_count + 1)/2 + int(csf_count,int64) + 1

                call alloc (matrix_diagonal, csf_count, "matrix_diagonal", "diagonalize_hamiltonian_block")

                if (nstore > in_memory_matrix_element_limit) then
                    ! (2.2.1) - DVDSON --- Disk, load diagonal
                    write (summary_output_unit, *) " matrix stored on disc;"

                    ! Disk storage; necessarily sparse; one column of the matrix in
                    ! memory

                    ldisc = .true.
                    sparse = .true.
                    imv = 1

                    ! Load diagonal - Each node will have the same, complete copy
                    ! after this if block

                    read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                    if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                        & stop &
                        "diagonalize_hamiltonian_block:2"

                    do i = mpi_rank + 1, csf_count, mpi_process_count
                        read (restart_matrix_unit) nelc, elsto, (dummy,ir=2,nelc), diatmp
                        matrix_diagonal(i) = diatmp - hamiltonian_energy_shift
                    end do

                else
                    ! (2.2.2) - DVDSON --- Memory, load all
                    !
                    !   Core storage; load matrix into memory
                    ldisc = .false.
                    if (sparse) then

                        ! (2.2.2.1) - DVDSON --- Memory, load all, sparse
                        if (mpi_rank == 0) write (summary_output_unit, *) &
                            " matrix stored in sparse representation in core;"

                        imv = 2
                        write (output_unit, *) "nelmnt = ", stored_matrix_element_count
                        allocate(hamiltonian_values(stored_matrix_element_count))
                        allocate(matrix_partner_indices(stored_matrix_element_count))
                        allocate(stored_segment_end_offsets(0:csf_count))
                        matrix_offset = 0_int64
                        stored_segment_end_offsets = 0_int64

                        read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                        if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                            stop "diagonalize_hamiltonian_block:3"
                        do i = mpi_rank + 1, csf_count, mpi_process_count
                            read (restart_matrix_unit) nelc, elsto, (hamiltonian_values(ir+matrix_offset),ir=1,nelc), (&
                                matrix_partner_indices(ir + matrix_offset),ir=1,nelc)
                            matrix_offset = matrix_offset + nelc
                            hamiltonian_values(matrix_offset) = hamiltonian_values(matrix_offset) - hamiltonian_energy_shift
                            matrix_diagonal(i) = hamiltonian_values(matrix_offset)
                            stored_segment_end_offsets(i) = matrix_offset

                        end do
                    else
                        ! (2.2.2.2) - DVDSON --- Memory, load all, dense
                        write (summary_output_unit, *) &
                            " matrix stored in full representation in core;"

                        imv = 3

                        ! Find NDENSE_L, the number of elements on the node (dense form)

                        ndense_l = 0
                        do i = mpi_rank + 1, csf_count, mpi_process_count
                            ndense_l = ndense_l + i
                        end do

                        allocate(hamiltonian_values(ndense_l))
                        hamiltonian_values = 0.0_real64
                        call alloc (real_workspace, csf_count, "real_workspace", "diagonalize_hamiltonian_block")
                        allocate(matrix_partner_indices(csf_count))

                        read (restart_matrix_unit) file_csf_count, iccutdum, file_mpi_rank, file_mpi_process_count
                        if (csf_count/=file_csf_count .or. mpi_rank/=file_mpi_rank .or. file_mpi_process_count/=mpi_process_count) &
                            stop "diagonalize_hamiltonian_block:4"

                        matrix_offset = 0_int64
                        do i = mpi_rank + 1, csf_count, mpi_process_count
                            read (restart_matrix_unit) nelc, elsto, (real_workspace(ir),ir=1,nelc), (matrix_partner_indices(ir),&
                                ir=1,nelc)
                            real_workspace(nelc) = real_workspace(nelc) - hamiltonian_energy_shift
                            matrix_diagonal(i) = real_workspace(nelc)
                            do ir = 1, nelc
                                hamiltonian_values(matrix_offset+matrix_partner_indices(ir)) = real_workspace(ir)
                            end do
                            matrix_offset = matrix_offset + i
                        end do
                        call dalloc (real_workspace, "real_workspace", "diagonalize_hamiltonian_block")
                        deallocate(matrix_partner_indices)

                    end if
                    !              ...Memory mode - sparse or dense
                    !
                    ! (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
                end if
                !           ...Disk or Memory
                !
                ! (2.2.3e)      *** E n d   o f   D V D S O N
                !
                !  Allocate storage for workspace; see the header of DVDSON for
                !  the expression below; the value of LIM can be reduced to NVECMX
                !  plus a smaller number if storage is severely constrained
                subspace_limit = min(csf_count,2*highest_active_level_index + 60)
                ! lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
                real_workspace_size = 2*csf_count*subspace_limit + subspace_limit*subspace_limit*2 + 11*subspace_limit + &
                    highest_active_level_index
                call alloc (real_workspace, real_workspace_size, "real_workspace", "diagonalize_hamiltonian_block")
                integer_workspace_size = 6*subspace_limit + highest_active_level_index
                call alloc (integer_workspace, integer_workspace_size, "integer_workspace", "diagonalize_hamiltonian_block")
                ! *changed by Misha 02/12/97
                eigenvalue_tolerance = 1.0e-17_real64
                coefficient_tolerance = 1.0e-9_real64
                residual_tolerance = 1.0e-9_real64
                ! end of changes

                ! maxitr = MAX (nvecmx*100, ncf/10)
                iteration_limit = max(highest_active_level_index*200,csf_count/10)
                call alloc (selected_pair_indices, subspace_limit,"selected_pair_indices", "diagonalize_hamiltonian_block" )

                allocate(shifted_eigenvalues(highest_active_level_index))
                allocate(csf_mixing_coefficients(int(csf_count,int64)*highest_active_level_index))

                dmungo = 10.0e99_real64
                call initialize_real_vector (highest_active_level_index, dmungo, shifted_eigenvalues, 1)

                ! Compute the eigenpairs in each block

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
                ! Call Davidson eigensolver
                select case (imv)
                case (1)
                    ! sparse and matrix on disk **********************
                    write (output_unit, *) " Sparse - Disk,  iniestsd"
                    call position_calculation_file (0, restart_matrix_unit, nposition)! was within iniestsd before
                    call initialize_sparse_eigenvectors (2000, csf_count, mpi_rank, mpi_process_count, initial_vector_count, &
                        & real_workspace, &
                        & restart_matrix_unit, hamiltonian_energy_shift)

                    ! NIV = 0   ! Why equal 0 ???
                    call solve_davidson_eigenpairs (multiply_restart_hamiltonian, csf_count, subspace_limit, matrix_diagonal, &
                        & lowest_pair_index, &
                        & highest_pair_index, &
                        & selected_pair_indices, initial_vector_count, &
                        expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                        & iteration_limit, real_workspace, real_workspace_size, &
                        integer_workspace, integer_workspace_size, solve_high_end, iteration_count, &
                        & matrix_vector_product_count, ierr)

                case (2)
                    ! sparse and matrix in memory ********************
                    write (output_unit, *) " Sparse - Memory,  iniest2"
                    call estimate_sparse_initial_basis (2000, csf_count, initial_vector_count, real_workspace, hamiltonian_values, &
                        & stored_segment_end_offsets, matrix_partner_indices)

                    call solve_davidson_eigenpairs (multiply_sparse_symmetric_matrix, csf_count, subspace_limit, &
                        & matrix_diagonal, lowest_pair_index, highest_pair_index, &
                        & selected_pair_indices, initial_vector_count, &
                        expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                        & iteration_limit, real_workspace, real_workspace_size, &
                        integer_workspace, integer_workspace_size, solve_high_end, iteration_count, &
                        & matrix_vector_product_count, ierr)

                    deallocate(hamiltonian_values, matrix_partner_indices, stored_segment_end_offsets)

                case (3)
                    ! dense and in memory **********************
                    write (output_unit, *) " Dense - Memory,  iniestdm"
                    call initialize_dense_eigenvectors (2000, csf_count, initial_vector_count, real_workspace, hamiltonian_values)
                    call solve_davidson_eigenpairs (multiply_sparse_hamiltonian, csf_count, subspace_limit, matrix_diagonal, &
                        & lowest_pair_index, &
                        & highest_pair_index, &
                        & selected_pair_indices, initial_vector_count, &
                        expansion_block_size, eigenvalue_tolerance, coefficient_tolerance, residual_tolerance, &
                        & iteration_limit, real_workspace, real_workspace_size, &
                        integer_workspace, integer_workspace_size, solve_high_end, iteration_count, &
                        & matrix_vector_product_count, ierr)
                    deallocate(hamiltonian_values)
                end select
                call dalloc (matrix_diagonal, "matrix_diagonal", "diagonalize_hamiltonian_block")
                call dalloc (integer_workspace, "integer_workspace", "diagonalize_hamiltonian_block")
                call dalloc (selected_pair_indices, "selected_pair_indices", "diagonalize_hamiltonian_block")

                write (summary_output_unit, *) " ", iteration_count, " iterations;"
                write (summary_output_unit, *) " ", matrix_vector_product_count, " matrix-vector multiplies."

                if (ierr /= 0) then
                    write (stderr_unit, *) "diagonalize_hamiltonian_block: Returned from solve_davidson_eigenpairs with"
                    write (stderr_unit, *) " IERR = ", ierr, "."
                    stop "diagonalize_hamiltonian_block: solve_davidson_eigenpairs wrong"
                end if

                ! Put the eigenpairs in order, overwriting as necessary

                call dcopy (nvex, real_workspace(nend+1:nend+(nvex)), 1, shifted_eigenvalues, 1)
                call dcopy (csf_count*nvex, real_workspace(1:csf_count*nvex), 1, csf_mixing_coefficients, 1)
                call dalloc (real_workspace, "real_workspace", "diagonalize_hamiltonian_block")

                ! Rearrange and reallocate storage for the eigenpairs
                ! as necessary

                if (active_level_count < highest_active_level_index) then
                    call alloc (integer_workspace, highest_active_level_index, "integer_workspace", "diagonalize_hamiltonian_block")
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
                    call dalloc (integer_workspace, "integer_workspace", "diagonalize_hamiltonian_block")
                    shifted_eigenvalues = shifted_eigenvalues(:active_level_count)
                    csf_mixing_coefficients = csf_mixing_coefficients(:int(csf_count,int64)*active_level_count)

                end if

            end if
            ! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

        end if
        ! (3e)               *** E N D   O F    A L L

        ! Only the following quantities are needed after this routine is
        ! finished:
        !   eval(), evec(), iatjpo, iaspar
        !
        !   Clean up eigenvectors; determine their J/P values
        do j = 1, active_level_count

            ! Find the dominant component of each eigenvector

            iofset = (j - 1)*csf_count

            amax = 0.0_real64
            do i = 1, csf_count
                wa = abs(csf_mixing_coefficients(i+iofset))
                if (wa <= amax) cycle
                amax = wa
                ia = i
            end do

            ! Find the angular momentum and parity of the dominant component

            iatjpo = csf_total_j_code(ia)
            iaspar = csf_parity_code(ia)

            ! Change sign of eigenvactor if dominant component is negative

            if (csf_mixing_coefficients(ia+iofset) >= 0.0_real64) cycle
            csf_mixing_coefficients(1+iofset:csf_count+iofset) = -csf_mixing_coefficients(1+iofset:csf_count+iofset)
        end do

        return
    end subroutine diagonalize_hamiltonian_block
end module diagonalize_hamiltonian_block_interface
