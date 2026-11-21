module assemble_orbital_potential_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_orbital_potential_coefficients
contains

    !***********************************************************************
    !                                                                      *
    subroutine assemble_orbital_potential_coefficients(extended_optimal_level, orbital_index)
        ! Module dependencies
        !                                                                      *
        !   This  subroutine  sets  up the coefficients and orbital pointers   *
        !   for the direct and exchange potentials for orbital  J .  It also   *
        !   sets  up  the  coefficients  and  pointers for the inhomogeneous   *
        !   terms arising from off-diagonal I (a,b) integrals.                 *
        !                                                                      *
        !   Call(s) to: [LIB92]: alloc, dalloc.                                *
        !               [RSCF92]: alcsca, dsubrs.                              *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 21 Dec 1992   *
        !   Modified by Xinghong He                 Last update: 21 Dec 1997   *
        !   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
        !      Spin-angular coefficiants have been put on memory               *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use grasp_parameters,    only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc, ralloc
        use orbital_csf_basis, only: orbital_count, orbital_kappa, orbital_two_j
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, matrix_partner_indices
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use serial_mpi_context, only: mpi_process_count, mpi_rank
        use symmetry_block_offsets, only: csf_offsets_by_block
        use scf_iteration_state, only: orbital_occupation_weights, direct_coefficient_count, exchange_coefficient_count, &
            y_potential_coefficient_count, direct_coefficient_capacity, exchange_coefficient_capacity, &
            y_potential_coefficient_capacity, direct_coefficients, exchange_coefficients, y_potential_coefficients, &
            direct_orbital_indices, exchange_coefficient_labels, y_potential_coefficient_labels
        use sparse_hamiltonian_storage,           only: stored_matrix_element_count
        use rmcdhf_mcp_memory,     only: indx, iclmn, coeff, imin2_30, imax2_30

        ! Procedure interfaces
        !      USE cons_C,           ONLY: EPS
        use csf_density_matrix_element_interface, only: csf_density_matrix_element
        use direct_angular_coefficient_interface, only: direct_angular_coefficient
        use exchange_angular_coefficient_interface, only: exchange_angular_coefficient
        use load_mcp_records_into_memory_interface, only: load_mcp_records_into_memory
        use read_mcp_coefficients_into_memory_interface, only: read_mcp_coefficients_into_memory
        use read_mcp_structure_from_memory_interface, only: read_mcp_structure_from_memory
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        logical  :: extended_optimal_level

        ! Local constants
        character(len=*), parameter :: routine_name = "assemble_orbital_potential_coefficients"
        real(real64), parameter :: eps = 1.0e-10_real64

        ! Local declarations
        integer  :: decoded_integral_orbital_indices(4)
        integer :: coefficient_buffer_size, current_orbital_kappa, current_orbital_two_j, packed_orbital_label_offset, &
            packed_orbital_label_stride, partner_orbital_index, multipole_rank, symmetry_block_index, csf_index, &
            partner_orbital_two_j, maximum_multipole_rank, minimum_multipole_rank, nfile, i, jblock, lab, contribution_count, &
            ia, ic, ithis, irank, iind, il, iorb, iyo1, iyo2, ifound, loc1, loc2, ithis2, selected_integral_orbital_index, &
                & channel_record_count
        real(real64) :: current_orbital_weight, sumr, ykab, xkab, sum,  contr
        integer(int64) :: total_contribution_count, total_secondary_contribution_count
        !-----------------------------------------------
        !** Locals
        !     POINTER (PCOEFF,COEFF(1))
        !     POINTER (PICLMN,ICLMN(1))
        !     POINTER (PINDEX,indx(1))
        !
        !     POINTER (PNIROW,IROW(1))
        !     POINTER (PNTRDA,DA(1))
        !     POINTER (PNTRXA,XA(1))
        !     POINTER (PNTRYA,YA(1))
        !     POINTER (PNTNDA,NDA(1))
        !     POINTER (PNTNXA,NXA(1))
        !     POINTER (PNTNYA,NYA(1))
        !
        !
        !     POINTER (pncfblk, ncfblk(0:*))
        !
        !     POINTER (pnevblk, nevblk(1))
        !     POINTER (pncmaxblk, ncmaxblk(1))
        !
        !
        !     POINTER (pncfpast, ncfpast(1))
        !     POINTER (pncminpast, ncminpast(1))
        !     POINTER (pnevecpast, nevecpast(1))
        !
        !     INTEGER, DIMENSION(:), pointer :: nda, nxa, nya
        !     INTEGER, DIMENSION(:), pointer :: ncfblk, nevblk, ncmaxblk, &
        !                                       ncfpast, ncminpast, nevecpast

        !=======================================================================
        !   Initializations
        !=======================================================================

        coefficient_buffer_size = 1
        call alloc (coeff, coefficient_buffer_size, "COEFF", "assemble_orbital_potential_coefficients" )
        call alloc (iclmn, coefficient_buffer_size, "ICLMN", "assemble_orbital_potential_coefficients")
        call alloc (indx, coefficient_buffer_size, "INDX", "assemble_orbital_potential_coefficients")

        direct_coefficient_count = 0
        exchange_coefficient_count = 0
        y_potential_coefficient_count = 0

        current_orbital_kappa = orbital_kappa(orbital_index)
        current_orbital_two_j = orbital_two_j(orbital_index)
        current_orbital_weight = orbital_occupation_weights(orbital_index)

        !=======================================================================
        !   Generate YA coefficients that do not require MCP output list.
        !   Computation distributed and then collected.
        !=======================================================================

        packed_orbital_label_offset = 0
        packed_orbital_label_stride = orbital_index_radix*(orbital_index_radix + 1)
        do partner_orbital_index = 1, orbital_count
            packed_orbital_label_offset = packed_orbital_label_offset + packed_orbital_label_stride
            if (partner_orbital_index == orbital_index) then
                maximum_multipole_rank = current_orbital_two_j - 1
            else
                maximum_multipole_rank = 0
            end if
            do multipole_rank = 0, maximum_multipole_rank, 2

                !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
                sumr = 0.0_real64
                do symmetry_block_index = 1, symmetry_block_count
                    do csf_index = mpi_rank + 1, csf_counts_by_block(symmetry_block_index), mpi_process_count
                        sumr = sumr + csf_density_matrix_element(extended_optimal_level,csf_index,csf_index, &
                            & symmetry_block_index)*direct_angular_coefficient(multipole_rank, &
                            & csf_index + csf_offsets_by_block(symmetry_block_index),orbital_index,&
                            partner_orbital_index)
                    end do
                end do
                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                if (partner_orbital_index == orbital_index) then
                    ykab = 2.0_real64*sumr/current_orbital_weight
                else
                    ykab = sumr/current_orbital_weight
                end if

                !*** The following IF has to be removed since YKAB
                !*** is incomplete in multiprocessor case
                !IF (ABS (YKAB) .GT. EPS) THEN
                y_potential_coefficient_count = y_potential_coefficient_count + 1
                if (y_potential_coefficient_count > y_potential_coefficient_capacity ) then
                    if (y_potential_coefficient_capacity > 0) then
                        y_potential_coefficient_capacity = 2*y_potential_coefficient_capacity
                        !cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                        call ralloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                            "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                        call ralloc(y_potential_coefficients, y_potential_coefficient_capacity, "y_potential_coefficients", &
                            "assemble_orbital_potential_coefficients")
                    !                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
                    !                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                    else
                        !cjb similar here
                        y_potential_coefficient_capacity = 64
                        call alloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                            "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                        call alloc(y_potential_coefficients, y_potential_coefficient_capacity, "y_potential_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        !                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
                        !                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
                        !                 NYDIM = NYCOF
                        !cjb
                    end if
                end if
                y_potential_coefficients(y_potential_coefficient_count) = ykab
                y_potential_coefficient_labels(y_potential_coefficient_count) = multipole_rank + packed_orbital_label_offset
                !ENDIF

            end do
        end do

        !=======================================================================
        !   Generate XA coefficients that do not require MCP output list
        !   Computation distributed and then collected.
        !=======================================================================

        packed_orbital_label_offset = orbital_index_radix*orbital_index
        packed_orbital_label_stride = orbital_index_radix*orbital_index_radix*(orbital_index_radix + 1)
        do partner_orbital_index = 1, orbital_count
            packed_orbital_label_offset = packed_orbital_label_offset + packed_orbital_label_stride
            if (partner_orbital_index == orbital_index) cycle
            partner_orbital_two_j = orbital_two_j(partner_orbital_index)
            if (current_orbital_kappa*orbital_kappa(partner_orbital_index) > 0) then
                minimum_multipole_rank = abs((current_orbital_two_j - partner_orbital_two_j)/2)
            else
                minimum_multipole_rank = abs((current_orbital_two_j - partner_orbital_two_j)/2) + 1
            end if
            maximum_multipole_rank = (current_orbital_two_j + partner_orbital_two_j)/2

            do multipole_rank = minimum_multipole_rank, maximum_multipole_rank, 2

                !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
                sumr = 0.0_real64
                do symmetry_block_index = 1, symmetry_block_count
                    do csf_index = mpi_rank + 1, csf_counts_by_block(symmetry_block_index), mpi_process_count
                        sumr = sumr + csf_density_matrix_element(extended_optimal_level,csf_index,csf_index, &
                            & symmetry_block_index)*exchange_angular_coefficient(multipole_rank, &
                            & csf_index + csf_offsets_by_block(symmetry_block_index),orbital_index,&
                            partner_orbital_index)
                    end do
                end do
                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                xkab = sumr/current_orbital_weight

                !*** The following IF has to be removed since YKAB
                !*** is incomplete in multiprocessor case
                !IF (ABS (XKAB) .GT. EPS) THEN
                exchange_coefficient_count = exchange_coefficient_count + 1
                if (exchange_coefficient_count > exchange_coefficient_capacity) then
                    if (exchange_coefficient_capacity > 0) then
                        exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                        !cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                        call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, "exchange_coefficient_labels", &
                            "assemble_orbital_potential_coefficients")
                    !                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
                    !                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                    else
                        exchange_coefficient_capacity = 64
                        call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, "exchange_coefficient_labels", &
                            "assemble_orbital_potential_coefficients")
                        !                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
                        !                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                        !                  NXDIM = NXCOF
                        !cjb
                    end if
                end if
                exchange_coefficients(exchange_coefficient_count) = xkab
                exchange_coefficient_labels(exchange_coefficient_count) = multipole_rank + packed_orbital_label_offset
                !ENDIF
            end do
        end do

        !=======================================================================
        !   Generate DA coefficients; these arise from the one-electron
        !   integrals
        !=======================================================================

        total_contribution_count = 0
        total_secondary_contribution_count = 0
        do jblock = 1, symmetry_block_count

            !*** Read in IROW from file mcp.30 ***
            channel_record_count = imax2_30(jblock) - imin2_30(jblock)
            stored_matrix_element_count = int(channel_record_count,kind=int64)
            call alloc (matrix_partner_indices, stored_matrix_element_count, "matrix_partner_indices", &
                "assemble_orbital_potential_coefficients")
            call read_mcp_structure_from_memory(jblock)

            !*** Loop over labels having non-zero coefficients
            !*** it exits when no more labels for the block

            do
                total_secondary_contribution_count = total_secondary_contribution_count + 1
                call load_mcp_records_into_memory(31,total_secondary_contribution_count,lab, contribution_count)

                ! 0, 0 marks the end of a block. This is the normal exit

                if (lab==0 .and. contribution_count==0) then
                    call dalloc (matrix_partner_indices, "matrix_partner_indices", "assemble_orbital_potential_coefficients")
                    exit                              ! Actually to next block
                end if

                !*** Decode the labels of I(ab) ***
                ia = mod(lab,orbital_index_radix)
                partner_orbital_index = lab/orbital_index_radix

                ! At least one orbital should be J in order to have
                ! non-zero value; otherwise, skip to the next label.

                if (ia/=orbital_index .and. partner_orbital_index/=orbital_index) then
                    total_contribution_count = total_contribution_count + contribution_count  ! No contributions from this integral; skip
                    cycle                             ! to next label
                end if

                if (contribution_count > coefficient_buffer_size) then
                    call dalloc (coeff, "COEFF", "assemble_orbital_potential_coefficients")
                    call dalloc (iclmn, "ICLMN", "assemble_orbital_potential_coefficients")
                    call dalloc (indx, "INDX", "assemble_orbital_potential_coefficients")
                    coefficient_buffer_size = contribution_count
                    call alloc (coeff, coefficient_buffer_size, "COEFF", "assemble_orbital_potential_coefficients" )
                    call alloc (iclmn, coefficient_buffer_size, "ICLMN", "assemble_orbital_potential_coefficients")
                    call alloc (indx, coefficient_buffer_size, "INDX", "assemble_orbital_potential_coefficients")
                end if

                ! Read the column index, the sparse matrix index, and the
                ! coefficient for all contributions from this integral

                call read_mcp_coefficients_into_memory(31,total_contribution_count,contribution_count)
                total_contribution_count = total_contribution_count + contribution_count

                ! Add up all the contributions from this integral;
                ! off-diagonal contributions have double the weight

                sum = 0.0_real64
                do i = 1, contribution_count
                    csf_index = matrix_partner_indices(indx(i))
                    ic = iclmn(i)

                    contr = csf_density_matrix_element(extended_optimal_level,csf_index,ic,jblock)*coeff(i)

                    if (csf_index /= ic) contr = contr + contr
                    sum = sum + contr
                end do
                sum = 0.5_real64*sum/current_orbital_weight

                ! Put coefficients in the list. Since there is always
                ! (almost) some repetition from different blocks, a check
                ! and merge is performed. This will significantly reduce
                ! the number of direct coefficients and calls to the Hartree Y and Z evaluator.

                !*** Find the right counting parameter ***
                if (ia == orbital_index) then
                    ithis = partner_orbital_index
                else
                    ithis = ia
                end if

                !*** Check it against the previously recorded ***
                if (jblock > 1) then
                    do i = 1, direct_coefficient_count
                        if (direct_orbital_indices(i) /= ithis) cycle     ! found, add the value
                        direct_coefficients(i) = direct_coefficients(i) + sum
                        exit
                    end do
                    if (i <= direct_coefficient_count) cycle
                end if

                !*** Not found in the record, add an item ***
                direct_coefficient_count = direct_coefficient_count + 1
                if (direct_coefficient_count > direct_coefficient_capacity) then
                    if (direct_coefficient_capacity > 0) then
                        direct_coefficient_capacity = 2*direct_coefficient_capacity
                        call ralloc (direct_coefficients, direct_coefficient_capacity, "direct_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call ralloc (direct_orbital_indices, direct_coefficient_capacity, "direct_orbital_indices", &
                            "assemble_orbital_potential_coefficients")
                    !                  CALL RALLOC (DA, NDCOF, 'DA', 'SETCOF')
                    !                  CALL RALLOC (NDA, NDCOF, 'NDA', 'SETCOF')
                    else
                        direct_coefficient_capacity = 64
                        call alloc (direct_coefficients, direct_coefficient_capacity, "direct_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call alloc (direct_orbital_indices, direct_coefficient_capacity, "direct_orbital_indices", &
                            "assemble_orbital_potential_coefficients")
                        !                  CALL ALLOC (DA, NDCOF, 'DA', 'SETCOF')
                        !                  CALL ALLOC (NDA, NDCOF, 'NDA', 'SETCOF')
                        !                  NDDIM = NDCOF
                    end if
                end if
                direct_coefficients(direct_coefficient_count) = sum                      ! print*, DA(NDCOF), ndcof, myid, 'myid'
                direct_orbital_indices(direct_coefficient_count) = ithis

            end do                             ! For labels
        end do                                     ! For blocks

        !=======================================================================
        !   Generate YA and XA coefficients; these arise from the two-electron
        !   integrals
        !=======================================================================

        do nfile = 32, 32 + max_mcp_file_multipole_rank

            !=======================================================================
            !   Loop over blocks again, this time, for V-coefficients
            !=======================================================================

            total_contribution_count = 0
            total_secondary_contribution_count = 0
            do jblock = 1, symmetry_block_count
                !           ...Read in IROW from file mcp.30
                channel_record_count = imax2_30(jblock) - imin2_30(jblock)
                stored_matrix_element_count = int(channel_record_count,kind=int64)
                call alloc (matrix_partner_indices, stored_matrix_element_count, "matrix_partner_indices", &
                    "assemble_orbital_potential_coefficients")
                call read_mcp_structure_from_memory(jblock)

                multipole_rank = nfile - 32                       ! multipolarity of the integral

                !=======================================================================
                !   Attempt to read another block of data
                !=======================================================================
                do
                    total_secondary_contribution_count = total_secondary_contribution_count + 1
                    call load_mcp_records_into_memory(nfile,total_secondary_contribution_count,lab, contribution_count)
                    !
                    if (lab==0 .and. contribution_count==0) then
                        call dalloc (matrix_partner_indices, "matrix_partner_indices", "assemble_orbital_potential_coefficients")
                        exit
                    end if
                    !***                       k
                    !*** Decode the labels of R (abcd)
                    decoded_integral_orbital_indices(4) = mod(lab,orbital_index_radix)
                    lab = lab/orbital_index_radix
                    decoded_integral_orbital_indices(2) = mod(lab,orbital_index_radix)
                    lab = lab/orbital_index_radix
                    decoded_integral_orbital_indices(3) = mod(lab,orbital_index_radix)
                    decoded_integral_orbital_indices(1) = lab/orbital_index_radix

                    !*** Determine the number of indices that match
                    irank = 0
                    irank = irank + count(decoded_integral_orbital_indices==orbital_index)

                    if (irank == 0) then
                        total_contribution_count = total_contribution_count + contribution_count
                        cycle
                    end if

                    !*** At least one subshell index matches; allocate storage
                    !*** for reading in the rest of this block
                    if (contribution_count > coefficient_buffer_size) then
                        call dalloc (coeff, "COEFF", "assemble_orbital_potential_coefficients")
                        call dalloc (iclmn, "ICLMN", "assemble_orbital_potential_coefficients")
                        call dalloc (indx, "INDX", "assemble_orbital_potential_coefficients")
                        coefficient_buffer_size = contribution_count
                        call alloc (coeff, coefficient_buffer_size, "COEFF", "assemble_orbital_potential_coefficients")
                        call alloc (iclmn, coefficient_buffer_size, "ICLMN", "assemble_orbital_potential_coefficients")
                        call alloc (indx, coefficient_buffer_size, "INDX", "assemble_orbital_potential_coefficients")
                    end if

                    !*** Read column index, sparse matrix index, and
                    !*** coefficient for all contributions from this integral
                    call read_mcp_coefficients_into_memory(nfile,total_contribution_count,contribution_count)
                    total_contribution_count = total_contribution_count + contribution_count

                    !*** Add up all the contributions from this integral;
                    !*** off-diagonal contributions have double the weight
                    sum = 0.0_real64
                    do i = 1, contribution_count
                        csf_index = matrix_partner_indices(indx(i))
                        ic = iclmn(i)
                        contr = csf_density_matrix_element(extended_optimal_level,csf_index,ic,jblock)*coeff(i)
                        if (csf_index /= ic) contr = contr + contr
                        sum = sum + contr
                    end do
                    sum = 0.5_real64*sum/current_orbital_weight

                    select case (irank)
                    case (1)

                        !=======================================================================
                        !   One matching index: exchange potential contribution
                        !=======================================================================

                        !*** Similar to DA, find ithis ***
                        ithis = -911                      ! initialize to an impossible value
                        ! though not necessary
                        do iind = 1, 4
                            if (decoded_integral_orbital_indices(iind) /= orbital_index) cycle   ! at least one
                            il = iind + 2
                            if (il > 4) il = il - 4
                            iorb = decoded_integral_orbital_indices(il)
                            il = iind + 1
                            if (il > 4) il = il - 4
                            iyo1 = decoded_integral_orbital_indices(il)
                            il = iind + 3
                            if (il > 4) il = il - 4
                            iyo2 = decoded_integral_orbital_indices(il)
                            ithis = ((iorb*orbital_index_radix + iyo2)*orbital_index_radix + iyo1)*orbital_index_radix + &
                                & multipole_rank
                            exit
                        end do

                        if (ithis == (-911)) stop "ithis .EQ. -911"

                        !*** Check ithis against the previously recorded ***
                        if (jblock > 1) then
                            do i = 1, exchange_coefficient_count
                                if (exchange_coefficient_labels(i) /= ithis) cycle
                                exchange_coefficients(i) = exchange_coefficients(i) + sum
                                exit
                            end do
                            if (i <= exchange_coefficient_count) cycle
                        end if

                        !*** Not found in records, add an item ***
                        exchange_coefficient_count = exchange_coefficient_count + 1
                        if (exchange_coefficient_count > exchange_coefficient_capacity) then
                            if (exchange_coefficient_capacity > 0) then
                                exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                                !cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                                call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                            !                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
                            !                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                            else
                                exchange_coefficient_capacity = 64
                                call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                                !                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
                                !                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                                !                  NXDIM = NXCOF
                                !cjb
                            end if
                        end if
                        exchange_coefficients(exchange_coefficient_count) = sum
                        exchange_coefficient_labels(exchange_coefficient_count) = ithis

                    case (2)

                        !=======================================================================
                        !   Two matching indices: either direct or exchange potential
                        !   contribution
                        !=======================================================================

                        ifound = 0
                        do iind = 1, 4
                            if (decoded_integral_orbital_indices(iind) /= orbital_index) cycle
                            if (ifound == 0) then
                                loc1 = iind
                                ifound = ifound + 1
                            else if (ifound == 1) then
                                loc2 = iind
                                exit
                            end if
                        end do

                        if (loc2 - loc1 == 2) then
                            !
                            !   Direct contribution
                            !
                            !*** Find ithis ***
                            il = loc1 + 3
                            if (il > 4) il = il - 4
                            iyo2 = decoded_integral_orbital_indices(il)
                            il = loc1 + 1
                            if (il > 4) il = il - 4
                            iyo1 = decoded_integral_orbital_indices(il)
                            ithis = (iyo2*orbital_index_radix + iyo1)*orbital_index_radix + multipole_rank

                            !*** Check it against the previously recorded ***
                            if (jblock > 1) then
                                do i = 1, y_potential_coefficient_count
                                    if (y_potential_coefficient_labels(i) /= ithis) cycle
                                    y_potential_coefficients(i) = y_potential_coefficients(i) + sum + sum
                                    exit
                                end do
                                if (i <= y_potential_coefficient_count) cycle
                            end if

                            !*** Not found, add an item ***
                            y_potential_coefficient_count = y_potential_coefficient_count + 1
                            if (y_potential_coefficient_count > y_potential_coefficient_capacity) then
                                if (y_potential_coefficient_capacity > 0) then
                                    y_potential_coefficient_capacity = 2*y_potential_coefficient_capacity
                                    !cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                                    call ralloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    call ralloc(y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                !                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
                                !                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                                else
                                    y_potential_coefficient_capacity = 64
                                    call alloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    call alloc(y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                    !                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
                                    !                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
                                    !                 NYDIM = NYCOF
                                    !cjb
                                end if
                            end if
                            y_potential_coefficients(y_potential_coefficient_count) = sum + sum
                            y_potential_coefficient_labels(y_potential_coefficient_count) = ithis

                        else
                            !
                            !   Exchange contribution
                            !
                            !*** Find ithis ***
                            il = loc1 + 2
                            if (il > 4) il = il - 4
                            iorb = decoded_integral_orbital_indices(il)
                            il = loc1 + 1
                            if (il > 4) il = il - 4
                            iyo1 = decoded_integral_orbital_indices(il)
                            il = loc1 + 3
                            if (il > 4) il = il - 4
                            iyo2 = decoded_integral_orbital_indices(il)
                            ithis = ((iorb*orbital_index_radix + iyo2)*orbital_index_radix + iyo1)*orbital_index_radix + &
                                & multipole_rank

                            !*** Check it against the previously recorded ***
                            if (jblock > 1) then
                                do i = 1, exchange_coefficient_count
                                    if (exchange_coefficient_labels(i) /= ithis) cycle
                                    exchange_coefficients(i) = exchange_coefficients(i) + sum + sum
                                    exit
                                end do
                                if (i <= exchange_coefficient_count) cycle
                            end if

                            !*** Not found, add an item ***
                            exchange_coefficient_count = exchange_coefficient_count + 1
                            if (exchange_coefficient_count > exchange_coefficient_capacity) then
                                if (exchange_coefficient_capacity > 0) then
                                    exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                                    !cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                                    call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                        "assemble_orbital_potential_coefficients")
                                    call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                        "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                                !                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
                                !                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                                else
                                    exchange_coefficient_capacity = 64
                                    call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                        "assemble_orbital_potential_coefficients")
                                    call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                        "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    !                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
                                    !                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                                    !                  NXDIM = NXCOF
                                    !cjb
                                end if
                            end if
                            exchange_coefficients(exchange_coefficient_count) = sum + sum
                            exchange_coefficient_labels(exchange_coefficient_count) = ithis

                        end if

                    case (3)

                        !=======================================================================
                        !   Three matching indices: direct and exchange potential contributions
                        !=======================================================================

                        !*** Find ithis AND ithis2
                        ithis = -911
                        ithis2 = -911
                        do iind = 1, 4
                            if (decoded_integral_orbital_indices(iind) == orbital_index) cycle
                            selected_integral_orbital_index = decoded_integral_orbital_indices(iind)
                            iyo2 = selected_integral_orbital_index
                            iyo1 = orbital_index
                            ithis = (iyo2*orbital_index_radix + iyo1)*orbital_index_radix + multipole_rank
                            iorb = selected_integral_orbital_index
                            iyo1 = orbital_index
                            iyo2 = orbital_index
                            ithis2 = ((iorb*orbital_index_radix + iyo2)*orbital_index_radix + iyo1)*orbital_index_radix + &
                                & multipole_rank
                        end do

                        if (ithis==(-911) .or. ithis2==(-911)) stop "ithis2"

                        !*** Check the previously recorded for YA
                        if (jblock > 1) then
                            do i = 1, y_potential_coefficient_count
                                if (y_potential_coefficient_labels(i) /= ithis) cycle
                                y_potential_coefficients(i) = y_potential_coefficients(i) + sum + sum
                                exit
                            end do
                        end if

                        ! Add a direct contribution if the search found no match.
                        if (jblock <= 1) i = y_potential_coefficient_count + 1
                        if (i > y_potential_coefficient_count) then
                            y_potential_coefficient_count = y_potential_coefficient_count + 1
                            if (y_potential_coefficient_count > y_potential_coefficient_capacity) then
                                if (y_potential_coefficient_capacity > 0) then
                                    y_potential_coefficient_capacity = 2*y_potential_coefficient_capacity
                                    call ralloc (y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    call ralloc (y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                !                     CALL RALLOC (NYA, NYCOF, 'NYA', 'SETCOF')
                                !                     CALL RALLOC (YA, NYCOF, 'YA', 'SETCOF')
                                else
                                    y_potential_coefficient_capacity = 64
                                    call alloc (y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    call alloc (y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                    !                     CALL ALLOC (NYA, NYCOF, 'NYA', 'SETCOF')
                                    !                     CALL ALLOC (YA, NYCOF, 'YA', 'SETCOF')
                                    !                     NYDIM = NYCOF
                                end if
                            end if
                            y_potential_coefficients(y_potential_coefficient_count) = sum + sum
                            y_potential_coefficient_labels(y_potential_coefficient_count) = ithis

                        end if

                        !*** Check the previously recorded for XA
                        if (jblock > 1) then
                            do i = 1, exchange_coefficient_count
                                if (exchange_coefficient_labels(i) /= ithis2) cycle
                                exchange_coefficients(i) = exchange_coefficients(i) + sum
                                exit
                            end do
                            if (i <= exchange_coefficient_count) cycle
                        end if

                        ! Not found, add an item ***
                        exchange_coefficient_count = exchange_coefficient_count + 1
                        if (exchange_coefficient_count > exchange_coefficient_capacity) then
                            if (exchange_coefficient_capacity > 0) then
                                exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                                !cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                                call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                            !                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
                            !                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                            else
                                call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                                !                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
                                !                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                                !                  NXDIM = NXCOF
                                !cjb
                            end if
                        end if
                        exchange_coefficients(exchange_coefficient_count) = sum
                        exchange_coefficient_labels(exchange_coefficient_count) = ithis2

                    case (4)

                        !=======================================================================
                        !   Four matching indices: direct potential contribution
                        !=======================================================================

                        !*** Find ithis AND ithis2
                        iyo2 = orbital_index
                        iyo1 = orbital_index
                        ithis = (iyo2*orbital_index_radix + iyo1)*orbital_index_radix + multipole_rank

                        !*** Check the previously recorded for YA
                        if (jblock > 1) then
                            do i = 1, y_potential_coefficient_count
                                if (y_potential_coefficient_labels(i) /= ithis) cycle
                                y_potential_coefficients(i) = y_potential_coefficients(i) + 4.0_real64*sum
                                exit
                            end do
                            if (i <= y_potential_coefficient_count) cycle
                        end if

                        ! Not found, add an item ***
                        y_potential_coefficient_count = y_potential_coefficient_count + 1
                        if (y_potential_coefficient_count > y_potential_coefficient_capacity) then
                            if (y_potential_coefficient_capacity > 0) then
                                y_potential_coefficient_capacity = 2*y_potential_coefficient_capacity
                                !cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                                call ralloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                    "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                call ralloc(y_potential_coefficients, y_potential_coefficient_capacity, &
                                    "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                            !                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
                            !                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                            else
                                !cjb similar here
                                y_potential_coefficient_capacity = 64
                                call alloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                    "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                call alloc(y_potential_coefficients, y_potential_coefficient_capacity, "y_potential_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                !                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
                                !                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
                                !                 NYDIM = NYCOF
                                !cjb
                            end if
                        end if
                        y_potential_coefficients(y_potential_coefficient_count) = 4.0_real64*sum
                        y_potential_coefficient_labels(y_potential_coefficient_count) = ithis

                    end select

                end do
            end do                                  ! loop for V-Coefficients
        end do

        !=======================================================================
        !   Deallocate storage for arrays local to this routine
        !=======================================================================
        call dalloc (coeff, "COEFF", "assemble_orbital_potential_coefficients")
        call dalloc (iclmn, "ICLMN", "assemble_orbital_potential_coefficients")
        call dalloc (indx, "INDX", "assemble_orbital_potential_coefficients")

        return
    end subroutine assemble_orbital_potential_coefficients
end module assemble_orbital_potential_coefficients_interface
