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
    use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_two_j
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
        matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
        nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use mcp_file_metadata, only: max_mcp_file_multipole_rank
    use mpi_state, only: mpi_process_count, mpi_rank
    use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
        csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
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
    integer, parameter :: orbital_label_radix = orbital_index_radix
    character(len=6), parameter :: routine_name = "SETCOF"
    real(real64), parameter :: eps = 1.0e-10_real64

    ! Local declarations
    integer  :: indexs(4)
    integer :: coefficient_buffer_size, current_orbital_kappa, current_orbital_two_j, packed_orbital_label_offset, &
        & packed_orbital_label_stride, partner_orbital_index, multipole_rank, symmetry_block_index, csf_index, &
        & partner_orbital_two_j, maximum_multipole_rank, &
        minimum_multipole_rank, nfile, i, jblock, jblockt, iendcdum, coefficient_count, io_status, lab, contribution_count, &
        ia, ic, ithis, irank, iind, il, iorb, iyo1, iyo2, ifound, loc1, loc2, &
        ithis2, indind, channel_record_count
    real(real64) :: current_orbital_weight, sumr, ykab, xkab, sum,  contr
    character :: mcp_channel_label*3, block_label_text*3, msg*128
    integer(int64) :: total_contribution_count, total_secondary_contribution_count
    integer :: iclmndum
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
    call alloc (coeff, coefficient_buffer_size, "COEFF", "SETCOF" )
    call alloc (iclmn, coefficient_buffer_size, "ICLMN", "SETCOF")
    call alloc (indx, coefficient_buffer_size, "INDX", "SETCOF")

    ndcof = 0
    nxcof = 0
    nycof = 0

    current_orbital_kappa = orbital_kappa(orbital_index)
    current_orbital_two_j = orbital_two_j(orbital_index)
    current_orbital_weight = ucf(orbital_index)

    !=======================================================================
    !   Generate YA coefficients that do not require MCP output list.
    !   Computation distributed and then collected.
    !=======================================================================

    packed_orbital_label_offset = 0
    packed_orbital_label_stride = orbital_label_radix*(orbital_label_radix + 1)
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
            nycof = nycof + 1
            if (nycof > nydim ) then
                if (nydim > 0) then
                    nydim = 2*nydim
                    call ralloc(nya, nydim, "NYA", "SETCOF")
                    call ralloc(ya, nydim, "YA", "SETCOF")
                else
                    nydim = 64
                    call alloc(nya, nydim, "NYA", "SETCOF")
                    call alloc(ya, nydim, "YA", "SETCOF")
                end if
            end if
            ya(nycof) = ykab
            nya(nycof) = multipole_rank + packed_orbital_label_offset
            !ENDIF

        end do
    end do

    !=======================================================================
    !   Generate XA coefficients that do not require MCP output list
    !   Computation distributed and then collected.
    !=======================================================================

    packed_orbital_label_offset = orbital_label_radix*orbital_index
    packed_orbital_label_stride = orbital_label_radix*orbital_label_radix*(orbital_label_radix + 1)
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
            nxcof = nxcof + 1
            if (nxcof > nxdim) then
                if (nxdim > 0) then
                    nxdim = 2*nxdim
                    call ralloc (xa, nxdim, "XA", "SETCOF")
                    call ralloc (nxa, nxdim, "NXA", "SETCOF")
                else
                    nxdim = 64
                    call alloc (xa, nxdim, "XA", "SETCOF")
                    call alloc (nxa, nxdim, "NXA", "SETCOF")
                end if
            end if
            xa(nxcof) = xkab
            nxa(nxcof) = multipole_rank + packed_orbital_label_offset
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
        call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "SETCOF")
        call read_mcp_structure_from_memory(jblock)

        !*** Loop over labels having non-zero coefficients
        !*** it exits when no more labels for the block

        do
            total_secondary_contribution_count = total_secondary_contribution_count + 1
            call load_mcp_records_into_memory(31,total_secondary_contribution_count,lab, contribution_count)

            ! 0, 0 marks the end of a block. This is the normal exit

            if (lab==0 .and. contribution_count==0) then
                call dalloc (matrix_partner_indices, "IROW", "SETCOF")
                exit                              ! Actually to next block
            end if

            !*** Decode the labels of I(ab) ***
            ia = mod(lab,orbital_label_radix)
            partner_orbital_index = lab/orbital_label_radix

            ! At least one orbital should be J in order to have
            ! non-zero value; otherwise, skip to the next label.

            if (ia/=orbital_index .and. partner_orbital_index/=orbital_index) then
                total_contribution_count = total_contribution_count + contribution_count  ! No contributions from this integral; skip
                cycle                             ! to next label
            end if

            if (contribution_count > coefficient_buffer_size) then
                call dalloc (coeff, "COEFF", "SETCOF")
                call dalloc (iclmn, "ICLMN", "SETCOF")
                call dalloc (indx, "INDX", "SETCOF")
                coefficient_buffer_size = contribution_count
                call alloc (coeff, coefficient_buffer_size, "COEFF", "SETCOF" )
                call alloc (iclmn, coefficient_buffer_size, "ICLMN", "SETCOF")
                call alloc (indx, coefficient_buffer_size, "INDX", "SETCOF")
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
            ! (almost) some repeatence from different blocks, a check
            ! and merge is performed. This will significantly reduce
            ! the NDCOF and thus the number of calls to YZK later.

            !*** Find the right counting parameter ***
            if (ia == orbital_index) then
                ithis = partner_orbital_index
            else
                ithis = ia
            end if

            !*** Check it against the previously recorded ***
            if (jblock > 1) then
                do i = 1, ndcof
                    if (nda(i) /= ithis) cycle     ! found, add the value
                    da(i) = da(i) + sum
                    exit
                end do
                if (i <= ndcof) cycle
            end if

            !*** Not found in the record, add an item ***
            ndcof = ndcof + 1
            if (ndcof > nddim) then
                if (nddim > 0) then
                    nddim = 2*nddim
                    call ralloc (da, nddim, "DA", "SETCOF")
                    call ralloc (nda, nddim, "NDA", "SETCOF")
                else
                    nddim = 64
                    call alloc (da, nddim, "DA", "SETCOF")
                    call alloc (nda, nddim, "NDA", "SETCOF")
                end if
            end if
            da(ndcof) = sum                      ! print*, DA(NDCOF), ndcof, myid, 'myid'
            nda(ndcof) = ithis

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
            call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "SETCOF")
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
                    call dalloc (matrix_partner_indices, "IROW", "SETCOF")
                    exit
                end if
                !***                       k
                !*** Decode the labels of R (abcd)
                indexs(4) = mod(lab,orbital_label_radix)
                lab = lab/orbital_label_radix
                indexs(2) = mod(lab,orbital_label_radix)
                lab = lab/orbital_label_radix
                indexs(3) = mod(lab,orbital_label_radix)
                indexs(1) = lab/orbital_label_radix

                !*** Determine the number of indices that match
                irank = 0
                irank = irank + count(indexs==orbital_index)

                if (irank == 0) then
                    total_contribution_count = total_contribution_count + contribution_count
                    cycle
                end if

                !*** At least one subshell index matches; allocate storage
                !*** for reading in the rest of this block
                if (contribution_count > coefficient_buffer_size) then
                    call dalloc (coeff, "COEFF", "SETCOF")
                    call dalloc (iclmn, "ICLMN", "SETCOF")
                    call dalloc (indx, "INDX", "SETCOF")
                    coefficient_buffer_size = contribution_count
                    call alloc (coeff, coefficient_buffer_size, "COEFF", "SETCOF")
                    call alloc (iclmn, coefficient_buffer_size, "ICLMN", "SETCOF")
                    call alloc (indx, coefficient_buffer_size, "INDX", "SETCOF")
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
                        if (indexs(iind) /= orbital_index) cycle   ! at least one
                        il = iind + 2
                        if (il > 4) il = il - 4
                        iorb = indexs(il)
                        il = iind + 1
                        if (il > 4) il = il - 4
                        iyo1 = indexs(il)
                        il = iind + 3
                        if (il > 4) il = il - 4
                        iyo2 = indexs(il)
                        ithis = ((iorb*orbital_label_radix + iyo2)*orbital_label_radix + iyo1)*orbital_label_radix + multipole_rank
                        exit
                    end do

                    if (ithis == (-911)) stop "ithis .EQ. -911"

                    !*** Check ithis against the previously recorded ***
                    if (jblock > 1) then
                        do i = 1, nxcof
                            if (nxa(i) /= ithis) cycle
                            xa(i) = xa(i) + sum
                            exit
                        end do
                        if (i <= nxcof) cycle
                    end if

                    !*** Not found in records, add an item ***
                    nxcof = nxcof + 1
                    if (nxcof > nxdim) then
                        if (nxdim > 0) then
                            nxdim = 2*nxdim
                            call ralloc (xa, nxdim, "XA", "SETCOF")
                            call ralloc (nxa, nxdim, "NXA", "SETCOF")
                        else
                            nxdim = 64
                            call alloc (xa, nxdim, "XA", "SETCOF")
                            call alloc (nxa, nxdim, "NXA", "SETCOF")
                        end if
                    end if
                    xa(nxcof) = sum
                    nxa(nxcof) = ithis

                case (2)

                    !=======================================================================
                    !   Two matching indices: either direct or exchange potential
                    !   contribution
                    !=======================================================================

                    ifound = 0
                    do iind = 1, 4
                        if (indexs(iind) /= orbital_index) cycle
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
                        iyo2 = indexs(il)
                        il = loc1 + 1
                        if (il > 4) il = il - 4
                        iyo1 = indexs(il)
                        ithis = (iyo2*orbital_label_radix + iyo1)*orbital_label_radix + multipole_rank

                        !*** Check it against the previously recorded ***
                        if (jblock > 1) then
                            do i = 1, nycof
                                if (nya(i) /= ithis) cycle
                                ya(i) = ya(i) + sum + sum
                                exit
                            end do
                            if (i <= nycof) cycle
                        end if

                        !*** Not found, add an item ***
                        nycof = nycof + 1
                        if (nycof > nydim) then
                            if (nydim > 0) then
                                nydim = 2*nydim
                                call ralloc(nya, nydim, "NYA", "SETCOF")
                                call ralloc(ya, nydim, "YA", "SETCOF")
                            else
                                nydim = 64
                                call alloc(nya, nydim, "NYA", "SETCOF 2b")
                                call alloc(ya, nydim, "YA", "SETCOF")
                            end if
                        end if
                        ya(nycof) = sum + sum
                        nya(nycof) = ithis

                    else
                        !
                        !   Exchange contribution
                        !
                        !*** Find ithis ***
                        il = loc1 + 2
                        if (il > 4) il = il - 4
                        iorb = indexs(il)
                        il = loc1 + 1
                        if (il > 4) il = il - 4
                        iyo1 = indexs(il)
                        il = loc1 + 3
                        if (il > 4) il = il - 4
                        iyo2 = indexs(il)
                        ithis = ((iorb*orbital_label_radix + iyo2)*orbital_label_radix + iyo1)*orbital_label_radix + multipole_rank

                        !*** Check it against the previously recorded ***
                        if (jblock > 1) then
                            do i = 1, nxcof
                                if (nxa(i) /= ithis) cycle
                                xa(i) = xa(i) + sum + sum
                                exit
                            end do
                            if (i <= nxcof) cycle
                        end if

                        !*** Not found, add an item ***
                        nxcof = nxcof + 1
                        if (nxcof > nxdim) then
                            if (nxdim > 0) then
                                nxdim = 2*nxdim
                                call ralloc (xa, nxdim, "XA", "SETCOF")
                                call ralloc (nxa, nxdim, "NXA", "SETCOF")
                            else
                                nxdim = 64
                                call alloc (xa, nxdim, "XA", "SETCOF")
                                call alloc (nxa, nxdim, "NXA", "SETCOF")
                            end if
                        end if
                        xa(nxcof) = sum + sum
                        nxa(nxcof) = ithis

                    end if

                case (3)

                    !=======================================================================
                    !   Three matching indices: direct and exchange potential contributions
                    !=======================================================================

                    !*** Find ithis AND ithis2
                    ithis = -911
                    ithis2 = -911
                    do iind = 1, 4
                        if (indexs(iind) == orbital_index) cycle
                        indind = indexs(iind)
                        iyo2 = indind
                        iyo1 = orbital_index
                        ithis = (iyo2*orbital_label_radix + iyo1)*orbital_label_radix + multipole_rank
                        iorb = indind
                        iyo1 = orbital_index
                        iyo2 = orbital_index
                        ithis2 = ((iorb*orbital_label_radix + iyo2)*orbital_label_radix + iyo1)*orbital_label_radix + multipole_rank
                    end do

                    if (ithis==(-911) .or. ithis2==(-911)) stop "ithis2"

                    !*** Check the previously recorded for YA
                    if (jblock > 1) then
                        do i = 1, nycof
                            if (nya(i) /= ithis) cycle
                            ya(i) = ya(i) + sum + sum
                            exit
                        end do
                    end if

                    ! Add a direct contribution if the search found no match.
                    if (jblock <= 1) i = nycof + 1
                    if (i > nycof) then
                        nycof = nycof + 1
                        if (nycof > nydim) then
                            if (nydim > 0) then
                                nydim = 2*nydim
                                call ralloc (nya, nydim, "NYA", "SETCOF")
                                call ralloc (ya, nydim, "YA", "SETCOF")
                            else
                                nydim = 64
                                call alloc (nya, nydim, "NYA", "SETCOF")
                                call alloc (ya, nydim, "YA", "SETCOF")
                            end if
                        end if
                        ya(nycof) = sum + sum
                        nya(nycof) = ithis

                    end if

                    !*** Check the previously recorded for XA
                    if (jblock > 1) then
                        do i = 1, nxcof
                            if (nxa(i) /= ithis2) cycle
                            xa(i) = xa(i) + sum
                            exit
                        end do
                        if (i <= nxcof) cycle
                    end if

                    ! Not found, add an item ***
                    nxcof = nxcof + 1
                    if (nxcof > nxdim) then
                        if (nxdim > 0) then
                            nxdim = 2*nxdim
                            call ralloc (xa, nxdim, "XA", "SETCOF")
                            call ralloc (nxa, nxdim, "NXA", "SETCOF")
                        else
                            nxdim = 64
                            call alloc (xa, nxdim, "XA", "SETCOF")
                            call alloc (nxa, nxdim, "NXA", "SETCOF")
                        end if
                    end if
                    xa(nxcof) = sum
                    nxa(nxcof) = ithis2

                case (4)

                    !=======================================================================
                    !   Four matching indices: direct potential contribution
                    !=======================================================================

                    !*** Find ithis AND ithis2
                    iyo2 = orbital_index
                    iyo1 = orbital_index
                    ithis = (iyo2*orbital_label_radix + iyo1)*orbital_label_radix + multipole_rank

                    !*** Check the previously recorded for YA
                    if (jblock > 1) then
                        do i = 1, nycof
                            if (nya(i) /= ithis) cycle
                            ya(i) = ya(i) + 4.0_real64*sum
                            exit
                        end do
                        if (i <= nycof) cycle
                    end if

                    ! Not found, add an item ***
                    nycof = nycof + 1
                    if (nycof > nydim) then
                        if (nydim > 0) then
                            nydim = 2*nydim
                            call ralloc(nya, nydim, "NYA", "SETCOF")
                            call ralloc(ya, nydim, "YA", "SETCOF")
                        else
                            nydim = 64
                            call alloc(nya, nydim, "NYA", "SETCOF 4b")
                            call alloc(ya, nydim, "YA", "SETCOF")
                        end if
                    end if
                    ya(nycof) = 4.0_real64*sum
                    nya(nycof) = ithis

                end select

            end do
        end do                                  ! loop for V-Coefficients
    end do

    !=======================================================================
    !   Deallocate storage for arrays local to this routine
    !=======================================================================
    call dalloc (coeff, "COEFF", "SETCOF")
    call dalloc (iclmn, "ICLMN", "SETCOF")
    call dalloc (indx, "INDX", "SETCOF")

    return
end subroutine assemble_orbital_potential_coefficients
