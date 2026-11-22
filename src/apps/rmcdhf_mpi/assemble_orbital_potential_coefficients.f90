! Contract MCP contributions with the current mixing coefficients to form orbital-equation
! coefficients.
module assemble_orbital_potential_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_orbital_potential_coefficients
contains

    ! This  subroutine  sets  up the coefficients and orbital pointers
    ! for the direct and exchange potentials for orbital  J .  It also
    ! sets  up  the  coefficients  and  pointers for the inhomogeneous
    ! terms arising from off-diagonal I (a,b) integrals.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 21 Dec 1992
    ! Modified by Xinghong He. Last update: 21 Dec 1997
    subroutine assemble_orbital_potential_coefficients(extended_optimal_level, orbital_index)
        ! Module dependencies
        use mpi_indexed_mcp, only: decode_mcp_label
        use mcp_density_contractions, only: prepare_contractions, contraction_bounds, contracted_labels, contracted_values
        use rmcdhf_mcp_io, only: mcp_now, profile_stage, log_bits, verification_enabled
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use grasp_parameters,    only: orbital_index_radix
        use array_memory_management, only: alloc, ralloc
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_two_j
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use mpi_state, only: mpi_process_count, mpi_rank
        use symmetry_block_offsets, only: csf_offsets_by_block
        use scf_iteration_state, only: orbital_occupation_weights, direct_coefficient_count, exchange_coefficient_count, &
            y_potential_coefficient_count, direct_coefficient_capacity, exchange_coefficient_capacity, &
            y_potential_coefficient_capacity, direct_coefficients, exchange_coefficients, y_potential_coefficients, &
            direct_orbital_indices, exchange_coefficient_labels, y_potential_coefficient_labels

        ! Procedure interfaces
        use csf_density_matrix_element_interface, only: csf_density_matrix_element
        use direct_angular_coefficient_interface, only: direct_angular_coefficient
        use exchange_angular_coefficient_interface, only: exchange_angular_coefficient
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        logical, intent(in) :: extended_optimal_level

        real(real64), parameter :: eps = 1.0e-10_real64
        integer :: decoded_integral_orbital_indices(4)
        integer :: current_orbital_kappa, current_orbital_two_j, packed_orbital_label_offset, packed_orbital_label_stride, &
            partner_orbital_index, multipole_rank, symmetry_block_index, csf_index, partner_orbital_two_j, &
            maximum_multipole_rank, minimum_multipole_rank, nfile, i, jblock, lab, ia, ithis, irank, iind, il, iorb, &
            iyo1, iyo2, ifound, loc1, loc2, ithis2, selected_integral_orbital_index
        integer(int64) :: first, last, group_index
        real(real64) :: current_orbital_weight, sumr, ykab, xkab, sum, started

        started = mcp_now()
        call prepare_contractions(extended_optimal_level)
        direct_coefficient_count = 0
        exchange_coefficient_count = 0
        y_potential_coefficient_count = 0

        current_orbital_kappa = orbital_kappa(orbital_index)
        current_orbital_two_j = orbital_two_j(orbital_index)
        current_orbital_weight = orbital_occupation_weights(orbital_index)

        ! Generate YA coefficients that do not require MCP output list.
        ! Computation distributed and then collected.

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

                ! <<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
                sumr = 0.0_real64
                do symmetry_block_index = 1, symmetry_block_count
                    do csf_index = mpi_rank + 1, csf_counts_by_block(symmetry_block_index), mpi_process_count
                        sumr = sumr + csf_density_matrix_element(extended_optimal_level,csf_index,csf_index, &
                            & symmetry_block_index)*direct_angular_coefficient(multipole_rank, &
                            & csf_index + csf_offsets_by_block(symmetry_block_index),orbital_index,&
                            partner_orbital_index)
                    end do
                end do
                ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                if (partner_orbital_index == orbital_index) then
                    ykab = 2.0_real64*sumr/current_orbital_weight
                else
                    ykab = sumr/current_orbital_weight
                end if

                ! The following IF has to be removed since YKAB
                ! is incomplete in multiprocessor case
                y_potential_coefficient_count = y_potential_coefficient_count + 1
                if (y_potential_coefficient_count > y_potential_coefficient_capacity ) then
                    if (y_potential_coefficient_capacity > 0) then
                        y_potential_coefficient_capacity = 2*y_potential_coefficient_capacity
                        call ralloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                            "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                        call ralloc(y_potential_coefficients, y_potential_coefficient_capacity, "y_potential_coefficients", &
                            "assemble_orbital_potential_coefficients")
                    else
                        y_potential_coefficient_capacity = 64
                        call alloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                            "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                        call alloc(y_potential_coefficients, y_potential_coefficient_capacity, "y_potential_coefficients", &
                            "assemble_orbital_potential_coefficients")
                    end if
                end if
                y_potential_coefficients(y_potential_coefficient_count) = ykab
                y_potential_coefficient_labels(y_potential_coefficient_count) = multipole_rank + packed_orbital_label_offset

            end do
        end do

        ! Generate XA coefficients that do not require MCP output list
        ! Computation distributed and then collected.

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

                ! <<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
                sumr = 0.0_real64
                do symmetry_block_index = 1, symmetry_block_count
                    do csf_index = mpi_rank + 1, csf_counts_by_block(symmetry_block_index), mpi_process_count
                        sumr = sumr + csf_density_matrix_element(extended_optimal_level,csf_index,csf_index, &
                            & symmetry_block_index)*exchange_angular_coefficient(multipole_rank, &
                            & csf_index + csf_offsets_by_block(symmetry_block_index),orbital_index,&
                            partner_orbital_index)
                    end do
                end do
                ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                xkab = sumr/current_orbital_weight

                ! The following IF has to be removed since YKAB
                ! is incomplete in multiprocessor case
                exchange_coefficient_count = exchange_coefficient_count + 1
                if (exchange_coefficient_count > exchange_coefficient_capacity) then
                    if (exchange_coefficient_capacity > 0) then
                        exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                        call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, "exchange_coefficient_labels", &
                            "assemble_orbital_potential_coefficients")
                    else
                        exchange_coefficient_capacity = 64
                        call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, "exchange_coefficient_labels", &
                            "assemble_orbital_potential_coefficients")
                    end if
                end if
                exchange_coefficients(exchange_coefficient_count) = xkab
                exchange_coefficient_labels(exchange_coefficient_count) = multipole_rank + packed_orbital_label_offset
            end do
        end do

        ! Replay the one-electron groups in block and label order.
        do jblock = 1, symmetry_block_count
            call begin_groups(jblock,31)
            do while (advance_group())
                ia = mod(lab,orbital_index_radix)
                partner_orbital_index = lab/orbital_index_radix
                if (ia /= orbital_index .and. partner_orbital_index /= orbital_index) then
                    cycle
                end if
                sum = normalized_sum(jblock,31)

                ! Put coefficients in the list. Since there is always
                ! (almost) some repeatence from different blocks, a check
                ! and merge is performed. This will significantly reduce
                ! the number of direct coefficients and calls to the Hartree Y and Z evaluator.

                ! Find the right counting parameter ***
                if (ia == orbital_index) then
                    ithis = partner_orbital_index
                else
                    ithis = ia
                end if

                ! Check it against the previously recorded ***
                if (jblock > 1) then
                    do i = 1, direct_coefficient_count
                        if (direct_orbital_indices(i) /= ithis) cycle     ! found, add the value
                        direct_coefficients(i) = direct_coefficients(i) + sum
                        exit
                    end do
                    if (i <= direct_coefficient_count) cycle
                end if

                ! Not found in the record, add an item ***
                direct_coefficient_count = direct_coefficient_count + 1
                if (direct_coefficient_count > direct_coefficient_capacity) then
                    if (direct_coefficient_capacity > 0) then
                        direct_coefficient_capacity = 2*direct_coefficient_capacity
                        call ralloc (direct_coefficients, direct_coefficient_capacity, "direct_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call ralloc (direct_orbital_indices, direct_coefficient_capacity, "direct_orbital_indices", &
                            "assemble_orbital_potential_coefficients")
                    else
                        direct_coefficient_capacity = 64
                        call alloc (direct_coefficients, direct_coefficient_capacity, "direct_coefficients", &
                            "assemble_orbital_potential_coefficients")
                        call alloc (direct_orbital_indices, direct_coefficient_capacity, "direct_orbital_indices", &
                            "assemble_orbital_potential_coefficients")
                    end if
                end if
                direct_coefficients(direct_coefficient_count) = sum                      ! print*, DA(NDCOF), ndcof, myid, 'myid'
                direct_orbital_indices(direct_coefficient_count) = ithis

            end do
        end do

        ! Preserve the k, block, label order used for YA and XA additions.
        do nfile = 32, 32 + max_mcp_file_multipole_rank
            multipole_rank = nfile - 32
            do jblock = 1, symmetry_block_count
                call begin_groups(jblock,nfile)
                do while (advance_group())
                    call decode_mcp_label(lab,.false.,decoded_integral_orbital_indices)
                    irank = count(decoded_integral_orbital_indices == orbital_index)
                    if (irank == 0) then
                            cycle
                    end if
                    sum = normalized_sum(jblock,nfile)

                    select case (irank)
                    case (1)

                        ! One matching index: exchange potential contribution

                        ! Similar to DA, find ithis ***
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

                        if (ithis == (-911)) error stop "ithis .EQ. -911"

                        ! Check ithis against the previously recorded ***
                        if (jblock > 1) then
                            do i = 1, exchange_coefficient_count
                                if (exchange_coefficient_labels(i) /= ithis) cycle
                                exchange_coefficients(i) = exchange_coefficients(i) + sum
                                exit
                            end do
                            if (i <= exchange_coefficient_count) cycle
                        end if

                        ! Not found in records, add an item ***
                        exchange_coefficient_count = exchange_coefficient_count + 1
                        if (exchange_coefficient_count > exchange_coefficient_capacity) then
                            if (exchange_coefficient_capacity > 0) then
                                exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                                call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                            else
                                exchange_coefficient_capacity = 64
                                call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                            end if
                        end if
                        exchange_coefficients(exchange_coefficient_count) = sum
                        exchange_coefficient_labels(exchange_coefficient_count) = ithis

                    case (2)

                        ! Two matching indices: either direct or exchange potential
                        ! contribution

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
                            !    Direct contribution
                            !
                            ! Find ithis ***
                            il = loc1 + 3
                            if (il > 4) il = il - 4
                            iyo2 = decoded_integral_orbital_indices(il)
                            il = loc1 + 1
                            if (il > 4) il = il - 4
                            iyo1 = decoded_integral_orbital_indices(il)
                            ithis = (iyo2*orbital_index_radix + iyo1)*orbital_index_radix + multipole_rank

                            ! Check it against the previously recorded ***
                            if (jblock > 1) then
                                do i = 1, y_potential_coefficient_count
                                    if (y_potential_coefficient_labels(i) /= ithis) cycle
                                    y_potential_coefficients(i) = y_potential_coefficients(i) + sum + sum
                                    exit
                                end do
                                if (i <= y_potential_coefficient_count) cycle
                            end if

                            ! Not found, add an item ***
                            y_potential_coefficient_count = y_potential_coefficient_count + 1
                            if (y_potential_coefficient_count > y_potential_coefficient_capacity) then
                                if (y_potential_coefficient_capacity > 0) then
                                    y_potential_coefficient_capacity = 2*y_potential_coefficient_capacity
                                    call ralloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    call ralloc(y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                else
                                    y_potential_coefficient_capacity = 64
                                    call alloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients 2b")
                                    call alloc(y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                end if
                            end if
                            y_potential_coefficients(y_potential_coefficient_count) = sum + sum
                            y_potential_coefficient_labels(y_potential_coefficient_count) = ithis

                        else
                            !    Exchange contribution
                            !
                            ! Find ithis ***
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

                            ! Check it against the previously recorded ***
                            if (jblock > 1) then
                                do i = 1, exchange_coefficient_count
                                    if (exchange_coefficient_labels(i) /= ithis) cycle
                                    exchange_coefficients(i) = exchange_coefficients(i) + sum + sum
                                    exit
                                end do
                                if (i <= exchange_coefficient_count) cycle
                            end if

                            ! Not found, add an item ***
                            exchange_coefficient_count = exchange_coefficient_count + 1
                            if (exchange_coefficient_count > exchange_coefficient_capacity) then
                                if (exchange_coefficient_capacity > 0) then
                                    exchange_coefficient_capacity = 2*exchange_coefficient_capacity
                                    call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                        "assemble_orbital_potential_coefficients")
                                    call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                        "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                                else
                                    exchange_coefficient_capacity = 64
                                    call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                        "assemble_orbital_potential_coefficients")
                                    call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                        "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                                end if
                            end if
                            exchange_coefficients(exchange_coefficient_count) = sum + sum
                            exchange_coefficient_labels(exchange_coefficient_count) = ithis

                        end if

                    case (3)

                        ! Three matching indices: direct and exchange potential contributions

                        ! Find ithis AND ithis2
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

                        if (ithis==(-911) .or. ithis2==(-911)) error stop "ithis2"

                        ! Check the previously recorded for YA
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
                                else
                                    y_potential_coefficient_capacity = 64
                                    call alloc (y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                        "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                    call alloc (y_potential_coefficients, y_potential_coefficient_capacity, &
                                        "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                                end if
                            end if
                            y_potential_coefficients(y_potential_coefficient_count) = sum + sum
                            y_potential_coefficient_labels(y_potential_coefficient_count) = ithis

                        end if

                        ! Check the previously recorded for XA
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
                                call ralloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call ralloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                            else
                                exchange_coefficient_capacity = 64
                                call alloc (exchange_coefficients, exchange_coefficient_capacity, "exchange_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                                call alloc (exchange_coefficient_labels, exchange_coefficient_capacity, &
                                    "exchange_coefficient_labels", "assemble_orbital_potential_coefficients")
                            end if
                        end if
                        exchange_coefficients(exchange_coefficient_count) = sum
                        exchange_coefficient_labels(exchange_coefficient_count) = ithis2

                    case (4)

                        ! Four matching indices: direct potential contribution

                        ! Find ithis AND ithis2
                        iyo2 = orbital_index
                        iyo1 = orbital_index
                        ithis = (iyo2*orbital_index_radix + iyo1)*orbital_index_radix + multipole_rank

                        ! Check the previously recorded for YA
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
                                call ralloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                    "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients")
                                call ralloc(y_potential_coefficients, y_potential_coefficient_capacity, &
                                    "y_potential_coefficients", "assemble_orbital_potential_coefficients")
                            else
                                y_potential_coefficient_capacity = 64
                                call alloc(y_potential_coefficient_labels, y_potential_coefficient_capacity, &
                                    "y_potential_coefficient_labels", "assemble_orbital_potential_coefficients 4b")
                                call alloc(y_potential_coefficients, y_potential_coefficient_capacity, "y_potential_coefficients", &
                                    "assemble_orbital_potential_coefficients")
                            end if
                        end if
                        y_potential_coefficients(y_potential_coefficient_count) = 4.0_real64*sum
                        y_potential_coefficient_labels(y_potential_coefficient_count) = ithis

                    end select

                end do
                end do
        end do
        if (verification_enabled) then
            do i = 1, direct_coefficient_count
                call log_bits("DA",orbital_index,0,i,int(direct_orbital_indices(i),int64),direct_coefficients(i))
            end do
            do i = 1, y_potential_coefficient_count
                call log_bits("YA",orbital_index,0,i,int(y_potential_coefficient_labels(i),int64),y_potential_coefficients(i))
            end do
            do i = 1, exchange_coefficient_count
                call log_bits("XA",orbital_index,0,i,int(exchange_coefficient_labels(i),int64),exchange_coefficients(i))
            end do
        end if
        call profile_stage("potential_coefficients",orbital_index,0,started)

    contains

        subroutine begin_groups(block, channel)
            integer, intent(in) :: block, channel
            csf_count = csf_counts_by_block(block)
            call contraction_bounds(block,channel,first,last)
            group_index = first-1
        end subroutine begin_groups

        logical function advance_group() result(found)
            group_index = group_index + 1
            found = group_index <= last
            if (found) lab = contracted_labels(group_index)
        end function advance_group

        real(real64) function normalized_sum(block, channel) result(value)
            integer, intent(in) :: block, channel
            value = contracted_values(group_index)
            if (verification_enabled) call log_bits("group",orbital_index,block,channel,int(lab,int64),value)
            value = 0.5_real64*value/current_orbital_weight
        end function normalized_sum
    end subroutine assemble_orbital_potential_coefficients
end module assemble_orbital_potential_coefficients_interface
