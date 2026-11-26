! Build one symmetry block of the CI Hamiltonian from angular and radial integrals.
module assemble_hamiltonian_block_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_hamiltonian_block
contains

    ! Modified by G. Gaigalas. Last revision: 07 Sep 2021
    !    Spin-angular coefficiants have been put on memory
    subroutine assemble_hamiltonian_block(symmetry_block_index, mpi_rank, mpi_process_count, total_contribution_count, &
        & total_secondary_contribution_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use grasp_parameters,    only: orbital_index_radix
        use array_memory_management, only: alloc, dalloc
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_two_j
        use symmetry_block_offsets, only: csf_offsets_by_block
        use rmcdhf_mcp_memory,     only: indx, coeff
        use csf_quantum_queries, only: csf_occupation
        use one_electron_radial_integral, only: integrate_one_electron_radial

        ! Procedure interfaces
        use direct_angular_coefficient_interface, only: direct_angular_coefficient

        ! Module dependencies
        use slater_radial_integral, only: calculate_slater_radial_integral

        ! Procedure interfaces
        use exchange_angular_coefficient_interface, only: exchange_angular_coefficient
        use load_mcp_records_into_memory_interface, only: load_mcp_records_into_memory
        use read_mcp_structure_into_memory_interface, only: read_mcp_structure_into_memory
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: symmetry_block_index
        integer , intent(in) :: mpi_rank
        integer , intent(in) :: mpi_process_count
        integer(int64) :: total_contribution_count(20), total_secondary_contribution_count(20)

        ! Local constants
        character(len=*), parameter :: routine_name = "assemble_hamiltonian_block"

        ! Local declarations
        integer(int64) :: diagonal_element_index
        integer :: block_csf_offset, orbital_index, csf_index, itmp, km, k, k0, ib, first_orbital_two_j, &
            second_orbital_two_j, kmin, kmax, ndim, packed_orbital_label, contribution_count, file_csf_column, i, loc, &
            mcp_channel, id, ic
        integer(int64) :: lod      !  INTEGER*8 version of LOC
        real(real64) :: diaa, coef, f0aa, fkaa, f0ab, gkab, radial_integral_value
        logical :: set
        ! !*** needs EMT, IENDC.
        !
        ! !*** needs ncf, nw
        !
        ! !POINTER (PICLMN,ICLMNdum)
        block_csf_offset = csf_offsets_by_block(symmetry_block_index)
        ! Accumulate diagonal terms that do not require MCP coefficients

        ! Piece involving I(a,a) integrals

        do orbital_index = 1, orbital_count
            set = .false.
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                itmp = csf_occupation(orbital_index,csf_index + block_csf_offset)
                if (itmp <= 0) cycle
                ! Occupation number not zero ...
                if (.not.set) then
                    diaa = integrate_one_electron_radial(orbital_index,orbital_index,0)
                    set = .true.
                end if
                ! IDIAG = IENDC(IR-1)+1
                diagonal_element_index = stored_segment_end_offsets(csf_index)                    ! lower-triangle-by-rows mode
                hamiltonian_values(diagonal_element_index) = hamiltonian_values(diagonal_element_index) + itmp*diaa
            end do
        end do

        !                  0
        ! Piece involving F (a,a) integrals

        do orbital_index = 1, orbital_count
            set = .false.
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                coef = direct_angular_coefficient(0,csf_index + block_csf_offset,orbital_index,orbital_index)
                if (coef == 0.0_real64) cycle
                ! Angular coefficient not zero ...
                if (.not.set) then
                    f0aa = calculate_slater_radial_integral(orbital_index,orbital_index,orbital_index,orbital_index,0)
                    set = .true.
                end if
                ! IDIAG = IENDC(IR-1)+1
                diagonal_element_index = stored_segment_end_offsets(csf_index)
                hamiltonian_values(diagonal_element_index) = hamiltonian_values(diagonal_element_index) + coef*f0aa
            end do
        end do

        !                  k
        ! Piece involving F (a,a) integrals

        km = 0
        k = 0
        do
            k = k + 2
            do orbital_index = 1, orbital_count
                k0 = orbital_two_j(orbital_index) - 1
                km = max0(k0,km)
                if (k > k0) cycle
                set = .false.

                do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                    coef = direct_angular_coefficient(k,csf_index + block_csf_offset,orbital_index,orbital_index)
                    if (coef == 0.0_real64) cycle
                    if (.not.set) then
                        fkaa = calculate_slater_radial_integral(orbital_index,orbital_index,orbital_index,orbital_index,k)
                        set = .true.
                    end if
                    ! IDIAG = IENDC(IR-1)+1
                    diagonal_element_index = stored_segment_end_offsets(csf_index)
                    hamiltonian_values(diagonal_element_index) = hamiltonian_values(diagonal_element_index) + coef*fkaa
                end do
            end do
            if (k >= km) exit
        end do

        !                  0
        ! Piece involving F (a,b) integrals

        do orbital_index = 1, orbital_count - 1
            do ib = orbital_index + 1, orbital_count
                set = .false.
                do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                    coef = direct_angular_coefficient(0,csf_index + block_csf_offset,orbital_index,ib)
                    if (coef == 0.0_real64) cycle
                    if (.not.set) then
                        f0ab = calculate_slater_radial_integral(orbital_index,ib,orbital_index,ib,0)
                        set = .true.
                    end if
                    ! IDIAG = IENDC(IR-1)+1
                    diagonal_element_index = stored_segment_end_offsets(csf_index)
                    hamiltonian_values(diagonal_element_index) = hamiltonian_values(diagonal_element_index) + coef*f0ab
                end do
            end do
        end do

        !                  k
        ! Piece involving G (a,b) integrals

        km = 0
        k = -1
        do
            k = k + 1
            do orbital_index = 1, orbital_count - 1
                first_orbital_two_j = orbital_two_j(orbital_index)
                do ib = orbital_index + 1, orbital_count
                    second_orbital_two_j = orbital_two_j(ib)
                    set = .false.
                    if (orbital_kappa(orbital_index)*orbital_kappa(ib) > 0) then
                        kmin = abs((first_orbital_two_j - second_orbital_two_j)/2)
                    else
                        kmin = abs((first_orbital_two_j - second_orbital_two_j)/2) + 1
                    end if
                    if (mod(k - kmin,2) /= 0) cycle

                    kmax = (first_orbital_two_j + second_orbital_two_j)/2
                    km = max0(kmax,km)
                    if (k<kmin .or. k>kmax) cycle

                    do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                        coef = exchange_angular_coefficient(k,csf_index + block_csf_offset,orbital_index,ib)
                        if (coef == 0.0_real64) cycle
                        if (.not.set) then
                            gkab = calculate_slater_radial_integral(orbital_index,ib,ib,orbital_index,k)
                            set = .true.
                        end if
                        ! IDIAG = IENDC(IR-1)+1
                        diagonal_element_index = stored_segment_end_offsets(csf_index)
                        hamiltonian_values(diagonal_element_index) = hamiltonian_values(diagonal_element_index) + coef*gkab
                    end do
                end do
            end do
            if (k >= km) exit
        end do

        ! Local storage for reading mcpXXX files

        ndim = 1
        call alloc (coeff, ndim, "COEFF", "assemble_hamiltonian_block" )
        call alloc (indx, ndim, "INDX", "assemble_hamiltonian_block")

        ! Loop over non-zero labels which have non-zero elements

        total_secondary_contribution_count(1) = total_secondary_contribution_count(1) + 1
        ! Cached labels are memory accesses; there is no IOSTAT result.
        call load_mcp_records_into_memory(31,total_secondary_contribution_count(1),packed_orbital_label, contribution_count)
        do while(packed_orbital_label/=0 .or. contribution_count/=0)

            ! decode the label of I(ab)
            orbital_index = mod(packed_orbital_label,orbital_index_radix)
            ib = packed_orbital_label/orbital_index_radix

            ! Compute radial integral I(ab)
            radial_integral_value = integrate_one_electron_radial(orbital_index,ib,0)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            if (contribution_count > ndim) then
                call dalloc (coeff, "COEFF", "assemble_hamiltonian_block")
                call dalloc (indx, "INDX", "assemble_hamiltonian_block")
                ndim = contribution_count
                call alloc (coeff, ndim, "COEFF", "assemble_hamiltonian_block")
                call alloc (indx, ndim, "INDX", "assemble_hamiltonian_block")
            end if
            call read_mcp_structure_into_memory(31,total_contribution_count(1),contribution_count,file_csf_column)
            total_contribution_count(1) = total_contribution_count(1) + contribution_count

            ! Store all the contributions from this integral
            do i = 1, contribution_count
                loc = indx(i)
                lod = loc   !  (Convert type)
                if (lod > stored_matrix_element_count) then
                    write (output_unit, *) "  Error in computing 1-e contribution"
                    write (output_unit, *) "  LOC = ", loc , "  NELMNT = ", stored_matrix_element_count
                    stop
                end if
                hamiltonian_values(loc) = hamiltonian_values(loc) + radial_integral_value*coeff(i)
            end do

            total_secondary_contribution_count(1) = total_secondary_contribution_count(1) + 1
            call load_mcp_records_into_memory(31,total_secondary_contribution_count(1),packed_orbital_label, contribution_count)
        end do

        do mcp_channel = 32, 32 + max_mcp_file_multipole_rank
            k = mcp_channel - 32

            ! Loop over non-zero labels which have non-zero elements

            total_secondary_contribution_count(mcp_channel-30) = total_secondary_contribution_count(mcp_channel-30) + 1
            call load_mcp_records_into_memory(mcp_channel, total_secondary_contribution_count(mcp_channel-30), &
                    packed_orbital_label, &
                & contribution_count)
            do while(packed_orbital_label/=0 .or. contribution_count/=0)

                !                          k
                ! decode the label of R (abcd)
                id = mod(packed_orbital_label,orbital_index_radix)
                packed_orbital_label = packed_orbital_label/orbital_index_radix
                ib = mod(packed_orbital_label,orbital_index_radix)
                packed_orbital_label = packed_orbital_label/orbital_index_radix
                ic = mod(packed_orbital_label,orbital_index_radix)
                orbital_index = packed_orbital_label/orbital_index_radix

                ! Compute radial integral
                radial_integral_value = calculate_slater_radial_integral(orbital_index,ib,ic,id,k)

                ! Read column index, sparse matrix index, and coefficient
                ! for all contributions from this integral.
                if (contribution_count > ndim) then
                    call dalloc (coeff, "COEFF", "assemble_hamiltonian_block")
                    call dalloc (indx, "INDX", "assemble_hamiltonian_block")
                    ndim = contribution_count
                    call alloc (coeff, ndim, "COEFF", "assemble_hamiltonian_block")
                    call alloc (indx, ndim, "INDX", "assemble_hamiltonian_block")
                end if
                call read_mcp_structure_into_memory(mcp_channel, total_contribution_count(mcp_channel-30), &
                    contribution_count, file_csf_column)
                total_contribution_count(mcp_channel-30) = total_contribution_count(mcp_channel-30) + contribution_count

                ! Store all the contributions from this integral
                do i = 1, contribution_count
                    loc = indx(i)
                    lod = loc
                    if (lod > stored_matrix_element_count) then   !! NELMENT (from Hmat_C) is INT*8
                        write (output_unit, *) "  Error in computing 2-e contribution"
                        write (output_unit, *) "  LOC= ", loc , "  NELMNT = ", stored_matrix_element_count
                        stop
                    end if
                    hamiltonian_values(loc) = hamiltonian_values(loc) + radial_integral_value*coeff(i)
                end do

                total_secondary_contribution_count(mcp_channel-30) = total_secondary_contribution_count(mcp_channel-30) + 1
                call load_mcp_records_into_memory(mcp_channel, total_secondary_contribution_count(mcp_channel-30), &
                    packed_orbital_label, &
                    & contribution_count)
            end do
        end do

        ! Deallocate local storage

        call dalloc (coeff, "COEFF", "assemble_hamiltonian_block")
        call dalloc (indx, "INDX", "assemble_hamiltonian_block")

        return
    end subroutine assemble_hamiltonian_block
end module assemble_hamiltonian_block_interface
