! Build one symmetry block of the CI Hamiltonian from angular and radial integrals.
module assemble_hamiltonian_block_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_hamiltonian_block
contains

    subroutine assemble_hamiltonian_block(symmetry_block_index, mpi_rank, mpi_process_count)
        ! Module dependencies
        use mpi_indexed_mcp, only: mcp_segment_size, decode_mcp_label
        use rmcdhf_mcp_io, only: mcp_cursor, mcp_buffer, start_channel, next_group, read_segment, profile_channel
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_two_j
        use symmetry_block_offsets, only: csf_offsets_by_block
        use csf_quantum_queries, only: csf_occupation
        use one_electron_radial_integral, only: integrate_one_electron_radial

        ! Procedure interfaces
        use direct_angular_coefficient_interface, only: direct_angular_coefficient

        ! Module dependencies
        use slater_radial_integral, only: calculate_slater_radial_integral

        ! Procedure interfaces
        use exchange_angular_coefficient_interface, only: exchange_angular_coefficient
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: symmetry_block_index
        integer , intent(in) :: mpi_rank
        integer , intent(in) :: mpi_process_count

        integer(int64) :: diagonal_element_index
        integer :: block_csf_offset, orbital_index, csf_index, itmp, km, k, k0, ib, &
            first_orbital_two_j, second_orbital_two_j, kmin, kmax, i, loc, mcp_channel, length, orbitals(4)
        integer(int64) :: first
        real(real64) :: diaa, coef, f0aa, fkaa, f0ab, gkab, radial_integral_value
        logical :: set
        type(mcp_cursor) :: cursor
        type(mcp_buffer) :: buffer

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

        ! Read the groups in their original order and use the current radial integrals.
        do mcp_channel = 31, 32 + max_mcp_file_multipole_rank
            call start_channel(cursor,symmetry_block_index,mcp_channel)
            do while (next_group(cursor))
                call decode_mcp_label(int(cursor%group%label),mcp_channel==31,orbitals)
                if (mcp_channel == 31) then
                    radial_integral_value = integrate_one_electron_radial(orbitals(1),orbitals(2),0)
                else
                    radial_integral_value = calculate_slater_radial_integral(orbitals(1),orbitals(2),orbitals(3), &
                        orbitals(4),mcp_channel-32)
                end if
                do first = 1, cursor%group%count, mcp_segment_size
                    length = int(min(int(mcp_segment_size,int64),cursor%group%count-first+1))
                    call read_segment(cursor,buffer,first,length,stored_matrix_element_count)
                    do i = 1, length
                        loc = buffer%positions(i)
                        hamiltonian_values(loc) = hamiltonian_values(loc) + radial_integral_value*buffer%values(i)
                    end do
                end do
            end do
            call profile_channel("matrix",0,cursor)
        end do
    end subroutine assemble_hamiltonian_block
end module assemble_hamiltonian_block_interface
