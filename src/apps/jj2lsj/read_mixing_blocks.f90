! Load level energies and mixing coefficients from the selected block-structured mixing file.
module read_mixing_blocks_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mixing_blocks
contains

    ! Reads mixing coefficient file from block-structured format
    !
    ! Note:
    ! eav is not compatible with the non-block version if some blocks
    ! were not diagonalized
    !
    ! This is a modified version of cvtmix.f
    !
    ! Written by Per Jonsson,                           September 2003
    ! Modified by G. Gaigalas,. May 2011
    !
    ! Modified by C. Cychen. 2021
    subroutine read_mixing_blocks(state_name, ci_input_mode)
        ! Module dependencies
        use mixing_file_format, only: mixing_file, open_mixing_input, next_mixing_block, read_mixing_levels, read_mixing_vectors
        use csf_descriptor_reader, only: loaded_csf_basis_id
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use calculation_state_files, only: mixing_input_path
        use atomic_calculation_state, only: electron_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count, orbital_count
        use csf_quantum_queries, only: csf_total_j_code, csf_parity_code
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use transition_basis_blocks,            only: transition_level_counts_by_block, transition_csf_counts_by_block, &
            transition_two_j_by_block, transition_block_count, cumulative_csf_counts
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: mix_electron_count, mix_orbital_count
        integer , intent(in) :: ci_input_mode
        character(len=24), intent(in) :: state_name
        integer :: total_csf_count, nvectot, nblock, nvecpat, &
            ncfpat, neavsum, jb, nb, ncfblk, nevblk, iatjp, iaspa
        real(real64) :: eavsum
        integer(int64) :: nvecsizpat, ncftot_i8
        character(len=256) :: file_path
        type(mixing_file) :: input_file

        file_path = mixing_input_path(state_name, ci_input_mode == 0)
        call open_mixing_input(input_file, file_path, loaded_csf_basis_id)
        mix_electron_count = input_file%electron_count
        mix_orbital_count = input_file%orbital_count
        total_csf_count = input_file%csf_count
        nvectot = input_file%level_count
        nblock = input_file%block_count

        if (mix_electron_count /= electron_count .or. total_csf_count /= csf_count .or. &
            mix_orbital_count < orbital_count .or. nblock /= transition_block_count) then
            error stop "JJ2LSJ mixing header does not match descriptor basis"
        end if
        if (nvectot < 1 .or. nvectot > csf_count) error stop "JJ2LSJ invalid mixing level count"
        write (*, *) "   nelec  = ", electron_count
        write (*, *) "   ncftot = ", total_csf_count
        write (*, *) "   nw     = ", mix_orbital_count
        write (*, *) "   nblock = ", nblock
        write (*, *)
        ncftot_i8 = total_csf_count

        ! Allocate memory for old format data

        allocate (shifted_eigenvalues(nvectot))
        allocate (csf_mixing_coefficients(ncftot_i8*nvectot))
        allocate (active_level_indices(nvectot))
        allocate (level_two_j_plus_one(nvectot))
        allocate (level_parity_signs(nvectot))

        ! Initialize mixing coefficients to zero; others are fine
        csf_mixing_coefficients(:nvectot*ncftot_i8) = 0.0_real64

        ! Initialize counters and sum registers
        !
        !    nvecpat:    total number of eigenstates of the previous blocks
        !    ncfpat:     total number of CSF of the previous blocks
        !    nvecsizpat: vector size of the previous blocks
        !    eavsum:     sum of diagonal elements of the previous blocks where
        !                at least one eigenstate is calculated
        !    neavsum:    total number CSF contributing to eavsum

        nvecpat = 0
        ncfpat = 0
        nvecsizpat = 0
        neavsum = 0
        eavsum = 0.0_real64

        write (*, *) "  block     ncf     nev    2j+1  parity"
        do jb = 1, nblock

            call next_mixing_block(input_file, nb, ncfblk, nevblk, iatjp, iaspa)
            write (*, "(5I8)") nb, ncfblk, nevblk, iatjp, iaspa
            if (jb /= nb .or. ncfblk /= cumulative_csf_counts(jb) - ncfpat) then
                error stop "JJ2LSJ mixing block does not match descriptor CSF partition"
            end if
            if (nevblk < 0 .or. nevblk > ncfblk .or. nevblk > nvectot - nvecpat) then
                error stop "JJ2LSJ invalid mixing block level count"
            end if
            if (nevblk > 0) then
                if (iatjp /= csf_total_j_code(ncfpat+1) .or. iaspa /= csf_parity_code(ncfpat+1)) then
                    error stop "JJ2LSJ mixing block does not match descriptor symmetry"
                end if
            end if
            transition_level_counts_by_block(jb) = nevblk
            transition_csf_counts_by_block(jb) = ncfblk
            transition_two_j_by_block(jb) = iatjp - 1

            if (nevblk > 0) then

                call read_mixing_levels(input_file, active_level_indices(nvecpat+1:nvecpat+nevblk), &
                    hamiltonian_energy_shift, shifted_eigenvalues(nvecpat+1:nvecpat+nevblk))
                ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
                level_two_j_plus_one(nvecpat+1:nevblk+nvecpat) = iatjp
                level_parity_signs(nvecpat+1:nevblk+nvecpat) = iaspa

                ! ...Construct the true energy by adding up the average
                shifted_eigenvalues(nvecpat+1:nevblk+nvecpat) = shifted_eigenvalues(nvecpat+1:nevblk+nvecpat) + &
                    hamiltonian_energy_shift
                ! ...For overal (all blocks) average energy
                eavsum = eavsum + hamiltonian_energy_shift*ncfblk
                neavsum = neavsum + ncfblk

                call read_mixing_vectors(input_file, csf_mixing_coefficients(nvecsizpat+ncfpat+1:), total_csf_count)
            end if
            nvecpat = nvecpat + nevblk
            ncfpat = ncfpat + ncfblk
            nvecsizpat = nvecsizpat + nevblk*ncftot_i8
        end do

        ! ...Here eav is the average energy of the blocks where at least
        !    one eigenstate is calculated. It is not the averge of the
        !    total Hamiltonian.

        if (nvecpat /= nvectot) error stop "JJ2LSJ mixing blocks disagree with total level count"
        hamiltonian_energy_shift = eavsum/neavsum

        if (total_csf_count /= neavsum) write (output_unit, *) &
            "Not all blocks are diagonalized --- Average E ", "not correct"

        ! ...Substrct the overal average energy
        shifted_eigenvalues(:nvectot) = shifted_eigenvalues(:nvectot) - hamiltonian_energy_shift
        close(input_file%unit)
        active_level_count = nvectot
        return
    end subroutine read_mixing_blocks
end module read_mixing_blocks_interface
