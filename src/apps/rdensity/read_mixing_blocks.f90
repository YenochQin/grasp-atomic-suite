! Load level energies and mixing coefficients from the selected block-structured mixing file.
module read_mixing_blocks_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mixing_blocks
contains

    !    Reads mixing coefficient file from block-structured format
    !
    ! Note:
    !     eav is not compatible with the non-block version if some blocks
    !     were not diagonalized
    !
    !     This is a modified version of cvtmix.f
    ! written by Per Jonsson, September 2003
    subroutine read_mixing_blocks(state_name, ci_input_mode)
        ! Module dependencies
        use mixing_file_format, only: mixing_file, open_mixing_input, next_mixing_block, read_mixing_levels, read_mixing_vectors
        use mixing_file_format, only: validate_mixing_basis
        use csf_descriptor_reader, only: loaded_csf_basis_id
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use calculation_state_files, only: mixing_input_path
        use array_memory_management, only: alloc
        use density_parallel, only: density_fail, density_abort
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use, intrinsic :: iso_fortran_env, only: int64
        use atomic_calculation_state, only: electron_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis, only: csf_count, orbital_count
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        implicit none (type, external)

        ! Arguments and local declarations
        integer, intent(in)           :: ci_input_mode
        character(len=24), intent(in) :: state_name
        integer :: total_csf_count, nvectot, nblock, i, nvecpat, &
            ncfpat, nvecsizpat, neavsum, jb, nb, ncfblk, nevblk, iatjp, iaspa
        real(real64) :: eavsum
        character :: file_path*256
        type(mixing_file) :: input_file

        file_path = mixing_input_path(state_name, ci_input_mode == 0)
        call open_mixing_input(input_file, file_path, loaded_csf_basis_id)
        call validate_mixing_basis(input_file, density_abort)
        total_csf_count = input_file%csf_count
        nvectot = input_file%level_count
        nblock = input_file%block_count
        if (int(total_csf_count,int64)*nvectot > huge(csf_count)) &
            call density_fail("Mixing vector size exceeds integer index capacity")
        write (*, *) "   nelec  = ", electron_count
        write (*, *) "   ncftot = ", total_csf_count
        write (*, *) "   nw     = ", orbital_count
        write (*, *) "   nblock = ", nblock
        write (*, *)

        ! Allocate memory for old format data

        allocate(shifted_eigenvalues(nvectot))
        allocate(csf_mixing_coefficients(total_csf_count*nvectot))
        call alloc (active_level_indices, nvectot, "active_level_indices", "read_mixing_blocks")
        call alloc (level_two_j_plus_one, nvectot, "level_two_j_plus_one", "read_mixing_blocks")
        call alloc (level_parity_signs, nvectot, "level_parity_signs", "read_mixing_blocks")

        ! Initialize mixing coefficients to zero; others are fine
        csf_mixing_coefficients(:nvectot*total_csf_count) = 0.0_real64

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
            nvecsizpat = nvecsizpat + nevblk*total_csf_count

        end do

        ! ...Here eav is the average energy of the blocks where at least
        !    one eigenstate is calculated. It is not the averge of the
        !    total Hamiltonian.

        if (.not.all(ieee_is_finite(csf_mixing_coefficients))) call density_fail("Non-finite mixing coefficients")
        do i=1,nvectot
            if (abs(sum(csf_mixing_coefficients((i-1)*total_csf_count+1:i*total_csf_count)**2)-1.0_real64) > 1.0e-6_real64) &
                call density_fail("Mixing eigenvector is not normalized")
        end do
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
