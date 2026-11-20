!***********************************************************************
!
subroutine read_mixing_blocks(state_name, ci_input_mode)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !
    !    Reads mixing coefficient file from block-structured format
    !
    ! Note:
    !     eav is not compatible with the non-block version if some blocks
    !     were not diagonalized
    !
    !     This is a modified version of cvtmix.f
    !     written by Per Jonsson, September 2003
    !
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use calculation_state_files, only: mixing_input_path
    use array_memory_management, only: alloc
    use atomic_calculation_state, only: electron_count
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use orbital_csf_basis, only: csf_count, orbital_count
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments and local declarations
    integer :: mixing_input_unit
    integer, intent(in)           :: ci_input_mode
    character(len=24), intent(in) :: state_name
    integer :: open_status, io_status, total_csf_count, nvectot, vector_size, nblock, i, nvecpat, &
        ncfpat, nvecsizpat, neavsum, jb, nb, ncfblk, nevblk, iatjp, iaspa, j
    real(real64) :: eavsum
    character :: file_path*256, form*11, mixing_file_signature*6, status*3
    !-----------------------------------------------
    !
    !   The  .mix  file is UNFORMATTED; it must exist
    !
    file_path = mixing_input_path(state_name, ci_input_mode == 0)
    form = "UNFORMATTED"
    status = "OLD"
    !
    call openfl (mixing_input_unit, file_path, form, status, open_status)
    if (open_status == 1) then
        write (stderr_unit, *) "Error when opening", file_path
        stop
    end if
    !
    !   Check the header of the file; if not as expected, try again
    !
    read (mixing_input_unit, iostat=io_status) mixing_file_signature
    if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
        write (stderr_unit, *) "Not a GRASP92 MIXing Coefficients File;"
        close(mixing_input_unit)
        stop
    end if

    read (mixing_input_unit) electron_count, total_csf_count, orbital_count, nvectot, vector_size, nblock
    write (*, *) "   nelec  = ", electron_count
    write (*, *) "   ncftot = ", total_csf_count
    write (*, *) "   nw     = ", orbital_count
    write (*, *) "   nblock = ", nblock
    write (*, *)

    !***********************************************************************
    ! Allocate memory for old format data
    !***********************************************************************

    call alloc (shifted_eigenvalues, nvectot, "EVAL", "GETMIXBLOCK")
    call alloc (csf_mixing_coefficients, total_csf_count*nvectot, "EVEC", "GETMIXBLOCK")
    call alloc (ivec, nvectot, "IVEC", "GETMIXBLOCK")
    call alloc (iatjpo, nvectot, "IATJPO", "GETMIXBLOCK")
    call alloc (iaspar, nvectot, "IASPAR", "GETMIXBLOCK")

    !***********************************************************************
    ! Initialize mixing coefficients to zero; others are fine
    !***********************************************************************
    csf_mixing_coefficients(:nvectot*total_csf_count) = 0.0_real64

    !***********************************************************************
    ! Initialize counters and sum registers
    !
    !    nvecpat:    total number of eigenstates of the previous blocks
    !    ncfpat:     total number of CSF of the previous blocks
    !    nvecsizpat: vector size of the previous blocks
    !    eavsum:     sum of diagonal elements of the previous blocks where
    !                at least one eigenstate is calculated
    !    neavsum:    total number CSF contributing to eavsum
    !***********************************************************************

    nvecpat = 0
    ncfpat = 0
    nvecsizpat = 0
    neavsum = 0
    eavsum = 0.0_real64

    write (*, *) "  block     ncf     nev    2j+1  parity"
    do jb = 1, nblock

        read (mixing_input_unit) nb, ncfblk, nevblk, iatjp, iaspa
        write (*, "(5I8)") nb, ncfblk, nevblk, iatjp, iaspa
        if (jb /= nb) stop "jb .NE. nb"

        if (nevblk > 0) then

            read (mixing_input_unit) (ivec(nvecpat + i),i=1,nevblk)
            ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
            iatjpo(nvecpat+1:nevblk+nvecpat) = iatjp
            iaspar(nvecpat+1:nevblk+nvecpat) = iaspa

            read (mixing_input_unit) hamiltonian_energy_shift, (shifted_eigenvalues(nvecpat+i),i=1,nevblk)

            !           ...Construct the true energy by adding up the average
            shifted_eigenvalues(nvecpat+1:nevblk+nvecpat) = shifted_eigenvalues(nvecpat+1:nevblk+nvecpat) + &
                hamiltonian_energy_shift
            !           ...For overal (all blocks) average energy
            eavsum = eavsum + hamiltonian_energy_shift*ncfblk
            neavsum = neavsum + ncfblk

            read (mixing_input_unit) ((csf_mixing_coefficients(nvecsizpat+ncfpat+i+(j-1)*total_csf_count),i=1,ncfblk),j=1,&
                nevblk)
        end if

        nvecpat = nvecpat + nevblk
        ncfpat = ncfpat + ncfblk
        nvecsizpat = nvecsizpat + nevblk*total_csf_count

    end do

    !     ...Here eav is the average energy of the blocks where at least
    !        one eigenstate is calculated. It is not the averge of the
    !        total Hamiltonian.

    hamiltonian_energy_shift = eavsum/neavsum

    if (total_csf_count /= neavsum) write (output_unit, *) &
        "Not all blocks are diagonalized --- Average E ", "not correct"

    !     ...Substrct the overal average energy
    shifted_eigenvalues(:nvectot) = shifted_eigenvalues(:nvectot) - hamiltonian_energy_shift

    close(mixing_input_unit)

    csf_count = total_csf_count
    nvec = nvectot

    return
end subroutine read_mixing_blocks
