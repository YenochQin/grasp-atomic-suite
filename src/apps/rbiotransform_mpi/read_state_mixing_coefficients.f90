!***********************************************************************
!                                                                      *
subroutine read_state_mixing_coefficients(name, inpci, iblk)
    ! Module dependencies
    use calculation_state_files, only: mixing_input_path
    use calculation_file_units, only: mixing_input_unit
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Open, check, load data from and close the rscf.mix file.           *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, LENGTH, OPENFL.                        *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 25 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use array_memory_management, only: alloc
    use atomic_calculation_state, only: electron_count
    use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
    use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
        nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: inpci
    integer, intent(in) :: iblk
    character(len=*), intent(in) :: name

    ! Local declarations
    integer :: j, io_status, nb, iatjp, iaspa, i
    character :: mixing_file_signature*6
    !-----------------------------------------------

    if (iblk == 1) then
        open(newunit=mixing_input_unit, file=mixing_input_path(name, inpci == 0), &
            form="UNFORMATTED", status="OLD", position="asis")

        read (mixing_input_unit, iostat=io_status) mixing_file_signature
        if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
            write (output_unit, *)                                              &
                "File", mixing_input_unit, "Not a GRASP MIXing Coefficients File;"
            close(mixing_input_unit)
            stop
        end if
        !
        read (mixing_input_unit) nelectot, ncftot, nwtot, nvectot, nvecsiztot, nblock1
        if (electron_count/=nelectot .or. orbital_count/=nwtot) then
            !     :    (NCF .NE. NCFT) .OR.
            write (output_unit, *) "File", mixing_input_unit, "is not an"
            write (output_unit, *) " appropriate to Coefficients file"
            close(mixing_input_unit)
            stop
        end if
    end if
    !
    !   Load data from the  rscf.mix  file
    !
    write (output_unit, *) "Loading MIXing Coefficients File ..."
    !
    read (mixing_input_unit) nb, csf_count, nvec, iatjp, iaspa
    call alloc (shifted_eigenvalues, nvec, "EVAL", "GETMIX")
    call alloc (csf_mixing_coefficients, csf_count*nvec, "EVEC", "GETMIX")
    call alloc (ivec, nvec, "IVEC", "GETMIX")
    call alloc (iatjpo, nvec, "IATJPO", "GETMIX")
    call alloc (iaspar, nvec, "IASPAR", "GETMIX")
    !
    !   These arrays are deallocated in mcp
    !
    read (mixing_input_unit) (ivec(i),i=1,nvec)
    !     READ (mixing_input_unit) (IATJPO(I),IASPAR(I),I = 1,NVEC)
    iatjpo(:nvec) = iatjp
    iaspar(:nvec) = iaspa
    read (mixing_input_unit) hamiltonian_energy_shift, (shifted_eigenvalues(i),i=1,nvec)
    read (mixing_input_unit) ((csf_mixing_coefficients(i + (j - 1)*csf_count),i=1,csf_count),j=1,nvec)
    !
    write (output_unit, *) " ... load complete;"
    !
    !   Close the  rscf.mix  file
    !
    !     CLOSE (mixing_input_unit)
    !
    return
end subroutine read_state_mixing_coefficients
