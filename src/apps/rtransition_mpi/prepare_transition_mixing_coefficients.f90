module prepare_transition_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: prepare_transition_mixing_coefficients
contains

    !************************************************************************
    subroutine prepare_transition_mixing_coefficients(name, inpci)
        ! Module dependencies
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: copied_mixing_unit
        integer :: original_mixing_unit
        integer, intent(in) :: inpci
        character, intent(inout) :: name(2)*128
        integer :: io_status, nelec, total_csf_count, nw, nvectot, vector_size, nblock, iblk, ib, &
            ncf, nvec, iatjp, iaspa, i, j
        real(real64), pointer :: eval(:), evec(:)
        integer, pointer :: ivec(:)
        real(real64) :: hamiltonian_energy_shift
        character :: mixing_file_signature*6
        !-----------------------------------------------



        name(2) = trim(name(1))//"_CP"
        if (inpci == 0) then
            open(newunit=copied_mixing_unit, file=trim(name(2))//".cbm", form="UNFORMATTED", status=&
                "UNKNOWN", position="asis")
        else
            open(newunit=copied_mixing_unit, file=trim(name(2))//".bm", form="UNFORMATTED", status=&
                "UNKNOWN", position="asis")
        end if
        if (inpci == 0) then
            open(newunit=original_mixing_unit, file=trim(name(1))//".cbm", form="UNFORMATTED", status=&
                "OLD", position="asis")
        else
            open(newunit=original_mixing_unit, file=trim(name(1))//".bm", form="UNFORMATTED", status=&
                "OLD", position="asis")
        end if
        read (original_mixing_unit, iostat=io_status) mixing_file_signature
        if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
            write (*, *) "Not a GRASP mixing file"
            stop
        end if
        write (copied_mixing_unit) mixing_file_signature
        read (original_mixing_unit) nelec, total_csf_count, nw, nvectot, vector_size, nblock
        write (copied_mixing_unit) nelec, total_csf_count, nw, nvectot, vector_size, nblock

        do iblk = 1, nblock
            read (original_mixing_unit) ib, ncf, nvec, iatjp, iaspa
            write (copied_mixing_unit) ib, ncf, nvec, iatjp, iaspa
            call alloc (eval, nvec, "EVAL", "prepare_transition_mixing_coefficients")
            call alloc (evec, ncf*nvec, "EVEC", "prepare_transition_mixing_coefficients")
            call alloc (ivec, nvec, "IVEC", "prepare_transition_mixing_coefficients")
            read (original_mixing_unit) (ivec(i),i=1,nvec)
            read (original_mixing_unit) hamiltonian_energy_shift, (eval(i),i=1,nvec)
            read (original_mixing_unit) ((evec(i + (j - 1)*ncf),i=1,ncf),j=1,nvec)

            write (copied_mixing_unit) (ivec(i),i=1,nvec)
            write (copied_mixing_unit) hamiltonian_energy_shift, (eval(i),i=1,nvec)
            write (copied_mixing_unit) ((evec(i + (j - 1)*ncf),i=1,ncf),j=1,nvec)
            call dalloc (eval, "EVAL", "prepare_transition_mixing_coefficients")
            call dalloc (evec, "EVEC", "prepare_transition_mixing_coefficients")
            call dalloc (ivec, "IVEC", "prepare_transition_mixing_coefficients")
        end do
        name(2) = name(1)
        close(original_mixing_unit)
        close(copied_mixing_unit)
        return
    end subroutine prepare_transition_mixing_coefficients
end module prepare_transition_mixing_coefficients_interface
