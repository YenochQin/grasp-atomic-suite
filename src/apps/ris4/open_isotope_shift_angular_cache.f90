! Open or validate cached angular coefficients for the isotope-shift calculation.
module open_isotope_shift_angular_cache_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_isotope_shift_angular_cache
contains

    !    Checks if the angular files name.IOB (one-body) or
    !    name.ITB (Two-body) are available
    !
    ! Modified by C. Naz\'e  Oct. 2011
    !
    ! Translated by Gediminas Gaigalas 11/18/19
    subroutine open_isotope_shift_angular_cache(name,avail,whichone,status,cache_unit)
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: iostat_end
        implicit none (type, external)

        ! Arguments
        character (len = 24), intent(in) :: name
        logical, intent(out) :: avail
        integer, intent(in)  :: whichone
        integer, intent(out) :: status
        integer, intent(out) :: cache_unit

        ! Local declarations
        logical :: found
        integer :: cache_probe_unit, istat, iword, close_status
        character(len=28) :: filename
        character :: probe

        avail = .false.
        cache_unit = -1
        status = 0
        if (whichone==1) then
            filename = trim(name)//".IOB"
        else if (whichone==2) then
            filename = trim(name)//".ITB"
        else
            status = 1
            return
        end if
        inquire(file=trim(filename), exist=found, iostat=status)
        if (status /= 0) return
        if (.not.found) then
            if (whichone==1) print *, " One-body angular file not available"
            if (whichone==2) print *, " Two-body angular file not available"
            avail = .false.
            return
        end if
        ! Missing/incomplete optional caches may be recomputed. Access and
        ! positioning failures on an otherwise usable cache are fatal to callers.
        ! Check raw readability separately from sequential-record decoding so an
        ! OS read error is not mistaken for an incomplete record (e.g. a directory).
        open(newunit=cache_probe_unit,file=trim(filename),status="OLD",action="READ", &
            access="STREAM",form="UNFORMATTED",iostat=status)
        if (status /= 0) return
        read(cache_probe_unit,iostat=status) probe
        if (status == iostat_end) status = 0
        close(unit=cache_probe_unit,iostat=close_status)
        if (close_status /= 0 .and. status == 0) status = close_status
        if (status /= 0) return
        open(newunit=cache_unit,file=trim(filename),status="OLD",action="READ", &
            position="APPEND",form="UNFORMATTED",iostat=status)
        if (status /= 0) return
        backspace(unit=cache_unit,iostat=istat)
        if (istat == 0) then
            read(cache_unit,iostat=istat) iword
            if (istat == 0) avail = iword == -1
        end if
        if (avail) then
            rewind(unit=cache_unit,iostat=status)
            if (status == 0) then
                if (whichone==1) print *, " One-body angular file available"
                if (whichone==2) print *, " Two-body angular file available"
                return
            end if
        end if
        avail = .false.
        close(unit=cache_unit,iostat=close_status)
        if (close_status /= 0 .and. status == 0) status = close_status
        if (status /= 0) return
        cache_unit = -1
        print *, " Angular file incomplete: ", trim(filename), "; iostat=", istat
        return
    end subroutine open_isotope_shift_angular_cache
end module open_isotope_shift_angular_cache_interface
