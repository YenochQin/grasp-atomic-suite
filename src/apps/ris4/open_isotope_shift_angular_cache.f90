! Open or validate cached angular coefficients for the isotope-shift calculation.
module open_isotope_shift_angular_cache_interface
    use, intrinsic :: iso_fortran_env, only: int64
    use csf_descriptor_reader, only: loaded_csf_basis_id
    use debug_output_settings, only: matrix_element_screening_cutoff
    use grasp_parameters, only: orbital_index_radix
    implicit none (type, external)

    ! Module interface
    private
    public :: open_isotope_shift_angular_cache, finish_angular_cache, isotope_shift_cache_path
    character(8), parameter :: cache_version = "RISANG01"
contains
    function isotope_shift_cache_path(name, operator) result(path)
        character(*), intent(in) :: name
        integer, intent(in) :: operator
        character(:), allocatable :: path
        character(*), parameter :: suffixes(2) = [".ris4.one-body.angular.cache", ".ris4.two-body.angular.cache"]

        path = trim(name)//suffixes(operator)
    end function isotope_shift_cache_path

    subroutine finish_angular_cache(unit, operator)
        integer, intent(in) :: unit, operator

        ! Keep the short end record consumed by the coefficient readers.
        write(unit) -1
        write(unit) cache_version, loaded_csf_basis_id, operator, orbital_index_radix, &
            transfer(matrix_element_screening_cutoff, 0_int64)
    end subroutine finish_angular_cache

    ! Checks whether the one- or two-body angular cache is available.
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
        character(:), allocatable :: filename
        character(8) :: version
        character(64) :: basis
        integer :: operator, radix
        integer(int64) :: cutoff_bits
        character :: probe

        avail = .false.
        cache_unit = -1
        status = 0
        if (whichone < 1 .or. whichone > 2) then
            status = 1
            return
        end if
        filename = isotope_shift_cache_path(name, whichone)
        inquire(file=filename, exist=found, iostat=status)
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
            read(cache_unit,iostat=istat) version, basis, operator, radix, cutoff_bits
            if (istat == 0) avail = version == cache_version .and. basis == loaded_csf_basis_id .and. &
                operator == whichone .and. radix == orbital_index_radix .and. &
                cutoff_bits == transfer(matrix_element_screening_cutoff, 0_int64)
            if (avail) then
                backspace(cache_unit,iostat=istat)
                if (istat == 0) backspace(cache_unit,iostat=istat)
                if (istat == 0) read(cache_unit,iostat=istat) iword
                avail = .false.
                if (istat == 0) avail = iword == -1
            end if
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
        print *, " Angular file incomplete or incompatible: ", trim(filename), "; iostat=", istat
        return
    end subroutine open_isotope_shift_angular_cache
end module open_isotope_shift_angular_cache_interface
