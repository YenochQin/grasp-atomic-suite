! File handoff between the biorthogonal transform and radiative transitions.
module biorthogonal_state_files
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: iostat_end, output_unit
    use calculation_state_files, only: state_files, mixing_input_path, orbital_input_path
    use input_error_handling, only: abort_input, fail_input_io
    use posix_services, only: make_directory
    implicit none (type, external)
    private

    public :: biorthogonal_pair, files_for_biorthogonal_pair, begin_biorthogonal_output
    public :: complete_biorthogonal_output, validate_biorthogonal_input, release_biorthogonal_lock

    character(*), parameter :: manifest_signature = "GRASP-BIORTH-1"

    type :: biorthogonal_pair
        character(:), allocatable :: directory
        type(state_files) :: sides(2), originals(2)
        character(5), private :: source = ""
        character(64), private :: source_digests(4) = ""
        integer(c_int), private :: lock = -1
    end type

    interface
        integer(c_int) function file_sha256(path, digest) bind(c, name="grasp_file_sha256")
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            character(kind=c_char), intent(out) :: digest(64)
        end function
        integer(c_int) function lock_state_files(path, writing) bind(c, name="grasp_lock_state_files")
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: writing
        end function
        subroutine unlock_state_files(lock) bind(c, name="grasp_unlock_state_files")
            import :: c_int
            integer(c_int), value :: lock
        end subroutine
    end interface
contains
    function files_for_biorthogonal_pair(names, from_ci, startup_directory) result(pair)
        character(*), intent(in) :: names(2), startup_directory
        logical, intent(in) :: from_ci
        type(biorthogonal_pair) :: pair

        character(:), allocatable :: root, source_prefix, destination
        character(7), parameter :: sides(2) = [character(7) :: "initial", "final"]
        integer :: side

        root = trim(startup_directory)//"/"
        pair%source = merge("ci   ", "mcdhf", from_ci)
        do side = 1, 2
            source_prefix = root//trim(names(side))
            pair%originals(side)%orbitals = orbital_input_path(source_prefix)
            pair%originals(side)%mixing = mixing_input_path(source_prefix, from_ci)
        end do
        if (names(1) == names(2)) then
            pair%sides = pair%originals
            return
        end if

        pair%directory = root//"biorth/"//trim(names(1))//"/"//trim(names(2))//"/"//trim(pair%source)
        do side = 1, 2
            destination = pair%directory//"/"//trim(sides(side))//".biorth"
            pair%sides(side)%orbitals = orbital_input_path(destination)
            pair%sides(side)%mixing = mixing_input_path(destination, from_ci)
        end do
    end function

    subroutine begin_biorthogonal_output(pair, abort_handler)
        type(biorthogonal_pair), intent(inout) :: pair
        procedure(abort_input), optional :: abort_handler

        character(:), allocatable :: manifest
        integer :: unit, status
        logical :: exists

        call make_directory(pair%directory, status)
        if (status /= 0) call fail_input_io(status, -1, "Cannot create "//pair%directory, abort_handler)
        pair%lock = lock_state_files(pair%directory//"/.lock"//c_null_char, 1_c_int)
        if (pair%lock < 0) call fail_input_io(0, -1, "Cannot lock "//pair%directory, abort_handler)
        manifest = pair%directory//"/manifest"
        inquire(file=manifest, exist=exists)
        if (exists) then
            open(newunit=unit, file=manifest, status="old", iostat=status)
            if (status /= 0) call fail_input_io(status, -1, "Cannot open "//manifest, abort_handler)
            close(unit, status="delete", iostat=status)
            if (status /= 0) call fail_input_io(status, unit, "Cannot remove "//manifest, abort_handler)
        end if
        call pair_digests(pair%originals, pair%source_digests, abort_handler)
    end subroutine

    subroutine complete_biorthogonal_output(pair, names, basis_ids, abort_handler)
        type(biorthogonal_pair), intent(inout) :: pair
        character(*), intent(in) :: names(2), basis_ids(2)
        procedure(abort_input), optional :: abort_handler

        character(64) :: output_digests(4)
        character(:), allocatable :: manifest
        integer :: unit, status

        call pair_digests(pair%sides, output_digests, abort_handler)
        manifest = pair%directory//"/manifest"
        open(newunit=unit, file=manifest, status="new", action="write", iostat=status)
        if (status /= 0) call fail_input_io(status, -1, "Cannot create "//manifest, abort_handler)
        write(unit, "(a)", iostat=status) manifest_signature, trim(names(1)), trim(names(2)), trim(pair%source), basis_ids, &
            pair%source_digests, output_digests
        if (status /= 0) call fail_input_io(status, unit, "Cannot write "//manifest, abort_handler)
        close(unit, iostat=status)
        if (status /= 0) call fail_input_io(status, unit, "Cannot complete "//manifest, abort_handler)
        write(output_unit, "(a)") "Biorthogonal result directory: "//pair%directory
        call release_biorthogonal_lock(pair)
    end subroutine

    subroutine validate_biorthogonal_input(pair, names, basis_ids, abort_handler)
        type(biorthogonal_pair), intent(inout) :: pair
        character(*), intent(in) :: names(2), basis_ids(2)
        procedure(abort_input), optional :: abort_handler

        character(:), allocatable :: manifest, rerun
        character(64) :: signature, saved_names(2), saved_source, saved_basis_ids(2)
        character(64) :: source_digests(4), output_digests(4), actual(4)
        integer :: unit, status, end_status

        if (.not. allocated(pair%directory)) return
        manifest = pair%directory//"/manifest"
        rerun = ". Run biorthtransform or biorthtransform_mpi for this ordered pair and mixing source."
        pair%lock = lock_state_files(pair%directory//"/.lock"//c_null_char, 0_c_int)
        if (pair%lock < 0) call fail_input_io(0, -1, "Cannot lock biorthogonal results at "//pair%directory//rerun, abort_handler)
        open(newunit=unit, file=manifest, status="old", action="read", iostat=status)
        if (status /= 0) call fail_input_io(status, -1, "Missing biorthogonal manifest "//manifest//rerun, abort_handler)
        read(unit, "(a)", iostat=status) signature, saved_names, saved_source, saved_basis_ids, source_digests, output_digests
        if (status == 0) then
            read(unit, "(a)", iostat=end_status)
            if (end_status /= iostat_end) status = 1
        end if
        close(unit)
        if (status /= 0) &
            call fail_input_io(status, -1, "Invalid biorthogonal manifest "//manifest//rerun, abort_handler)
        if (signature /= manifest_signature .or. any(saved_names /= names) .or. saved_source /= pair%source .or. &
            any(saved_basis_ids /= basis_ids)) &
            call fail_input_io(0, -1, "Invalid biorthogonal manifest "//manifest//rerun, abort_handler)
        call pair_digests(pair%originals, actual, abort_handler)
        if (any(source_digests /= actual)) &
            call fail_input_io(0, -1, "Biorthogonal input digest mismatch at "//manifest//rerun, abort_handler)
        call pair_digests(pair%sides, actual, abort_handler)
        if (any(output_digests /= actual)) &
            call fail_input_io(0, -1, "Biorthogonal output digest mismatch at "//manifest//rerun, abort_handler)
    end subroutine

    subroutine pair_digests(paths, digests, abort_handler)
        type(state_files), intent(in) :: paths(2)
        character(64), intent(out) :: digests(4)
        procedure(abort_input), optional :: abort_handler

        integer :: side

        do side = 1, 2
            if (file_sha256(paths(side)%orbitals//c_null_char, digests(side)) /= 0) &
                call fail_input_io(0, -1, "Cannot hash biorthogonal orbital file "//paths(side)%orbitals, abort_handler)
            if (file_sha256(paths(side)%mixing//c_null_char, digests(2+side)) /= 0) &
                call fail_input_io(0, -1, "Cannot hash biorthogonal mixing file "//paths(side)%mixing, abort_handler)
        end do
    end subroutine

    subroutine release_biorthogonal_lock(pair)
        type(biorthogonal_pair), intent(inout) :: pair

        if (pair%lock < 0) return
        call unlock_state_files(pair%lock)
        pair%lock = -1
    end subroutine
end module biorthogonal_state_files
