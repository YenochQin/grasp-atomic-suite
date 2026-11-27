! Shared state-name input and paths for Suite calculations.
module calculation_state_files
    ! Module dependencies
    use iso_fortran_env, only: input_unit, error_unit, iostat_eor, iostat_end
    use yes_no_input, only: read_yes_no
    implicit none (type, external)

    ! Module interface
    private
    public :: state_files, state_name_length, read_state_name, read_state_settings, files_for_state
    public :: read_state_entry, read_state_pair, property_state_arguments, resolve_mixing_file
    public :: named_angular_output
    public :: orbital_input_path, mixing_input_path, ci_mixing_suffix
    public :: natural_orbital_path, isotope_shift_report_path, isotope_shift_input_path
    public :: lsj_label_path, lsj_label_input_path, transition_report_path
    public :: descriptor_input_path
    public :: rci_restart_path
    public :: check_orbital_output

    ! Module constants
    ! RCI historically accepts up to 127 characters before its trailing blank.
    integer, parameter :: state_name_length = 127
    character(*), parameter :: orbital_suffix = ".wfn", mixing_suffix = ".mix", ci_mixing_suffix = ".cimix"
    character(*), parameter :: initial_orbital_suffix = ".initial.wfn"
    character(*), parameter :: descriptor_suffix = "_desc.parquet"

    ! Derived types
    type :: state_files
        ! Components
        character(:), allocatable :: descriptors, header, initial_orbitals, orbitals, mixing
        character(:), allocatable :: mcp, indexed_mcp_prefix, temporary, angular_log, angular_debug
        character(:), allocatable :: estimate_summary, estimate_debug, scf_log, scf_summary, scf_debug
    end type
contains
    function files_for_state(name) result(files)
        ! Arguments
        character(*), intent(in) :: name
        type(state_files) :: files
        files%descriptors      = descriptor_input_path(name)
        files%header           = trim(name)//"_header.toml"
        files%initial_orbitals = trim(name)//initial_orbital_suffix
        files%orbitals         = orbital_input_path(name)
        files%mixing           = mixing_input_path(name, .false.)
        files%mcp              = trim(name)//".mcp"
        files%indexed_mcp_prefix = trim(name)//".mpi"
        files%temporary        = trim(name)//".tmp"
        files%angular_log      = trim(name)//".rangular.log"
        files%angular_debug    = trim(name)//".rangular.dbg"
        files%estimate_summary = trim(name)//".rwfnestimate.summary.txt"
        files%estimate_debug   = trim(name)//".rwfnestimate.dbg"
        files%scf_log          = trim(name)//".rmcdhf.log"
        files%scf_summary      = trim(name)//".rmcdhf.summary"
        files%scf_debug        = trim(name)//".rmcdhf.dbg"
    end function

    function rci_restart_path(name, rank) result(path)
        character(*), intent(in) :: name
        integer, optional, intent(in) :: rank
        character(:), allocatable :: path
        character(16) :: rank_text

        path = trim(name)//".rci"
        if (present(rank)) then
            write(rank_text,"(I0.3)") rank
            path = path//".rank"//trim(rank_text)
        end if
        path = path//".restart"
    end function

    function natural_orbital_path(name, from_ci) result(path)
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path

        path = source_report_stem(name, from_ci)//".natural-orbitals"
    end function natural_orbital_path

    function isotope_shift_report_path(name, from_ci) result(path)
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path

        path = source_report_stem(name, from_ci)//".isotope-shift"
    end function isotope_shift_report_path

    function isotope_shift_input_path(name, from_ci) result(path)
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path

        path = isotope_shift_report_path(name, from_ci)
        call use_legacy_path_if_missing(path, trim(name)//trim(merge(".ci", ".i ", from_ci)))
    end function isotope_shift_input_path

    function lsj_label_path(name, from_ci) result(path)
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path

        path = source_report_stem(name, from_ci)//".lsj.labels"
    end function lsj_label_path

    function lsj_label_input_path(name, from_ci) result(path)
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path

        path = lsj_label_path(name, from_ci)
        call use_legacy_path_if_missing(path, trim(name)//".lsj.lbl")
    end function lsj_label_input_path

    function transition_report_path(initial_state, final_state, from_ci, with_lsj, with_phase) result(path)
        character(*), intent(in) :: initial_state, final_state
        logical, intent(in) :: from_ci, with_lsj, with_phase
        character(:), allocatable :: path

        path = source_report_stem(trim(initial_state)//"."//trim(final_state), from_ci)//".transition"
        if (with_phase) path = path//".phase"
        if (with_lsj) path = path//".lsj"
    end function transition_report_path

    function source_report_stem(name, from_ci) result(stem)
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: stem

        stem = trim(name)//trim(merge(".ci   ", ".mcdhf", from_ci))
    end function source_report_stem

    subroutine use_legacy_path_if_missing(path, legacy_path)
        character(:), allocatable, intent(inout) :: path
        character(*), intent(in) :: legacy_path
        character(:), allocatable :: previous_path
        logical :: exists

        inquire(file=path, exist=exists)
        if (exists) return

        previous_path = path//".txt"
        inquire(file=previous_path, exist=exists)
        if (exists) then
            path = previous_path
        else
            path = legacy_path
        end if
    end subroutine use_legacy_path_if_missing

    ! Accept an explicit mixing file or resolve a state to one of its existing mixing files.
    function resolve_mixing_file(entry, interactive) result(path)
        ! Arguments
        character(*), intent(in) :: entry
        logical, intent(in) :: interactive
        character(:), allocatable :: path, scf_path, ci_path

        ! Local declarations
        logical :: exists, has_scf, has_ci
        integer :: suffix, status, basename_start
        path = entry
        suffix = index(path, ".", back=.true.)
        if (suffix > 0) then
            select case (path(suffix:))
            case (".m", ".M", ".cm", ".cM", ".Cm", ".CM")
                error stop "Use a .mix or .cimix mixing file: "//path
            end select
        end if
        inquire(file=path, exist=exists)
        if (exists) return
        basename_start = scan(path, "/"//achar(92), back=.true.) + 1
        call validate_state_name(path(basename_start:), state_name_length, status)
        if (status /= 0) error stop "Invalid state name or missing mixing file: "//entry
        scf_path = mixing_input_path(entry, .false.)
        ci_path = mixing_input_path(entry, .true.)
        inquire(file=scf_path, exist=has_scf)
        inquire(file=ci_path, exist=has_ci)
        if (has_ci .and. has_scf) then
            if (.not. interactive) error stop "Both mixing files exist: specify "//scf_path//" or "//ci_path
            write(error_unit,"(A)") "Mixing coefficients from a CI calc.? (y/n)"
            has_ci = read_yes_no()
        end if
        if (has_scf) path = scf_path
        if (has_ci) path = ci_path
    end function

    ! An archived state A.initial may still own the path A.initial.wfn.
    subroutine check_orbital_output(path, writing_initial)
        ! Arguments
        character(*), intent(in) :: path
        logical, intent(in) :: writing_initial

        ! Local declarations
        character(:), allocatable :: peer
        integer :: stem_length
        logical :: exists, peer_header, peer_descriptors
        stem_length = len_trim(path) - len(initial_orbital_suffix)
        if (stem_length < 1) return
        if (path(stem_length+1:) /= initial_orbital_suffix) return
        inquire(file=path, exist=exists)
        if (.not. exists) return
        peer = path(:stem_length)
        if (writing_initial) peer = peer//".initial"
        inquire(file=peer//"_header.toml", exist=peer_header)
        inquire(file=descriptor_input_path(peer), exist=peer_descriptors)
        if (peer_header .or. peer_descriptors) then
            write(error_unit,"(A)") "Orbital path also belongs to state "//peer//": "//trim(path)
            error stop "Use separate calculation directories for these states"
        end if
    end subroutine

    function descriptor_input_path(name) result(path)
        character(*), intent(in) :: name
        character(:), allocatable :: path

        path = trim(name)//descriptor_suffix
    end function

    function orbital_input_path(name) result(path)
        ! Arguments
        character(*), intent(in) :: name
        character(:), allocatable :: path
        path = trim(name)//orbital_suffix
    end function

    function mixing_input_path(name, from_ci) result(path)
        ! Arguments
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path
        if (from_ci) then
            path = trim(name)//ci_mixing_suffix
        else
            path = trim(name)//mixing_suffix
        end if
    end function

    ! Explicit CLI inputs and profiling files must not alias the angular writer.
    logical function named_angular_output(basename)
        ! Arguments
        character(*), intent(in) :: basename

        ! Local declarations
        character(:), allocatable :: stem, channel
        integer :: dot
        named_angular_output = .false.
        dot = index(basename, ".", back=.true.)
        if (dot == 0) return
        stem = basename(:dot-1)
        channel = stem(index(stem, ".", back=.true.)+1:)
        select case (basename(dot+1:))
        case ("mcpidx")
            named_angular_output = .true.
        case ("log", "dbg")
            named_angular_output = channel == "rangular"
        case default
            if (channel == "mcpidx" .and. index(basename(dot+1:), "part-") == 1) then
                named_angular_output = .true.
                return
            end if
            if (len(basename)-dot /= 2 .or. verify(basename(dot+1:), "0123456789") /= 0) return
            if (len(channel) == 6) then
                if (verify(channel(4:), "0123456789") == 0) channel = channel(:3)
            end if
            named_angular_output = channel == "tmp" .or. channel == "mcp"
        end select
    end function

    ! Read a complete record before applying the destination's name limit. RLEVELS
    ! also accepts explicit, arbitrarily long paths and uses blank or EOF to finish.
    subroutine read_state_entry(record, status, role, message)
        ! Arguments
        character(:), allocatable, intent(out) :: record
        integer, intent(out) :: status
        character(*), optional, intent(in) :: role
        character(*), optional, intent(out) :: message

        ! Local declarations
        character(256) :: chunk, io_message
        integer :: count
        if (present(role)) write(error_unit,"(A)",advance="no") trim(role)//" - "
        write(error_unit,"(A)") "Name of state:"
        record = ""
        do
            read(input_unit,"(A)",advance="no",size=count,iostat=status,iomsg=io_message) chunk
            if (count > 0) record = record//chunk(:count)
            if (status /= 0) exit
        end do
        if (status == iostat_eor) status = 0
        if (present(message) .and. status /= 0) message = io_message
    end subroutine

    subroutine validate_state_name(name, limit, status)
        ! Arguments
        character(*), intent(in) :: name
        integer, intent(in) :: limit
        integer, intent(out) :: status

        ! Local declarations
        integer :: i, code
        status = 2
        if (len(name) == 0 .or. len(name) > limit) return
        if (scan(name, "./"//achar(92)) /= 0) return
        do i = 1, len(name)
            code = iachar(name(i:i))
            if (code <= 32 .or. code == 127) return
        end do
        status = 0
    end subroutine

    ! Reserve the trailing blank used by existing filename consumers.
    subroutine assign_state_name(text, name, status)
        ! Arguments
        character(*), intent(in) :: text
        character(*), intent(out) :: name
        integer, intent(out) :: status

        ! Local declarations
        integer :: limit
        limit = min(len(name)-1, state_name_length)
        call validate_state_name(text, limit, status)
        if (status == 0) then
            name = text
        else
            name = ""
            write(error_unit,"(A,I0,A)") "Invalid state name. Enter 1 to ", limit, &
                " characters without dots, whitespace, control characters, or directory separators."
        end if
    end subroutine

    subroutine read_state_name(name, role, status)
        ! Arguments
        character(*), intent(out) :: name
        character(*), optional, intent(in) :: role
        integer, optional, intent(out) :: status

        ! Local declarations
        character(:), allocatable :: record
        character(256) :: io_message
        integer :: ios
        call read_state_entry(record, ios, role, io_message)
        if (ios == 0 .or. (ios == iostat_end .and. len(record) > 0)) then
            call assign_state_name(record, name, ios)
        else
            name = ""
            write(error_unit,"(A)") "Cannot read state name: "//trim(io_message)
            ios = 2
        end if
        if (present(status)) then
            status = ios
        else if (ios /= 0) then
            error stop 1
        end if
    end subroutine

    subroutine read_state_pair(first, second, roles, status)
        ! Arguments
        character(*), intent(out) :: first, second
        character(*), optional, intent(in) :: roles(2)
        integer, optional, intent(out) :: status

        ! Local declarations
        character(*), parameter :: default_roles(2) = [character(13) :: "Initial state", "Final state"]
        character(len(first)) :: first_name
        character(len(second)) :: second_name
        integer :: ios

        if (present(roles)) then
            call read_state_name(first_name, roles(1), ios)
            if (ios == 0) call read_state_name(second_name, roles(2), ios)
        else
            call read_state_name(first_name, default_roles(1), ios)
            if (ios == 0) call read_state_name(second_name, default_roles(2), ios)
        end if
        if (ios == 0) then
            first = first_name
            second = second_name
        else
            first = ""
            second = ""
        end if
        if (present(status)) then
            status = ios
        else if (ios /= 0) then
            error stop 1
        end if
    end subroutine

    subroutine read_state_settings(question, answer)
        ! Arguments
        character(*), intent(in) :: question
        logical, intent(out) :: answer
        write(error_unit,"(A)") question
        answer = read_yes_no()
    end subroutine

    ! Positional state and CI selection shared by property programs with a CLI.
    ! Only density calculations expose the optional ordering argument.
    subroutine property_state_arguments(program, name, mixing_file_kind, status, ordering)
        ! Arguments
        character(*), intent(in) :: program
        character(*), intent(out) :: name
        integer, intent(inout) :: mixing_file_kind
        integer, intent(out) :: status
        integer, optional, intent(inout) :: ordering

        ! Local declarations
        character(:), allocatable :: arg, usage
        integer :: i, length, ios
        name = ""
        status = 0
        do i = 1, command_argument_count()
            call get_command_argument(i, length=length)
            allocate(character(length) :: arg)
            call get_command_argument(i, value=arg, status=ios)
            if (ios /= 0 .or. length == 0) then
                status = 2
                exit
            end if
            select case (arg)
            case ("-h", "--help")
                status = 1
            case ("--ci")
                mixing_file_kind = 0
            case ("--nonci")
                mixing_file_kind = 1
            case ("--order=1", "--order=2")
                if (present(ordering)) then
                    ordering = merge(1, 2, arg == "--order=1")
                else
                    status = 2
                end if
            case default
                if (i == 1 .and. index(arg,"--") /= 1) then
                    call assign_state_name(arg, name, status)
                else
                    write(error_unit,"(A)") "Unrecognized argument: "//arg
                    status = 2
                end if
            end select
            if (status /= 0) exit
            deallocate(arg)
        end do
        if (status == 0 .and. name == "") status = 2
        if (status /= 0) then
            usage = "Usage: "//program//" name [--ci|--nonci]"
            if (present(ordering)) usage = usage//" [--order=1|--order=2]"
            write(error_unit,"(A)") usage
        end if
    end subroutine
end module calculation_state_files
