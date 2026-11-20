! Shared state-name input and filenames for the CSF -> initial orbitals -> SCF chain.
module calculation_state_files
    ! Module dependencies
    use iso_fortran_env, only: input_unit, error_unit, iostat_eor, iostat_end
    use yes_no_input, only: getyn
    implicit none (type, external)

    ! Module interface
    private
    public :: state_files, state_name_length, read_state_name, select_state, files_for_state
    public :: read_state_entry, read_state_pair, property_state_arguments, resolve_mixing_file
    public :: named_angular_output
    public :: orbital_input_path, mixing_input_path, ci_mixing_suffix
    public :: check_orbital_output

    ! Module constants
    ! RCI historically accepts up to 127 characters before its trailing blank.
    integer, parameter :: state_name_length = 127
    character(*), parameter :: orbital_suffix = ".wfn", mixing_suffix = ".mix", ci_mixing_suffix = ".cimix"
    character(*), parameter :: initial_orbital_suffix = ".initial.wfn"

    ! Derived types
    type :: state_files
        ! Components
        character(:), allocatable :: descriptors, header, initial_orbitals, orbitals, mixing
        character(:), allocatable :: mcp, temporary, angular_log, angular_debug
        character(:), allocatable :: estimate_summary, estimate_debug, scf_log, scf_summary, scf_debug
    end type
contains
    function files_for_state(name) result(files)
        ! Arguments
        character(*), intent(in) :: name
        type(state_files) :: files
        files%descriptors      = filename("rcsf_descriptors.parquet", "_descriptors.parquet")
        files%header           = filename("rcsf_header.toml", "_header.toml")
        files%initial_orbitals = filename("rwfn.inp", initial_orbital_suffix)
        files%orbitals         = filename("rwfn.out", orbital_suffix)
        files%mixing           = filename("rmix.out", mixing_suffix)
        files%mcp              = "mcp"
        files%temporary        = filename("tmp", ".tmp")
        files%angular_log      = filename("rangular.log", ".rangular.log")
        files%angular_debug    = filename("genmcp.dbg", ".rangular.dbg")
        files%estimate_summary = filename("erwf.sum", ".rwfnestimate.summary.txt")
        files%estimate_debug   = filename("erwf.dbg", ".rwfnestimate.dbg")
        files%scf_log          = filename("rmcdhf.log", ".rmcdhf.log")
        files%scf_summary      = filename("rmcdhf.sum", ".rmcdhf.summary")
        files%scf_debug        = filename("rscf92.dbg", ".rmcdhf.dbg")
    contains
        function filename(legacy, suffix) result(path)
            ! Arguments
            character(*), intent(in) :: legacy, suffix
            character(:), allocatable :: path
            if (name == "") then
                path = legacy
            else
                path = trim(name)//suffix
            end if
        end function
    end function

    ! Accept an explicit mixing file or resolve a state to one of its existing mixing files.
    function resolve_mixing_file(entry, interactive) result(path)
        ! Arguments
        character(*), intent(in) :: entry
        logical, intent(in) :: interactive
        character(:), allocatable :: path, scf_path, ci_path

        ! Local declarations
        logical :: exists, has_scf, has_ci
        path = entry
        inquire(file=path, exist=exists)
        if (exists) return
        scf_path = mixing_input_path(entry, .false.)
        ci_path = mixing_input_path(entry, .true.)
        inquire(file=scf_path, exist=has_scf)
        inquire(file=ci_path, exist=has_ci)
        if (has_ci .and. has_scf) then
            if (.not. interactive) error stop "Both mixing files exist: specify "//scf_path//" or "//ci_path
            write(error_unit,"(A)") "Mixing coefficients from a CI calc.? (y/n)"
            has_ci = getyn()
        end if
        if (has_scf) path = scf_path
        if (has_ci) path = ci_path
    end function

    ! A.initial.wfn can be the initial orbitals of A or the final orbitals of A.initial.
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
        inquire(file=peer//"_descriptors.parquet", exist=peer_descriptors)
        if (peer_header .or. peer_descriptors) then
            write(error_unit,"(A)") "Orbital path also belongs to state "//peer//": "//trim(path)
            error stop "Use separate calculation directories for these states"
        end if
    end subroutine

    function orbital_input_path(name) result(path)
        ! Arguments
        character(*), intent(in) :: name
        character(:), allocatable :: path
        path = existing_input(name, orbital_suffix, ".w")
    end function

    function mixing_input_path(name, from_ci) result(path)
        ! Arguments
        character(*), intent(in) :: name
        logical, intent(in) :: from_ci
        character(:), allocatable :: path
        if (from_ci) then
            path = existing_input(name, ci_mixing_suffix, ".cm")
        else
            path = existing_input(name, mixing_suffix, ".m")
        end if
    end function

    ! Prefer the canonical input, including when reporting a missing file.
    function existing_input(name, suffix, old_suffix) result(path)
        ! Arguments
        character(*), intent(in) :: name, suffix, old_suffix
        character(:), allocatable :: path

        ! Local declarations
        logical :: exists
        path = trim(name)//suffix
        inquire(file=path, exist=exists)
        if (exists) return
        inquire(file=trim(name)//old_suffix, exist=exists)
        if (exists) path = trim(name)//old_suffix
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
        case ("log", "dbg")
            named_angular_output = channel == "rangular"
        case default
            if (len(basename)-dot /= 2 .or. verify(basename(dot+1:), "0123456789") /= 0) return
            if (len(channel) == 6) then
                if (verify(channel(4:), "0123456789") == 0) channel = channel(:3)
            end if
            named_angular_output = channel == "tmp"
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

    ! Several consumers locate the end of a name by finding its first blank.
    subroutine assign_state_name(text, name, status)
        ! Arguments
        character(*), intent(in) :: text
        character(*), intent(out) :: name
        integer, intent(out) :: status

        ! Local declarations
        integer :: limit
        limit = min(len(name)-1, state_name_length)
        name = ""
        if (len_trim(text) == 0 .or. len_trim(text) > limit .or. scan(trim(text), " "//achar(9)) /= 0) then
            write(error_unit,"(A,I0,A)") "Enter a state name of 1 to ", limit, " characters without spaces."
            status = 2
            return
        end if
        name = text
        status = 0
    end subroutine

    subroutine read_state_name(name, legacy_answer, directory, role, status)
        ! Arguments
        character(*), intent(out) :: name
        character, optional, intent(out) :: legacy_answer
        character(*), optional, intent(in) :: directory, role
        integer, optional, intent(out) :: status

        ! Local declarations
        character(:), allocatable :: record, token, parent
        character(256) :: io_message
        integer :: ios, comment
        logical :: named_header
        name = ""
        if (present(status)) status = 0
        if (present(legacy_answer)) legacy_answer = " "
        do
            call read_state_entry(record, ios, role, io_message)
            if (ios /= 0 .and. .not. (ios == iostat_end .and. len(record) > 0)) then
                write(error_unit,"(A)") "Cannot read state name: "//trim(io_message)
                if (.not. present(status)) error stop 1
                status = 2
                return
            end if
            if (present(legacy_answer)) then
                if (record == "") return
                ! Existing answer files begin with the first y/n response, often with a log comment.
                comment = index(record, "!")
                token = trim(record)
                if (comment > 0) token = trim(record(:comment-1))
                select case (token)
                case ("y","Y","n","N")
                    parent = "."
                    if (present(directory)) parent = trim(directory)
                    inquire(file=parent//"/"//token//"_header.toml", exist=named_header)
                    if (.not. named_header .or. comment > 0) then
                        legacy_answer = token
                        return
                    end if
                end select
                if (scan(record, "/"//achar(92)) /= 0) then
                    write(error_unit,"(A)") "Enter a state name without a directory."
                    cycle
                end if
            end if
            call assign_state_name(record, name, ios)
            if (ios == 0) return
        end do
    end subroutine

    subroutine read_state_pair(names)
        ! Arguments
        character(*), intent(out) :: names(2)
        call read_state_name(names(1), role="Initial state")
        call read_state_name(names(2), role="Final state")
    end subroutine

    subroutine select_state(name, question, answer, directory)
        ! Arguments
        character(*), intent(out) :: name
        character(*), intent(in) :: question
        logical, intent(out) :: answer
        character(*), optional, intent(in) :: directory

        ! Local declarations
        character :: legacy_answer
        character(len(name)+1) :: padded_name
        call read_state_name(padded_name, legacy_answer, directory)
        name = trim(padded_name)
        if (legacy_answer == " ") then
            write(error_unit,"(A)") question
            answer = getyn()
        else
            answer = index("Yy", legacy_answer) > 0
        end if
    end subroutine

    ! Positional state and CI selection shared by property programs with a CLI.
    ! Only density calculations expose the optional ordering argument.
    subroutine property_state_arguments(program, name, nci, status, ordering)
        ! Arguments
        character(*), intent(in) :: program
        character(*), intent(out) :: name
        integer, intent(inout) :: nci
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
                nci = 0
            case ("--nonci")
                nci = 1
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
