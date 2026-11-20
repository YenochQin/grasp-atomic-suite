! Ordinary rwfnestimate input selection and the rCSFs header.toml adapter.
module rwfnestimate_input
    ! Module dependencies
    use iso_fortran_env, only: error_unit
    use iso_c_binding, only: c_char, c_int32_t, c_null_char, c_size_t
    use csf_descriptor_c_bindings, only: rust_descriptor_read_radial_header
    use grasp_parameters, only: max_orbitals
    use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
        & orbital_symmetry_labels
    use atomic_calculation_state, only: electron_count
    implicit none (type, external)

    ! Module interface
    private
    public :: radial_input_options, parse_radial_options, load_radial_header

    ! Derived types
    type :: radial_input_options
        ! Components
        character(len=:), allocatable :: filename
    end type radial_input_options

contains

    subroutine require_input(condition, message)
        ! Arguments
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            write(error_unit, "(A)") "RWFNESTIMATE input: "//message
            error stop 1
        end if
    end subroutine require_input

    subroutine read_argument(position, value)
        ! Arguments
        integer, intent(in) :: position
        character(len=:), allocatable, intent(out) :: value

        ! Local declarations
        integer :: length, status
        call get_command_argument(position, length=length, status=status)
        call require_input(status == 0 .and. length > 0, "missing command argument; use --help")
        allocate(character(len=length) :: value)
        call get_command_argument(position, value=value, status=status)
        call require_input(status == 0, "cannot read command argument")
    end subroutine read_argument

    subroutine parse_radial_options(options)
        ! Arguments
        type(radial_input_options), intent(out) :: options

        ! Local declarations
        character(len=:), allocatable :: argument, basename
        logical :: selected
        integer :: position
        options%filename = ""
        selected = .false.
        position = 1
        do while (position <= command_argument_count())
            call read_argument(position, argument)
            select case (argument)
            case ("--header")
                call require_input(.not. selected, "duplicate --header")
                call read_argument(position+1, options%filename)
                selected = .true.
            case ("--help")
                print "(A)", "rwfnestimate [--header FILE.toml]"
                print "(A)", "Enter a state to read NAME_header.toml and write NAME.initial.wfn; --header overrides the input."
                print "(A)", "isodata and the interactive initial-orbital choices are required."
                stop
            case default
                call require_input(.false., "unknown option "//argument//"; use --help")
            end select
            position = position + 2
        end do
        basename = options%filename(scan(options%filename, "/\", back=.true.)+1:)
        select case (basename)
        case ("rwfn.inp", "erwf.sum", "erwf.dbg", "tmp_734")
            call require_input(.false., "input filename uses a reserved output name")
        end select
    end subroutine parse_radial_options

    subroutine load_radial_header(options, core_count)
        ! Arguments
        type(radial_input_options), intent(in) :: options
        integer, intent(out) :: core_count

        ! Local constants
        character(len=*), parameter :: letters = "spdfghiklm"

        ! Local declarations
        character(kind=c_char,len=1024) :: message
        integer(c_int32_t) :: status, info(3), orbitals(2,max_orbitals)
        integer :: count, total_electrons, i, l
        integer :: kappa(max_orbitals), momentum(max_orbitals)

        status = rust_descriptor_read_radial_header(options%filename//c_null_char, orbitals, &
            int(max_orbitals,c_size_t), info, message, len(message,kind=c_size_t))
        if (status /= 0) call require_input(.false., message(:index(message,c_null_char)-1))

        count = int(info(1))
        core_count = int(info(2))
        total_electrons = int(info(3))
        kappa(:count) = int(orbitals(2,:count))
        momentum(:count) = merge(kappa(:count), -kappa(:count)-1, kappa(:count) > 0)
        call require_input(all(momentum(:count) < len(letters)), "unsupported orbital letter for radial estimation")
        call require_input(total_electrons >= sum(2*abs(kappa(:core_count))), &
            "electron_count is smaller than the closed-core population")
        call require_input(total_electrons <= sum(2*abs(kappa(:count))), &
            "electron_count exceeds the orbital table capacity")

        orbital_count = count
        electron_count = total_electrons
        orbital_principal_n(:count) = int(orbitals(1,:count))
        orbital_kappa(:count) = kappa(:count)
        orbital_l(:count) = momentum(:count)
        orbital_two_j(:count) = 2*abs(kappa(:count))-1
        do i = 1, count
            l = momentum(i)
            orbital_symmetry_labels(i) = letters(l+1:l+1)//merge("-", " ", kappa(i) > 0)
        end do
        write(*, "(A,I0,A,I0,A,I0)") "Loaded TOML header: core subshells=", core_count, &
            ", orbitals=", orbital_count, ", total electrons=", electron_count
    end subroutine load_radial_header
end module rwfnestimate_input
