! Ordinary rwfnestimate input selection and the rCSFs header.toml adapter.
module rwfnestimate_input
   use iso_fortran_env, only: int64, error_unit
   use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
   use parameter_def, only: NNNW
   use orb_C, only: NW, NP, NAK, NKL, NKJ, NH
   use def_C, only: NELEC
   use setcsh_I, only: setcsh
   implicit none
   private
   public :: radial_input_options, parse_radial_options, load_radial_header

   type :: radial_input_options
      character(len=:), allocatable :: filename
      logical :: use_header = .false.
   end type radial_input_options

contains

   subroutine require_input(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(error_unit, '(A)') 'RWFNESTIMATE input: '//message
         error stop 1
      end if
   end subroutine require_input

   subroutine read_argument(position, value)
      integer, intent(in) :: position
      character(len=:), allocatable, intent(out) :: value
      integer :: length, status
      call get_command_argument(position, length=length, status=status)
      call require_input(status == 0 .and. length > 0, 'missing command argument; use --help')
      allocate(character(len=length) :: value)
      call get_command_argument(position, value=value, status=status)
      call require_input(status == 0, 'cannot read command argument')
   end subroutine read_argument

   subroutine parse_radial_options(options)
      type(radial_input_options), intent(out) :: options
      character(len=:), allocatable :: argument, basename
      logical :: selected
      integer :: position
      options%filename = 'rcsf.inp'
      selected = .false.
      position = 1
      do while (position <= command_argument_count())
         call read_argument(position, argument)
         select case (argument)
         case ('--header', '--csf-text')
            call require_input(.not. selected, 'select only one of --header and --csf-text')
            options%use_header = argument == '--header'
            call read_argument(position+1, options%filename)
            selected = .true.
         case ('--help')
            print '(A)', 'rwfnestimate [--header FILE.toml | --csf-text FILE]'
            print '(A)', 'Default: rcsf.inp. Header mode requires header_info.electron_count including core electrons.'
            print '(A)', 'isodata and the interactive initial-orbital choices are required in either mode.'
            stop
         case default
            call require_input(.false., 'unknown option '//argument//'; use --help')
         end select
         position = position + 2
      end do
      basename = options%filename(scan(options%filename, '/\', back=.true.)+1:)
      select case (basename)
      case ('rwfn.inp', 'erwf.sum', 'erwf.dbg', 'tmp_734')
         call require_input(.false., 'input filename uses a reserved output name')
      end select
   end subroutine parse_radial_options

   subroutine load_radial_header(options, core_count)
      type(radial_input_options), intent(in) :: options
      integer, intent(out) :: core_count
      if (options%use_header) then
         call load_toml_header(options%filename, core_count)
      else
         call setcsh(21, options%filename, core_count)
      end if
   end subroutine load_radial_header

   subroutine load_toml_header(filename, core_count)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: core_count
      type(toml_table), allocatable :: document
      type(toml_table), pointer :: info
      type(toml_array), pointer :: lines
      type(toml_error), allocatable :: error
      character(len=:), allocatable :: line
      integer(int64) :: total_electrons
      integer :: status, index, core_electrons, capacity

      call toml_load(document, filename, error=error)
      if (allocated(error)) then
         write(error_unit, '(A)') error%message
         call require_input(.false., 'cannot parse TOML header '//filename)
      end if
      call get_value(document, 'header_info', info, requested=.false.)
      call require_input(associated(info), 'missing or invalid [header_info] table')
      call get_value(info, 'electron_count', total_electrons, stat=status)
      call require_input(status == 0, &
         'header_info.electron_count must be an integer including core electrons; regenerate old headers with rCSFs')
      call require_input(total_electrons > 0 .and. total_electrons <= int(huge(NELEC), int64), &
                         'electron_count must be positive and fit a suite integer')
      call get_value(info, 'header_lines', lines, requested=.false.)
      call require_input(associated(lines), 'header_info.header_lines must be an array')
      call require_input(len(lines) == 5, 'header_info.header_lines must contain five strings')
      NW = 0
      core_count = 0
      do index = 1, 5
         call get_value(lines, index, line, stat=status)
         call require_input(status == 0, 'header_info.header_lines must contain strings only')
         call require_input(scan(line, achar(10)//achar(13)//achar(0)) == 0, &
                            'header_lines entries must not contain newline, carriage return or NUL')
         select case (index)
         case (1)
            call require_input(trim(line) == 'Core subshells:', 'invalid core header label')
         case (2)
            call append_subshells(line)
            core_count = NW
         case (3)
            call require_input(trim(line) == 'Peel subshells:', 'invalid peel header label')
         case (4)
            call append_subshells(line)
         case (5)
            call require_input(trim(line) == 'CSF(s):', 'invalid CSF header label')
         end select
      end do
      call require_input(NW > 0, 'empty orbital table')
      core_electrons = sum(NKJ(:core_count) + 1)
      capacity = sum(NKJ(:NW) + 1)
      call require_input(total_electrons >= core_electrons, 'electron_count is smaller than the closed-core population')
      call require_input(total_electrons <= capacity, 'electron_count exceeds the orbital table capacity')
      NELEC = int(total_electrons)
      write(*, '(A,I0,A,I0,A,I0)') 'Loaded TOML header: core subshells=', core_count, &
         ', orbitals=', NW, ', total electrons=', NELEC
   end subroutine load_toml_header

   subroutine append_subshells(text)
      character(len=*), intent(in) :: text
      integer :: first, last, length
      length = len_trim(text)
      first = 1
      do while (first <= length)
         if (scan(text(first:first), ' ,'//achar(9)) > 0) then
            first = first + 1
            cycle
         end if
         last = first
         do while (last <= length)
            if (scan(text(last:last), ' ,'//achar(9)) > 0) exit
            last = last + 1
         end do
         call append_subshell(text(first:last-1))
         first = last + 1
      end do
   end subroutine append_subshells

   subroutine append_subshell(label)
      character(len=*), intent(in) :: label
      character(len=*), parameter :: letters = 'spdfghiklm'
      integer :: digits, principal, momentum, kappa, status, previous
      character(len=:), allocatable :: suffix
      call require_input(NW < NNNW, 'orbital table exceeds configured NNNW capacity')
      digits = 0
      do while (digits < len(label))
         if (scan(label(digits+1:digits+1), '0123456789') == 0) exit
         digits = digits + 1
      end do
      call require_input(digits > 0 .and. digits < len(label), 'invalid subshell label '//label)
      read(label(:digits), *, iostat=status) principal
      call require_input(status == 0, 'invalid principal quantum number in '//label)
      suffix = label(digits+1:)
      call require_input(len(suffix) <= 2, 'invalid subshell suffix in '//label)
      momentum = index(letters, suffix(1:1)) - 1
      call require_input(momentum >= 0, 'unsupported orbital letter in '//label)
      call require_input(principal > momentum, 'principal quantum number must exceed l in '//label)
      if (len(suffix) == 2) then
         call require_input(suffix(2:2) == '-' .and. momentum > 0, 'invalid relativistic suffix in '//label)
         kappa = momentum
      else
         kappa = -momentum - 1
      end if
      do previous = 1, NW
         call require_input(NP(previous) /= principal .or. NAK(previous) /= kappa, &
                            'duplicate orbital across core/peel lists: '//label)
      end do
      NW = NW + 1
      NP(NW) = principal
      NAK(NW) = kappa
      NKL(NW) = momentum
      NKJ(NW) = 2*abs(kappa) - 1
      NH(NW) = suffix
   end subroutine append_subshell
end module rwfnestimate_input
