module rangular_input
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                   install_descriptor_header, read_descriptor_block
   use rangular_profile, only: open_stage_profile
   use rangular_reuse, only: configure_angular_reuse
   implicit none (type, external)
   private
   public :: angular_input_options, parse_angular_options, load_angular_header, load_angular_block, close_angular_input
   type :: angular_input_options
      character(len=:), allocatable :: input_filename, timings_filename, header_filename
      type(descriptor_source) :: descriptors
   end type angular_input_options
contains
   subroutine reject_output_filename(filename)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: basename
      integer :: separator_position
      separator_position = scan(filename, "/\", back=.true.)
      basename = filename(separator_position+1:)
      if (basename == "rangular.log" .or. basename == "genmcp.dbg" .or. &
          index(basename, "mcp.") == 1 .or. index(basename, "tmp.") == 1) &
         error stop "Input/timings filename uses a reserved rangular output name"
   end subroutine reject_output_filename

   ! Read directly into the caller-owned string, with its complete argument length.
   subroutine read_argument_text(argument_index, text)
      integer, intent(in) :: argument_index
      character(len=:), allocatable, intent(out) :: text
      integer :: text_length, argument_status
      call get_command_argument(argument_index, length=text_length, status=argument_status)
      if (argument_status /= 0 .or. text_length == 0) error stop "Missing command argument"
      allocate(character(len=text_length) :: text)
      call get_command_argument(argument_index, value=text, status=argument_status)
      if (argument_status /= 0) error stop "Cannot read command argument"
   end subroutine read_argument_text

   subroutine parse_angular_options(options)
      type(angular_input_options), intent(out) :: options
      integer :: argument_index
      integer :: integer_status, option_value
      integer :: reuse_memory_mib
      logical :: input_selected, timings_selected, reuse_selected, reuse_stats_selected
      character(len=:), allocatable :: option_name, reuse_filename
      options%input_filename = "rcsf_descriptors.parquet"
      options%timings_filename = ""
      options%header_filename = ""
      input_selected = .false.
      timings_selected = .false.
      reuse_selected = .false.
      reuse_stats_selected = .false.
      reuse_memory_mib = 0
      reuse_filename = ""
      argument_index = 1
      do while (argument_index <= command_argument_count())
         call read_argument_text(argument_index, option_name)
         select case (option_name)
         case ("--csf-descriptors")
            if (input_selected) error stop "Select exactly one CSF input option"
            call read_argument_text(argument_index+1, options%input_filename)
            input_selected = .true.
         case ("--header")
            if (len(options%header_filename) > 0) error stop "Duplicate --header option"
            call read_argument_text(argument_index+1,options%header_filename)
         case ("--timings")
            if (timings_selected) error stop "Duplicate --timings option"
            call read_argument_text(argument_index+1, options%timings_filename)
            timings_selected = .true.
         case ("--reuse-stats")
            if (reuse_stats_selected) error stop "Duplicate --reuse-stats option"
            call read_argument_text(argument_index+1, reuse_filename)
            reuse_stats_selected = .true.
         case ("--reuse-cache-mib")
            block
               character(len=:), allocatable :: argument_value
               call read_argument_text(argument_index+1, argument_value)
               read(argument_value, *, iostat=integer_status) option_value
            end block
            if (integer_status /= 0) error stop "--reuse-cache-mib requires an integer"
            if (reuse_selected) error stop "Duplicate --reuse-cache-mib option"
            if (option_value < 0 .or. option_value > 4096) error stop "--reuse-cache-mib must be in 0..4096"
            reuse_memory_mib = option_value
            reuse_selected = .true.
         case ("--help")
            print "(A)", "rangular [--csf-descriptors FILE] [--header FILE.toml] [--timings NEW_FILE.csv]"
            print "(A)", "Default: rcsf_descriptors.parquet with matching rcsf_header.toml."
            print "(A)", "Serial only: --reuse-cache-mib N [--reuse-stats NEW_FILE.csv]; default 0 (off)."
            stop
         case default
            error stop "Unknown rangular option; use --help"
         end select
         argument_index = argument_index + 2
      end do
      call reject_output_filename(options%header_filename)
      call reject_output_filename(options%input_filename)
      call reject_output_filename(options%timings_filename)
      call reject_output_filename(reuse_filename)
      if (len(reuse_filename) > 0) then
         if (reuse_filename == options%input_filename .or. reuse_filename == options%timings_filename) &
            error stop "Reuse statistics must use a separate filename"
      end if
      call configure_angular_reuse(reuse_memory_mib,reuse_filename)
      call open_stage_profile(options%timings_filename)
   end subroutine parse_angular_options

   subroutine load_angular_header(options, core_count, maximum_blocks, block_labels)
      type(angular_input_options), intent(inout) :: options
      integer, intent(out) :: core_count
      integer, intent(in) :: maximum_blocks
      character(len=8), intent(out) :: block_labels(:)
      call open_descriptor_source(options%descriptors, options%input_filename, maximum_blocks,options%header_filename)
      call install_descriptor_header(options%descriptors, core_count, block_labels)
   end subroutine load_angular_header

   subroutine load_angular_block(options, core_count, block_index)
      use orb_c, only: occupations => csf_subshell_occupations
      use stat_c, only: subshell_quantum_numbers => csf_subshell_quantum_codes, couplings => csf_coupling_quantum_codes
      type(angular_input_options), intent(inout) :: options
      integer, intent(in) :: core_count, block_index
      call read_descriptor_block(options%descriptors, block_index, occupations, subshell_quantum_numbers, couplings)
   end subroutine load_angular_block

   subroutine close_angular_input(options)
      type(angular_input_options), intent(inout) :: options
      call close_descriptor_source(options%descriptors)
   end subroutine close_angular_input
end module rangular_input
