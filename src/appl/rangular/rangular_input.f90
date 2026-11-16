module rangular_input
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, close_descriptor_source, &
                                   install_descriptor_header, read_descriptor_block
   use rangular_profile, only: open_stage_profile
   use rangular_parallel, only: configure_pair_engine, parallel_pairs_available
   use rangular_reuse, only: configure_angular_reuse
   use cslh_i, only: cslh
   use lodcsh2_i, only: lodcsh2
   implicit none
   private
   public :: angular_input_options, parse_angular_options, load_angular_header, load_angular_block, close_angular_input
   type :: angular_input_options
      character(len=:), allocatable :: input_filename, timings_filename
      logical :: use_descriptors = .false.
      type(descriptor_source) :: descriptors
   end type angular_input_options
contains
   subroutine reject_output_filename(filename)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: basename
      integer :: separator_position
      separator_position = scan(filename, '/\', back=.true.)
      basename = filename(separator_position+1:)
      if (basename == 'rangular.log' .or. basename == 'genmcp.dbg' .or. &
          index(basename, 'mcp.') == 1 .or. index(basename, 'tmp.') == 1) &
         error stop 'Input/timings filename uses a reserved rangular output name'
   end subroutine reject_output_filename

   ! Read directly into the caller-owned string, with its complete argument length.
   subroutine read_argument_text(argument_index, text)
      integer, intent(in) :: argument_index
      character(len=:), allocatable, intent(out) :: text
      integer :: text_length, argument_status
      call get_command_argument(argument_index, length=text_length, status=argument_status)
      if (argument_status /= 0 .or. text_length == 0) error stop 'Missing command argument'
      allocate(character(len=text_length) :: text)
      call get_command_argument(argument_index, value=text, status=argument_status)
      if (argument_status /= 0) error stop 'Cannot read command argument'
   end subroutine read_argument_text

   subroutine parse_angular_options(options)
      type(angular_input_options), intent(out) :: options
      integer :: argument_index
      integer :: requested_threads, pair_chunk_size, result_memory_mib, integer_status, option_value
      integer :: reuse_memory_mib
      logical :: input_selected, timings_selected, reuse_selected, reuse_stats_selected
      character(len=:), allocatable :: option_name, reuse_filename
      options%input_filename = 'rcsf.inp'
      options%timings_filename = ''
      input_selected = .false.
      timings_selected = .false.
      reuse_selected = .false.
      reuse_stats_selected = .false.
      reuse_memory_mib = 0
      reuse_filename = ''
      requested_threads = 0
      pair_chunk_size = 256
      result_memory_mib = 64
      argument_index = 1
      do while (argument_index <= command_argument_count())
         call read_argument_text(argument_index, option_name)
         select case (option_name)
         case ('--csf-text', '--csf-descriptors')
            if (input_selected) error stop 'Select exactly one CSF input option'
            call read_argument_text(argument_index+1, options%input_filename)
            options%use_descriptors = option_name == '--csf-descriptors'
            input_selected = .true.
         case ('--timings')
            if (timings_selected) error stop 'Duplicate --timings option'
            call read_argument_text(argument_index+1, options%timings_filename)
            timings_selected = .true.
         case ('--reuse-stats')
            if (reuse_stats_selected) error stop 'Duplicate --reuse-stats option'
            call read_argument_text(argument_index+1, reuse_filename)
            reuse_stats_selected = .true.
         case ('--threads', '--pair-chunk', '--pair-buffer-mib', '--reuse-cache-mib')
            block
               character(len=:), allocatable :: argument_value
               call read_argument_text(argument_index+1, argument_value)
               read(argument_value, *, iostat=integer_status) option_value
            end block
            if (integer_status /= 0) error stop 'Memory/thread option requires an integer'
            select case (option_name)
            case ('--threads')
               if (option_value < 1) error stop '--threads must be positive'
               requested_threads = option_value
            case ('--pair-chunk')
               if (option_value < 1) error stop '--pair-chunk must be positive'
               pair_chunk_size = option_value
            case ('--pair-buffer-mib')
               if (option_value < 0) error stop '--pair-buffer-mib must be nonnegative'
               result_memory_mib = option_value
            case ('--reuse-cache-mib')
               if (reuse_selected) error stop 'Duplicate --reuse-cache-mib option'
               if (option_value < 0 .or. option_value > 4096) error stop '--reuse-cache-mib must be in 0..4096'
               reuse_memory_mib = option_value
               reuse_selected = .true.
            end select
         case ('--help')
            print '(A)', 'rangular [--csf-text FILE | --csf-descriptors FILE] [--timings NEW_FILE.csv]'
            print '(A)', 'Default: original rcsf.inp text input; interactive scientific settings are unchanged.'
            print '(A)', 'Serial only: --reuse-cache-mib N [--reuse-stats NEW_FILE.csv]; default 0 (off).'
            if (parallel_pairs_available) then
               print '(A)', 'rangular_omp: --threads N --pair-chunk N --pair-buffer-mib N'
               print '(A)', 'Defaults: OMP_NUM_THREADS, 256 pairs/chunk, 64 MiB result buffers (0 forces spill).'
            end if
            stop
         case default
            error stop 'Unknown rangular option; use --help'
         end select
         argument_index = argument_index + 2
      end do
      call reject_output_filename(options%input_filename)
      call reject_output_filename(options%timings_filename)
      call reject_output_filename(reuse_filename)
      if (len(reuse_filename) > 0) then
         if (reuse_filename == options%input_filename .or. reuse_filename == options%timings_filename) &
            error stop 'Reuse statistics must use a separate filename'
      end if
      if (parallel_pairs_available .and. (reuse_selected .or. reuse_stats_selected)) &
         error stop 'Angular reuse is currently supported by serial rangular only'
      call configure_pair_engine(requested_threads, pair_chunk_size, result_memory_mib)
      call configure_angular_reuse(reuse_memory_mib,reuse_filename)
      call open_stage_profile(options%timings_filename)
   end subroutine parse_angular_options

   subroutine load_angular_header(options, core_count, maximum_blocks, block_labels)
      type(angular_input_options), intent(inout) :: options
      integer, intent(out) :: core_count
      integer, intent(in) :: maximum_blocks
      character(len=8), intent(out) :: block_labels(:)
      if (options%use_descriptors) then
         call open_descriptor_source(options%descriptors, options%input_filename, maximum_blocks)
         call install_descriptor_header(options%descriptors, core_count, block_labels)
      else
         call cslh(options%input_filename, core_count, maximum_blocks, block_labels)
      end if
   end subroutine load_angular_header

   subroutine load_angular_block(options, core_count, block_index)
      use orb_c, only: occupations => iqa
      use stat_c, only: subshell_quantum_numbers => jqsa, couplings => jcupa
      type(angular_input_options), intent(inout) :: options
      integer, intent(in) :: core_count, block_index
      if (options%use_descriptors) then
         call read_descriptor_block(options%descriptors, block_index, occupations, subshell_quantum_numbers, couplings)
      else
         call lodcsh2(21, core_count, block_index)
      end if
   end subroutine load_angular_block

   subroutine close_angular_input(options)
      type(angular_input_options), intent(inout) :: options
      if (options%use_descriptors) then
         call close_descriptor_source(options%descriptors)
      else
         close(21)
      end if
   end subroutine close_angular_input
end module rangular_input
