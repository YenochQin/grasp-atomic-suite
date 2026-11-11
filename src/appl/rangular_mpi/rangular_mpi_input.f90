! Rank zero owns the descriptor stream; all ranks retain the ordinary CSF arrays.
module rangular_mpi_input
   use iso_fortran_env, only: int64
   use mpi_C, only: rank_index => myid, MPI_Bcast, MPI_INTEGER, MPI_INTEGER1, MPI_CHARACTER, MPI_COMM_WORLD
   use parameter_def, only: orbital_capacity => NNNW
   use orb_C, only: current_csf_count => NCF, orbital_count => NW, occupations => IQA, &
                    principal_numbers => NP, relativistic_kappas => NAK, orbital_momenta => NKL, &
                    single_electron_two_j => NKJ, orbital_labels => NH
   use stat_C, only: subshell_quantum_numbers => JQSA, couplings => JCUPA
   use hblock_C, only: block_count => NBLOCK, block_lengths => NCFBLK
   use def_C, only: electron_count => NELEC
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, install_descriptor_header, &
                                   read_descriptor_block, close_descriptor_source
   use cslhmpi_I, only: cslhmpi
   use lodcslmpi_I, only: lodcslmpi
   use rangular_profile, only: open_stage_profile
   implicit none
   private
   public :: mpi_angular_input_options, parse_mpi_angular_options, open_mpi_angular_profile
   public :: load_mpi_angular_header, load_mpi_angular_block, close_mpi_angular_input

   type :: mpi_angular_input_options
      character(len=:), allocatable :: input_filename, timings_filename
      logical :: use_descriptors = .false.
      type(descriptor_source) :: descriptors
   end type mpi_angular_input_options
contains
   subroutine read_argument_text(argument_index, argument_text)
      integer, intent(in) :: argument_index
      character(len=:), allocatable, intent(out) :: argument_text
      integer :: text_length, argument_status
      call get_command_argument(argument_index, length=text_length, status=argument_status)
      if (argument_status /= 0 .or. text_length == 0) error stop 'Missing rangular_mpi argument'
      allocate(character(len=text_length) :: argument_text)
      call get_command_argument(argument_index, value=argument_text, status=argument_status)
      if (argument_status /= 0) error stop 'Cannot read rangular_mpi argument'
   end subroutine read_argument_text

   subroutine reject_reserved_filename(filename)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: basename
      integer :: separator_position
      separator_position = scan(filename, '/\', back=.true.)
      basename = filename(separator_position+1:)
      if (basename == 'rangular.log' .or. basename == 'genmcp.dbg' .or. &
          index(basename,'mcp') == 1 .or. index(basename,'tmp') == 1) &
         error stop 'Input/timings filename uses a reserved rangular_mpi output name'
   end subroutine reject_reserved_filename

   subroutine parse_mpi_angular_options(options)
      type(mpi_angular_input_options), intent(out) :: options
      integer :: argument_index
      logical :: input_selected, timings_selected
      character(len=:), allocatable :: option_name
      options%input_filename = 'rcsf.inp'
      options%timings_filename = ''
      input_selected = .false.
      timings_selected = .false.
      argument_index = 1
      do while (argument_index <= command_argument_count())
         call read_argument_text(argument_index,option_name)
         select case (option_name)
         case ('--csf-text','--csf-descriptors')
            if (input_selected) error stop 'Select exactly one CSF input option'
            call read_argument_text(argument_index+1,options%input_filename)
            options%use_descriptors = option_name == '--csf-descriptors'
            input_selected = .true.
         case ('--timings')
            if (timings_selected) error stop 'Duplicate --timings option'
            call read_argument_text(argument_index+1,options%timings_filename)
            timings_selected = .true.
         case ('--help')
            print '(A)', 'rangular_mpi [--csf-text FILE | --csf-descriptors FILE] [--timings NEW_FILE.csv]'
            print '(A)', 'Default: original rcsf.inp input in the MPI serial I/O directory.'
            print '(A)', 'Only rank zero reads CSFs. Relative input paths use the MPI serial I/O directory.'
            print '(A)', 'Timings: NEW_FILE.rank000.csv, etc.; relative paths use each MPI rank work directory.'
            stop
         case default
            error stop 'Unknown rangular_mpi option; use --help'
         end select
         argument_index = argument_index + 2
      end do
      call reject_reserved_filename(options%input_filename)
      call reject_reserved_filename(options%timings_filename)
   end subroutine parse_mpi_angular_options

   subroutine open_mpi_angular_profile(options)
      type(mpi_angular_input_options), intent(in) :: options
      character(len=:), allocatable :: profile_prefix
      character(len=16) :: rank_label
      integer :: prefix_length
      if (len(options%timings_filename) == 0) return
      profile_prefix = options%timings_filename
      prefix_length = len(profile_prefix)
      if (prefix_length >= 4) then
         if (profile_prefix(prefix_length-3:) == '.csv') profile_prefix = profile_prefix(:prefix_length-4)
      end if
      write(rank_label,'(I3.3)') rank_index
      call open_stage_profile(profile_prefix//'.rank'//trim(rank_label)//'.csv')
   end subroutine open_mpi_angular_profile

   subroutine load_mpi_angular_header(options, serial_directory, core_count, maximum_blocks, block_labels)
      type(mpi_angular_input_options), intent(inout) :: options
      character(len=*), intent(in) :: serial_directory
      integer, intent(out) :: core_count
      integer, intent(in) :: maximum_blocks
      character(len=8), intent(out) :: block_labels(:)
      character(len=:), allocatable :: input_path
      integer :: mpi_status, allocation_status
      input_path = options%input_filename
      if (rank_index == 0 .and. input_path(1:1) /= '/') input_path = trim(serial_directory)//'/'//input_path
      if (.not. options%use_descriptors) then
         call cslhmpi(input_path,core_count,maximum_blocks,block_labels)
         return
      end if
      if (rank_index == 0) then
         call open_descriptor_source(options%descriptors,input_path,maximum_blocks)
         call install_descriptor_header(options%descriptors,core_count,block_labels)
      end if
      call MPI_Bcast(block_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(current_csf_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      if (rank_index /= 0) then
         allocate(block_lengths(block_count+1),stat=allocation_status)
         if (allocation_status /= 0) error stop 'Cannot allocate MPI CSF block lengths'
         block_labels = ''
      end if
      call MPI_Bcast(block_lengths,block_count+1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(block_labels,8*block_count,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(core_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(electron_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(orbital_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(principal_numbers,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(relativistic_kappas,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(orbital_momenta,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(single_electron_two_j,orbital_count,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_status)
      call MPI_Bcast(orbital_labels,2*orbital_count,MPI_CHARACTER,0,MPI_COMM_WORLD,mpi_status)
   end subroutine load_mpi_angular_header

   subroutine load_mpi_angular_block(options, core_count, block_index)
      type(mpi_angular_input_options), intent(inout) :: options
      integer, intent(in) :: core_count, block_index
      integer :: first_csf, last_csf, columns_per_message, message_columns, mpi_status
      if (.not. options%use_descriptors) then
         call lodcslmpi(21,core_count,block_index)
         return
      end if
      if (rank_index == 0) call read_descriptor_block(options%descriptors,block_index, &
                                                     occupations,subshell_quantum_numbers,couplings)
      ! Bound the largest (three-channel) count before narrowing to the MPI integer interface.
      columns_per_message = huge(columns_per_message)
      columns_per_message = columns_per_message/(3*orbital_capacity)
      first_csf = 1
      do
         last_csf = int(min(int(current_csf_count,int64),int(first_csf,int64)+columns_per_message-1_int64))
         message_columns = last_csf-first_csf+1
         call MPI_Bcast(occupations(:,first_csf:last_csf),orbital_capacity*message_columns, &
                        MPI_INTEGER1,0,MPI_COMM_WORLD,mpi_status)
         call MPI_Bcast(subshell_quantum_numbers(:,:,first_csf:last_csf),3*orbital_capacity*message_columns, &
                        MPI_INTEGER1,0,MPI_COMM_WORLD,mpi_status)
         call MPI_Bcast(couplings(:,first_csf:last_csf),orbital_capacity*message_columns, &
                        MPI_INTEGER1,0,MPI_COMM_WORLD,mpi_status)
         if (last_csf == current_csf_count) exit
         first_csf = last_csf+1
      end do
   end subroutine load_mpi_angular_block

   subroutine close_mpi_angular_input(options)
      type(mpi_angular_input_options), intent(inout) :: options
      if (rank_index /= 0) return
      if (options%use_descriptors) then
         call close_descriptor_source(options%descriptors)
      else
         close(21)
      end if
   end subroutine close_mpi_angular_input
end module rangular_mpi_input
