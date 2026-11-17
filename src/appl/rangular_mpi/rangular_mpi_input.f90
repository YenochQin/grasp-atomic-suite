! Rank zero owns the descriptor stream; all ranks retain the ordinary CSF arrays.
module rangular_mpi_input
   use csf_descriptor_collectives, only: broadcast_descriptor_arrays, check_input_status => check_descriptor_status
   use mpi_c, only: mpi_allreduce, mpi_max, mpi_abort
   use mpi_c, only: rank_index => mpi_rank, mpi_bcast, mpi_integer, mpi_integer1, mpi_character, mpi_comm_world
   use parameter_def, only: orbital_capacity => max_orbitals
   use orb_c, only: current_csf_count => csf_count, orbital_count => orbital_count, occupations => csf_subshell_occupations, &
                    principal_numbers => orbital_principal_n, relativistic_kappas => orbital_kappa, orbital_momenta => orbital_l, &
                    single_electron_two_j => orbital_two_j, orbital_labels => orbital_symmetry_labels
   use stat_c, only: subshell_quantum_numbers => csf_subshell_quantum_codes, couplings => csf_coupling_quantum_codes
   use hblock_c, only: block_count => symmetry_block_count, block_lengths => csf_counts_by_block
   use def_c, only: electron_count => electron_count
   use csf_descriptor_input, only: descriptor_source, open_descriptor_source, install_descriptor_header, &
                                   read_descriptor_block, close_descriptor_source
   use rangular_profile, only: open_stage_profile
   implicit none (type, external)
   private
   public :: mpi_angular_input_options, parse_mpi_angular_options, open_mpi_angular_profile
   public :: load_mpi_angular_header, load_mpi_angular_block, close_mpi_angular_input

   type :: mpi_angular_input_options
      character(len=:), allocatable :: input_filename, timings_filename, header_filename
      type(descriptor_source) :: descriptors
   end type mpi_angular_input_options
contains


   subroutine read_argument_text(argument_index, argument_text)
      integer, intent(in) :: argument_index
      character(len=:), allocatable, intent(out) :: argument_text
      integer :: text_length, argument_status
      call get_command_argument(argument_index, length=text_length, status=argument_status)
      if (argument_status /= 0 .or. text_length == 0) error stop "Missing rangular_mpi argument"
      allocate(character(len=text_length) :: argument_text)
      call get_command_argument(argument_index, value=argument_text, status=argument_status)
      if (argument_status /= 0) error stop "Cannot read rangular_mpi argument"
   end subroutine read_argument_text

   subroutine reject_reserved_filename(filename)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: basename
      integer :: separator_position
      separator_position = scan(filename, "/\", back=.true.)
      basename = filename(separator_position+1:)
      if (basename == "rangular.log" .or. basename == "genmcp.dbg" .or. &
          index(basename,"mcp") == 1 .or. index(basename,"tmp") == 1) &
         error stop "Input/timings filename uses a reserved rangular_mpi output name"
   end subroutine reject_reserved_filename

   subroutine parse_mpi_angular_options(options)
      type(mpi_angular_input_options), intent(out) :: options
      integer :: argument_index
      logical :: input_selected, timings_selected
      character(len=:), allocatable :: option_name
      options%input_filename = "rcsf_descriptors.parquet"
      options%timings_filename = ""
      options%header_filename = ""
      input_selected = .false.
      timings_selected = .false.
      argument_index = 1
      do while (argument_index <= command_argument_count())
         call read_argument_text(argument_index,option_name)
         select case (option_name)
         case ("--csf-descriptors")
            if (input_selected) error stop "Select exactly one CSF input option"
            call read_argument_text(argument_index+1,options%input_filename)
            input_selected = .true.
         case ("--header")
            if (len(options%header_filename) > 0) error stop "Duplicate --header option"
            call read_argument_text(argument_index+1,options%header_filename)
         case ("--timings")
            if (timings_selected) error stop "Duplicate --timings option"
            call read_argument_text(argument_index+1,options%timings_filename)
            timings_selected = .true.
         case ("--help")
            print "(A)", "rangular_mpi [--csf-descriptors FILE] [--header FILE.toml] [--timings NEW_FILE.csv]"
            print "(A)", "Default: rcsf_descriptors.parquet in the MPI serial I/O directory."
            print "(A)", "Only rank zero reads CSFs. Relative input paths use the MPI serial I/O directory."
            print "(A)", "Timings: NEW_FILE.rank000.csv, etc.; relative paths use each MPI rank work directory."
            stop
         case default
            error stop "Unknown rangular_mpi option; use --help"
         end select
         argument_index = argument_index + 2
      end do
      call reject_reserved_filename(options%header_filename)
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
         if (profile_prefix(prefix_length-3:) == ".csv") profile_prefix = profile_prefix(:prefix_length-4)
      end if
      write(rank_label,"(I3.3)") rank_index
      call open_stage_profile(profile_prefix//".rank"//trim(rank_label)//".csv")
   end subroutine open_mpi_angular_profile

   subroutine load_mpi_angular_header(options, serial_directory, core_count, maximum_blocks, block_labels)
      type(mpi_angular_input_options), intent(inout) :: options
      character(len=*), intent(in) :: serial_directory
      integer, intent(out) :: core_count
      integer, intent(in) :: maximum_blocks
      character(len=8), intent(out) :: block_labels(:)
      character(len=:), allocatable :: input_path, header_path
      integer :: mpi_status, allocation_status, input_status
      header_path = options%header_filename
      if (rank_index == 0 .and. len(header_path) > 0) then
         if (header_path(1:1) /= "/") header_path = trim(serial_directory)//"/"//header_path
      end if
      input_path = options%input_filename
      if (rank_index == 0 .and. input_path(1:1) /= "/") input_path = trim(serial_directory)//"/"//input_path
      input_status = 0
      if (rank_index == 0) then
         call open_descriptor_source(options%descriptors,input_path,maximum_blocks,header_path,status=input_status)
         if (input_status == 0) call install_descriptor_header(options%descriptors,core_count,block_labels,input_status)
      end if
      call check_input_status(input_status)
      allocation_status = 0
      call mpi_bcast(block_count,1,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(current_csf_count,1,mpi_integer,0,mpi_comm_world,mpi_status)
      if (rank_index /= 0) then
         allocate(block_lengths(block_count+1),stat=allocation_status)

         block_labels = ""
      end if
      call check_input_status(allocation_status)
      call mpi_bcast(block_lengths,block_count+1,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(block_labels,8*block_count,mpi_character,0,mpi_comm_world,mpi_status)
      call mpi_bcast(core_count,1,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(electron_count,1,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(orbital_count,1,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(principal_numbers,orbital_count,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(relativistic_kappas,orbital_count,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(orbital_momenta,orbital_count,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(single_electron_two_j,orbital_count,mpi_integer,0,mpi_comm_world,mpi_status)
      call mpi_bcast(orbital_labels,2*orbital_count,mpi_character,0,mpi_comm_world,mpi_status)
   end subroutine load_mpi_angular_header

   subroutine load_mpi_angular_block(options, core_count, block_index)
      type(mpi_angular_input_options), intent(inout) :: options
      integer, intent(in) :: core_count, block_index
      integer :: input_status
      input_status = 0
      if (rank_index == 0) call read_descriptor_block(options%descriptors,block_index, &
                                                     occupations,subshell_quantum_numbers,couplings,input_status)
      call check_input_status(input_status)
      call broadcast_descriptor_arrays(occupations,subshell_quantum_numbers,couplings)
   end subroutine load_mpi_angular_block

   subroutine close_mpi_angular_input(options)
      type(mpi_angular_input_options), intent(inout) :: options
      if (rank_index /= 0) return
      call close_descriptor_source(options%descriptors)
   end subroutine close_mpi_angular_input
end module rangular_mpi_input
