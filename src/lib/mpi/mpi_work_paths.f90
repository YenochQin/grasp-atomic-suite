module mpi_work_paths
    implicit none (type, external)
    private
    public :: cpath
contains
      ! Select rank-local scratch directories and a serial I/O directory for rank 0.
      ! Return each rank's startup directory and the scratch root above its numbered subdirectory.
      subroutine cpath (startup_directory, serial_io_directory, rank_scratch_root)
      use mpi_state, only: mpi_abort, mpi_allreduce, mpi_barrier, mpi_character, mpi_comm_world, &
          directory_status => mpi_error_code, &
         & mpi_in_place, mpi_integer, mpi_process_count, mpi_rank, mpi_recv, mpi_send, mpi_status_buffer, mpi_sum
      use mpi_runtime, only: stopmpi
      use mpi_change_directory, only: sys_chdir
      use mpi_current_directory, only: sys_getwd
      use mpi_make_directory, only: sys_mkdir
      implicit none (type, external)
      character(len=*), intent(out):: startup_directory, serial_io_directory, rank_scratch_root
      integer :: directory_status_sum, mpi_return_code
      ! The rank subdirectory width matches the I3.3 format used below.
      integer, parameter:: scratch_path_capacity = 128, rank_directory_name_width = 3
      character(len=scratch_path_capacity) :: scratch_path
      character(len=rank_directory_name_width) :: rank_directory_name
      integer ::   scratch_path_length, destination_rank, startup_directory_length
      integer :: retry_delay_value, retry_attempt, delay_iteration
      logical :: disks_file_exists
      character(:), allocatable :: scratch_environment_path, user_name

      call sys_getwd (startup_directory,startup_directory_length)

      if (mpi_rank == 0) then
         inquire(file="disks",exist=disks_file_exists)
         if(disks_file_exists) then
            ! Read the serial I/O path, the rank 0 scratch root, then one scratch root per worker in rank order.
            open (unit=1001, file="disks", status="old")
            read (1001,*) serial_io_directory
            read (1001,*) rank_scratch_root
            do destination_rank = 1, mpi_process_count - 1
               read (1001,*) scratch_path
               call mpi_send (scratch_path, scratch_path_capacity, mpi_character, destination_rank, destination_rank, &
                              mpi_comm_world, mpi_return_code)
            end do
            scratch_path = rank_scratch_root
            close (1001)
         else
            serial_io_directory = startup_directory
            ! Use MPI_TMP for all ranks when disks is absent. Fall back to /scratch/$USER.
            if(.not.getenv_allocating("MPI_TMP", scratch_environment_path)) then
                print *, " warning041 cpath failed at getenv MPI_TMP, myid = ", mpi_rank
                print *, " warning041 cpath will now default to /scratch/USER "

                if(.not.getenv_allocating("USER", user_name)) then
                    print *, " error061 cpath failed at getenv USER, myid = ", mpi_rank
                    print *, " error061 cpath fatal error, calling MPI_ABORT"
                    call mpi_abort(mpi_comm_world, directory_status, mpi_return_code)
                end if

                scratch_environment_path = "/scratch/"//user_name
            end if
            if(len(scratch_environment_path) > scratch_path_capacity) then
                print *, " error062 cpath MPI_TMP too long, myid = ", mpi_rank
                print *, " error062 cpath fatal error, calling MPI_ABORT"
                call mpi_abort(mpi_comm_world, directory_status, mpi_return_code)
            end if

            rank_scratch_root = scratch_environment_path
            scratch_path = scratch_environment_path
            do destination_rank = 1, mpi_process_count - 1
               call mpi_send (scratch_path, scratch_path_capacity, mpi_character, destination_rank, destination_rank, &
                              mpi_comm_world, mpi_return_code)
            end do
         end if
      else
         if (mpi_process_count > 1) then
           call mpi_recv (scratch_path, scratch_path_capacity, mpi_character, 0, mpi_rank, &
                          mpi_comm_world, mpi_status_buffer, mpi_return_code)
         end if
      end if
      scratch_path_length = len_trim (scratch_path)

      ! Create the scratch root if the initial directory change fails.
      call sys_chdir (scratch_path, scratch_path_length, directory_status)
      if (directory_status/=0) call sys_mkdir (scratch_path, scratch_path_length, directory_status)
      if (directory_status/=0) then
        ! Retry the failed directory operation up to twice.
        do retry_attempt = 1, 2
          do delay_iteration = 1, mpi_rank+1
            retry_delay_value = 2
            retry_delay_value = retry_delay_value * int(sqrt(real(3*delay_iteration)))
            retry_delay_value = retry_delay_value / int(sqrt(real(5*delay_iteration)))
          end do

          call sys_mkdir (scratch_path, scratch_path_length, directory_status)
          if (directory_status==0) exit
        end do
        if (directory_status /= 0) then
        print *, " error141 cpath failed at sys_mkdir(disk), myid = ",&
                                                             mpi_rank
         goto 999
        end if
      end if
      call sys_chdir (scratch_path, scratch_path_length, directory_status)
      if (directory_status /= 0) then
        ! Retry the failed directory operation up to twice.
        do retry_attempt = 1, 2
          do delay_iteration = 1, mpi_rank+1
            retry_delay_value = 2
            retry_delay_value = retry_delay_value * int(sqrt(real(3*delay_iteration)))
            retry_delay_value = retry_delay_value / int(sqrt(real(5*delay_iteration)))
          end do

          retry_delay_value = int(1 + 4 * retry_attempt * real(mpi_rank) / real(mpi_process_count))
          call sys_chdir (scratch_path, scratch_path_length, directory_status)
          if (directory_status==0) exit
        end do
        if (directory_status /= 0) then
       print *, " error161 cpath failed at sys_chdir(disk); myid = ",&
                                                            mpi_rank
          goto 999
        end if
      end if

      ! Enter a numbered subdirectory for this rank's local files.
      write (rank_directory_name,"(I3.3)") mpi_rank
      call sys_chdir (rank_directory_name, rank_directory_name_width, directory_status)
      if (directory_status /= 0) then
         call sys_mkdir (rank_directory_name, rank_directory_name_width, directory_status)
        if (directory_status /= 0) then
          print *, " error181 cpath failed to make sub-dir " // rank_directory_name
          goto 999
        end if
        call sys_chdir (rank_directory_name, rank_directory_name_width, directory_status)
        if (directory_status /= 0) then
         print *, " error183 cpath failed to go to sub-dir " // rank_directory_name
          goto 999
        end if
      end if
  999 continue

      ! Use the sum of directory statuses to select the collective error path.
      directory_status_sum = directory_status
      call mpi_barrier(mpi_comm_world, mpi_return_code)
      call mpi_allreduce (mpi_in_place,directory_status_sum,1, &
                          mpi_integer, mpi_sum, mpi_comm_world, mpi_return_code)
      if (directory_status_sum /= 0) then
        print *, " error201 cpath ierr = ", directory_status, " myid = ", mpi_rank
        if (mpi_rank == 0) then
          print *, " error201 cpath terminated at 999, myid = ", mpi_rank
        end if
        call stopmpi (" error in cpath ", mpi_rank)
      end if

      rank_scratch_root = scratch_path
      return

contains

    ! Read an environment value into an exactly sized allocation. An absent value leaves it intact.
    function getenv_allocating(variable_name, value) result(is_defined)
        logical :: is_defined
        character(*), intent(in) :: variable_name
        character(:), allocatable, intent(inout) :: value
        integer :: value_length, environment_status
        character(1) :: probe_value

        ! The one-character probe obtains the length even if the value is truncated.
        call get_environment_variable(variable_name, probe_value, value_length, environment_status)
        if(environment_status > 0) then
            is_defined = .false.
            return
        end if
        if(allocated(value)) then
            if(len(value) /= value_length) deallocate(value)
        end if
        if(.not.allocated(value)) allocate(character(value_length) :: value)
        if(value_length > 0) then
            call get_environment_variable(variable_name, value, value_length, environment_status)
        end if
        is_defined = .true.
    end function getenv_allocating

end subroutine cpath
end module mpi_work_paths
