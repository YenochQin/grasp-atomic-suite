module rangular_profile
   use iso_fortran_env, only: int64, real64
   implicit none
   private
   public :: open_stage_profile, close_stage_profile, stage_clock, record_stage
   integer, save :: profile_unit = 0
   logical, save :: profile_enabled = .false.
contains
   subroutine open_stage_profile(filename)
      character(len=*), intent(in) :: filename
      integer :: open_status
      if (len_trim(filename) == 0) return
      open(newunit=profile_unit, file=filename, status='new', action='write', iostat=open_status)
      if (open_status /= 0) error stop 'Cannot create timings file (it must not already exist)'
      profile_enabled = .true.
      write(profile_unit, '(A)') 'stage,block,channel,wall_seconds,records,storage'
   end subroutine open_stage_profile

   function stage_clock() result(clock_count)
      integer(int64) :: clock_count
      clock_count = 0_int64
      if (profile_enabled) call system_clock(count=clock_count)
   end function stage_clock

   subroutine record_stage(stage_name, block_index, start_count, record_count, channel, storage)
      character(len=*), intent(in) :: stage_name
      integer, intent(in) :: block_index
      integer(int64), intent(in) :: start_count
      integer(int64), intent(in), optional :: record_count
      integer, intent(in), optional :: channel
      character(len=*), intent(in), optional :: storage
      integer(int64) :: finish_count, clock_rate, maximum_count, elapsed_count, reported_records
      integer :: reported_channel
      character(len=24) :: reported_storage
      if (.not. profile_enabled) return
      call system_clock(count=finish_count, count_rate=clock_rate, count_max=maximum_count)
      if (clock_rate <= 0) error stop 'Monotonic wall clock unavailable'
      elapsed_count = finish_count - start_count
      if (finish_count < start_count) elapsed_count = (maximum_count-start_count) + finish_count + 1_int64
      reported_records = -1_int64
      if (present(record_count)) reported_records = record_count
      reported_channel = -1
      if (present(channel)) reported_channel = channel
      reported_storage = 'unspecified'
      if (present(storage)) reported_storage = storage
      write(profile_unit, '(A,",",I0,",",I0,",",ES24.16E3,",",I0,",",A)') &
         stage_name, block_index, reported_channel, real(elapsed_count,real64)/real(clock_rate,real64), &
         reported_records, trim(reported_storage)
      flush(profile_unit)
   end subroutine record_stage

   subroutine close_stage_profile()
      if (profile_enabled) close(profile_unit)
      profile_enabled = .false.
   end subroutine close_stage_profile
end module rangular_profile
