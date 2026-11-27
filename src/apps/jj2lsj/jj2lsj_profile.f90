! Collect optional wall-clock timings for the conversion and label-writing stages.
! Call profile_emit from the coordinator thread.
module jj2lsj_profile
    ! Module dependencies
    use iso_fortran_env, only: int64, real64, error_unit
    implicit none (type, external)

    ! Module interface
    private
    public :: profile_start, profile_open, profile_finish, profile_time, profile_emit
    public :: profile_file_size

    ! Module state
    logical :: enabled = .false.
    integer :: stream
    integer(int64) :: rate
    real(real64) :: beginning
contains
    subroutine profile_start()
        ! Local declarations
        character(len=16) :: setting
        integer :: status
        call get_environment_variable("JJ2LSJ_PROFILE", setting, status=status)
        enabled = .false.
        if (status == 0) enabled = trim(setting) == "1"
        if (.not. enabled) return
        call system_clock(count_rate=rate)
        if (rate <= 0) error stop "profile_start: wall clock unavailable"
        beginning = profile_time()
    end subroutine profile_start

    subroutine profile_open(path)
        character(*), intent(in) :: path
        integer :: status

        if (.not. enabled) return
        open(newunit=stream, file=path, status="replace", action="write", iostat=status)
        if (status /= 0) then
            write(error_unit,"(A)") "Cannot write "//path
            error stop 1
        end if
        write(stream,"(A)") "phase,block,level,items,seconds"
    end subroutine profile_open

    function profile_time() result(seconds)
        ! Arguments
        real(real64) :: seconds

        ! Local declarations
        integer(int64) :: count
        seconds = 0.0_real64
        if (.not. enabled) return
        call system_clock(count)
        seconds = real(count,real64)/real(rate,real64)
    end function profile_time

    subroutine profile_emit(phase, started, block, level, items, elapsed)
        ! Arguments
        character(len=*), intent(in) :: phase
        real(real64), intent(in) :: started
        integer, intent(in), optional :: block, level, items
        real(real64), intent(in), optional :: elapsed

        ! Local declarations
        integer :: b, l, n
        real(real64) :: duration
        if (.not. enabled) return
        duration = profile_time()-started
        if (present(elapsed)) duration = elapsed
        b=0; l=0; n=0
        if (present(block)) b=block
        if (present(level)) l=level
        if (present(items)) n=items
        write(stream,'(A,3(",",I0),",",F18.9)') phase,b,l,n,duration
        flush(stream)
    end subroutine profile_emit

    subroutine profile_file_size(phase, unit)
        ! Arguments
        character(len=*), intent(in) :: phase
        integer, intent(in) :: unit

        ! Local declarations
        integer(int64) :: bytes
        integer :: status
        real(real64) :: started
        if (.not. enabled) return
        started=profile_time()
        flush(unit,iostat=status)
        if (status/=0) error stop "profile_file_size: output flush failed"
        call profile_emit("output_flush",started)
        inquire(unit=unit,size=bytes,iostat=status)
        if (status/=0 .or. bytes<0) error stop "profile_file_size: output size unavailable"
        write(stream,'(A,3(",",I0),",",F18.9)') phase,0,0,bytes,0.0_real64
    end subroutine profile_file_size

    subroutine profile_finish()
        if (.not. enabled) return
        call profile_emit("total",beginning)
        close(stream)
        enabled=.false.
    end subroutine profile_finish
end module jj2lsj_profile
