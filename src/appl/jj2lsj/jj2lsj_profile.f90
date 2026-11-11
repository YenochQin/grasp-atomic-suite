! Optional wall-clock diagnostics. Calls to profile_emit are coordinator-only.
module jj2lsj_profile
    use iso_fortran_env, only: int64, real64, error_unit
    implicit none
    private
    public :: profile_start, profile_finish, profile_time, profile_emit
    public :: profile_file_size
    logical :: enabled = .false.
    integer :: stream
    integer(int64) :: rate
    real(real64) :: beginning
contains
    subroutine profile_start()
        character(len=16) :: setting
        integer :: status
        call get_environment_variable('JJ2LSJ_PROFILE', setting, status=status)
        enabled = .false.
        if (status == 0) enabled = trim(setting) == '1'
        if (.not. enabled) return
        call system_clock(count_rate=rate)
        if (rate <= 0) error stop 'jj2lsj: wall clock unavailable'
        open(newunit=stream, file='jj2lsj_omp.profile.csv', status='replace', action='write', iostat=status)
        if (status /= 0) then
            write(error_unit,*) 'Cannot write jj2lsj_omp.profile.csv'
            error stop 1
        end if
        write(stream,'(A)') 'phase,block,level,items,seconds'
        beginning = profile_time()
    end subroutine

    function profile_time() result(seconds)
        real(real64) :: seconds
        integer(int64) :: count
        seconds = 0.0_real64
        if (.not. enabled) return
        call system_clock(count)
        seconds = real(count,real64)/real(rate,real64)
    end function

    subroutine profile_emit(phase, started, block, level, items, elapsed)
        character(len=*), intent(in) :: phase
        real(real64), intent(in) :: started
        integer, intent(in), optional :: block, level, items
        real(real64), intent(in), optional :: elapsed
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
    end subroutine

    subroutine profile_file_size(phase, unit)
        character(len=*), intent(in) :: phase
        integer, intent(in) :: unit
        integer(int64) :: bytes
        integer :: status
        real(real64) :: started
        if (.not. enabled) return
        started=profile_time()
        flush(unit,iostat=status)
        if (status/=0) error stop 'jj2lsj: output flush failed'
        call profile_emit('output_flush',started)
        inquire(unit=unit,size=bytes,iostat=status)
        if (status/=0 .or. bytes<0) error stop 'jj2lsj: output size unavailable'
        write(stream,'(A,3(",",I0),",",F18.9)') phase,0,0,bytes,0.0_real64
    end subroutine

    subroutine profile_finish()
        if (.not. enabled) return
        call profile_emit('total',beginning)
        close(stream)
        enabled=.false.
    end subroutine
end module
