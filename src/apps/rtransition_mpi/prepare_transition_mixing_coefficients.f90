! Make an independent copy for same-state transitions, preserving the file encoding.
module prepare_transition_mixing_coefficients_interface
    implicit none (type, external)
    private
    public :: prepare_transition_mixing_coefficients
contains
    subroutine prepare_transition_mixing_coefficients(name, inpci)
        use, intrinsic :: iso_fortran_env, only: int8, int64
        character(*), intent(in) :: name(2)
        integer, intent(in) :: inpci
        integer :: source, copied, count
        integer(int64) :: bytes, offset
        integer(int8), allocatable :: buffer(:)
        character(:), allocatable :: suffix

        suffix = trim(merge(".cbm", ".bm ", inpci == 0))
        open(newunit=source, file=trim(name(1))//suffix, form="unformatted", access="stream", action="read", status="old")
        open(newunit=copied, file=trim(name(1))//"_CP"//suffix, form="unformatted", access="stream", &
            action="write", status="replace")
        inquire(unit=source, size=bytes)
        allocate(buffer(1024*1024))
        do offset = 1_int64, bytes, size(buffer)
            count = int(min(int(size(buffer),int64), bytes-offset+1))
            read(source) buffer(:count)
            write(copied) buffer(:count)
        end do
        close(source)
        close(copied)
    end subroutine
end module prepare_transition_mixing_coefficients_interface
