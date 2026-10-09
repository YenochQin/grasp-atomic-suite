! Serial implementation of the rdensity process boundary.
module density_parallel
    use vast_kind_param, only: DOUBLE
    implicit none
    integer :: rank = 0, workers = 1
contains
    subroutine parallel_start
        write(*,*) 'RDENSITY: serial density and natural orbitals'
    end subroutine

    subroutine parallel_finish
        write(*,*) 'RDENSITY: Execution complete.'
    end subroutine

    subroutine broadcast_options(name, nci, ndef, ordering, status)
        character(*), intent(inout) :: name
        integer, intent(inout) :: nci, ndef, ordering, status
    end subroutine

    subroutine sum_density(matrix)
        real(DOUBLE), contiguous, intent(inout) :: matrix(:,:)
    end subroutine

    subroutine density_fail(message)
        character(*), intent(in) :: message
        write(*,*) 'RDENSITY: ', message
        error stop 1
    end subroutine
end module density_parallel
