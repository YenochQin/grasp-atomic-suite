! Serial implementation of the rdensity process boundary.
module density_parallel
    use setcsla_i, only: load_density_basis => setcsla
    use, intrinsic :: iso_fortran_env, only: double => real64
    implicit none (type, external)
    integer :: rank = 0, workers = 1
contains
    subroutine parallel_start
        write(*,*) "RDENSITY: serial density and natural orbitals"
    end subroutine parallel_start

    subroutine parallel_finish
        write(*,*) "RDENSITY: Execution complete."
    end subroutine parallel_finish

    subroutine broadcast_options(name, nci, ndef, ordering, status)
        character(*), intent(inout) :: name
        integer, intent(inout) :: nci, ndef, ordering, status
    end subroutine broadcast_options

    subroutine sum_density(matrix)
        use, intrinsic :: iso_fortran_env, only: real64
        real(real64), contiguous, intent(inout) :: matrix(:,:)
    end subroutine sum_density

    subroutine density_fail(message)
        character(*), intent(in) :: message
        write(*,*) "RDENSITY: ", message
        error stop 1
    end subroutine density_fail
end module density_parallel
