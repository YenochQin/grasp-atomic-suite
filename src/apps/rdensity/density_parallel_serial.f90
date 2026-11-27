! Provide the serial counterparts of the rdensity parallel operations.
! Serial implementation of the rdensity process boundary.
module density_parallel
    use, intrinsic :: iso_fortran_env, only: int64
    use run_timer, only: begin_run_timer, finish_run_timer
    implicit none (type, external)

    ! Module state
    integer :: mpi_rank = 0, mpi_process_count = 1
    integer(int64), private :: start_clock_count
contains
    subroutine load_density_basis(state_name, core_subshell_count)
        ! Module dependencies
        use csf_basis_loader, only: load_angular_csf_basis

        ! Arguments
        character(len=24), intent(in) :: state_name
        integer, intent(out) :: core_subshell_count

        call load_angular_csf_basis(state_name, core_subshell_count)
    end subroutine load_density_basis

    subroutine parallel_start
        call begin_run_timer (start_clock_count, "RDENSITY")
        write(*,*) "RDENSITY: serial density and natural orbitals"
    end subroutine parallel_start

    subroutine parallel_finish
        call finish_run_timer (start_clock_count, "RDENSITY")
    end subroutine parallel_finish

    subroutine broadcast_options(name, mixing_file_kind, nondefault_options_flag, ordering, status)
        ! Arguments
        character(*), intent(inout) :: name
        integer, intent(inout) :: mixing_file_kind, nondefault_options_flag, ordering, status
    end subroutine broadcast_options

    subroutine sum_density(matrix)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        real(real64), contiguous, intent(inout) :: matrix(:,:)
    end subroutine sum_density

    subroutine density_abort()
        error stop 1
    end subroutine density_abort

    subroutine density_fail(message)
        ! Arguments
        character(*), intent(in) :: message
        write(*,*) "RDENSITY: ", message
        error stop 1
    end subroutine density_fail
end module density_parallel
