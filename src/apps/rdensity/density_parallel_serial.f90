! Serial implementation of the rdensity process boundary.
module density_parallel
    implicit none (type, external)

    ! Module state
    integer :: mpi_rank = 0, mpi_process_count = 1
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
        write(*,*) "RDENSITY: serial density and natural orbitals"
    end subroutine parallel_start

    subroutine parallel_finish
        write(*,*) "RDENSITY: Execution complete."
    end subroutine parallel_finish

    subroutine broadcast_options(name, nci, ndef, ordering, status)
        ! Arguments
        character(*), intent(inout) :: name
        integer, intent(inout) :: nci, ndef, ordering, status
    end subroutine broadcast_options

    subroutine sum_density(matrix)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        real(real64), contiguous, intent(inout) :: matrix(:,:)
    end subroutine sum_density

    subroutine density_fail(message)
        ! Arguments
        character(*), intent(in) :: message
        write(*,*) "RDENSITY: ", message
        error stop 1
    end subroutine density_fail
end module density_parallel
