! This is the only NP module with an MPI dependency.
module npol_mpi
    use mpi, only: mpi_comm_rank, mpi_bcast, mpi_allreduce, mpi_abort, &
        mpi_integer, mpi_double_precision, mpi_logical, mpi_character, mpi_min
    use vast_kind_param, only: double
    use parameter_def, only: max_radial_points
    use npol_c
    use npol_input, only: npol_read_config, npol_validate_config
    use npol_kernel, only: npol_check_grid, npol_build_grid, npol_compose
    use npol_initialization, only: npol_write_summary, npol_summary_path
    implicit none
    private
    public :: npol_init_mpi, npol_abort_mpi
contains
    subroutine npol_init_mpi(path, z, n, r, zz, comm, status, message)
        character(*), intent(in) :: path
        real(double), intent(in) :: z, r(:)
        integer, intent(in) :: n, comm
        real(double), intent(inout) :: zz(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message
        type(npol_config) :: cfg
        real(double) :: a(max_radial_points), b(max_radial_points), v(max_radial_points), candidate(max_radial_points)
        integer :: rank, ierr, root_n
        ! Fixed-size wire message, independent of caller's character length.
        character(npol_message_len) :: detail
        call mpi_comm_rank(comm, rank, ierr)
        call npol_reset()
        call npol_check_grid(n, r, status, detail)
        if (size(zz) < n) then
            status = 1
            detail = 'npol: Coulomb array shorter than N'
        end if
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        root_n = n
        call mpi_bcast(root_n, 1, mpi_integer, 0, comm, ierr)
        if (root_n /= n) then
            status = 1
            detail = 'npol: inconsistent active grid length across ranks'
        end if
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        call npol_capture(n, zz(:n))

        if (rank == 0) then
            call npol_read_config(path, cfg, status, detail)
            if (status == 0) call npol_validate_config(cfg, z, status, detail)
        end if
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        call broadcast_config(cfg, comm)
        ! Also reject a mismatched nuclear model on a worker.
        call npol_validate_config(cfg, z, status, detail)
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        if (rank == 0) call npol_build_grid(cfg, n, r, a, b, v, status, detail)
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        call mpi_bcast(a, n, mpi_double_precision, 0, comm, ierr)
        call mpi_bcast(b, n, mpi_double_precision, 0, comm, ierr)
        call mpi_bcast(v, n, mpi_double_precision, 0, comm, ierr)
        call npol_compose(n, r, zz_coul, v, candidate, status, detail)
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        if (rank == 0) call npol_write_summary(npol_summary_path(path), cfg, n, status, detail)
        call sync_status(comm, rank, status, detail)
        message = detail
        if (status /= 0) return
        ! Every rank keeps the full external potential. YPOT alone apportions
        ! its rank-local contribution for the existing COFPOT reduction.
        zz(:n) = candidate(:n)
        call npol_publish(cfg, n, a(:n), b(:n), v(:n))
    end subroutine

    subroutine sync_status(comm, rank, status, message)
        integer, intent(in) :: comm, rank
        integer, intent(inout) :: status
        character(npol_message_len), intent(inout) :: message
        integer :: failed_rank, first_rank, ierr
        failed_rank = huge(0)
        if (status /= 0) failed_rank = rank
        call mpi_allreduce(failed_rank, first_rank, 1, mpi_integer, mpi_min, comm, ierr)
        status = 0
        if (first_rank == huge(0)) then
            message = ''
            return
        end if
        status = 1
        call mpi_bcast(message, len(message), mpi_character, first_rank, comm, ierr)
    end subroutine

    subroutine broadcast_config(cfg, comm)
        type(npol_config), intent(inout) :: cfg
        integer, intent(in) :: comm
        integer :: ierr
        call mpi_bcast(cfg%schema_version, 1, mpi_integer, 0, comm, ierr)
        call mpi_bcast(cfg%isotope_z, 1, mpi_integer, 0, comm, ierr)
        call mpi_bcast(cfg%isotope_a, 1, mpi_integer, 0, comm, ierr)
        call mpi_bcast(cfg%enabled, 1, mpi_logical, 0, comm, ierr)
        call mpi_bcast(cfg%use_e1, 1, mpi_logical, 0, comm, ierr)
        call mpi_bcast(cfg%use_e2, 1, mpi_logical, 0, comm, ierr)
        call mpi_bcast(cfg%k1, 1, mpi_double_precision, 0, comm, ierr)
        call mpi_bcast(cfg%b1, 1, mpi_double_precision, 0, comm, ierr)
        call mpi_bcast(cfg%k2, 1, mpi_double_precision, 0, comm, ierr)
        call mpi_bcast(cfg%b2, 1, mpi_double_precision, 0, comm, ierr)
        call mpi_bcast(cfg%units, len(cfg%units), mpi_character, 0, comm, ierr)
        call mpi_bcast(cfg%model_id, len(cfg%model_id), mpi_character, 0, comm, ierr)
        call mpi_bcast(cfg%parameter_source, len(cfg%parameter_source), mpi_character, 0, comm, ierr)
    end subroutine

    subroutine npol_abort_mpi(comm, message)
        integer, intent(in) :: comm
        character(*), intent(in) :: message
        integer :: rank, ierr
        call mpi_comm_rank(comm, rank, ierr)
        if (rank == 0) write(*,'(A)') trim(message)
        call mpi_abort(comm, 1, ierr)
    end subroutine
end module npol_mpi
