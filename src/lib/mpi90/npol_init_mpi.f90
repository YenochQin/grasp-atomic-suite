! This is the only NP module with an MPI dependency.
module npol_mpi
    use mpi, only: MPI_Comm_rank, MPI_Bcast, MPI_Allreduce, MPI_Abort, &
        MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_CHARACTER, MPI_MIN
    use vast_kind_param, only: double
    use parameter_def, only: NNNP
    use npol_C
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
        real(double) :: a(NNNP), b(NNNP), v(NNNP), candidate(NNNP)
        integer :: rank, ierr, root_n
        ! Fixed-size wire message, independent of caller's character length.
        character(NPOL_MESSAGE_LEN) :: detail
        call MPI_Comm_rank(comm, rank, ierr)
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
        call MPI_Bcast(root_n, 1, MPI_INTEGER, 0, comm, ierr)
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
        call MPI_Bcast(a, n, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(b, n, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(v, n, MPI_DOUBLE_PRECISION, 0, comm, ierr)
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
        character(NPOL_MESSAGE_LEN), intent(inout) :: message
        integer :: failed_rank, first_rank, ierr
        failed_rank = huge(0)
        if (status /= 0) failed_rank = rank
        call MPI_Allreduce(failed_rank, first_rank, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
        status = 0
        if (first_rank == huge(0)) then
            message = ''
            return
        end if
        status = 1
        call MPI_Bcast(message, len(message), MPI_CHARACTER, first_rank, comm, ierr)
    end subroutine

    subroutine broadcast_config(cfg, comm)
        type(npol_config), intent(inout) :: cfg
        integer, intent(in) :: comm
        integer :: ierr
        call MPI_Bcast(cfg%schema_version, 1, MPI_INTEGER, 0, comm, ierr)
        call MPI_Bcast(cfg%isotope_z, 1, MPI_INTEGER, 0, comm, ierr)
        call MPI_Bcast(cfg%isotope_a, 1, MPI_INTEGER, 0, comm, ierr)
        call MPI_Bcast(cfg%enabled, 1, MPI_LOGICAL, 0, comm, ierr)
        call MPI_Bcast(cfg%use_e1, 1, MPI_LOGICAL, 0, comm, ierr)
        call MPI_Bcast(cfg%use_e2, 1, MPI_LOGICAL, 0, comm, ierr)
        call MPI_Bcast(cfg%k1, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(cfg%b1, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(cfg%k2, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(cfg%b2, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(cfg%units, len(cfg%units), MPI_CHARACTER, 0, comm, ierr)
        call MPI_Bcast(cfg%model_id, len(cfg%model_id), MPI_CHARACTER, 0, comm, ierr)
        call MPI_Bcast(cfg%parameter_source, len(cfg%parameter_source), MPI_CHARACTER, 0, comm, ierr)
    end subroutine

    subroutine npol_abort_mpi(comm, message)
        integer, intent(in) :: comm
        character(*), intent(in) :: message
        integer :: rank, ierr
        call MPI_Comm_rank(comm, rank, ierr)
        if (rank == 0) write(*,'(A)') trim(message)
        call MPI_Abort(comm, 1, ierr)
    end subroutine
end module npol_mpi
