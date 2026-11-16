program npol_consumers
    use, intrinsic :: iso_fortran_env, only: real64
    use parameter_def, only: NNNP
    use npol_test_support
    use npol_C
    use npol_initialization
    use debug_C, only: LDBPR
    use def_c, only: C => speed_of_light_au
    use grid_c, only: N => radial_point_count, R => radial_grid, RP => radial_grid_derivative, RPOR => &
       & radial_derivative_over_radius, H => transformed_grid_step
    use npot_c, only: ZZ => negative_radius_times_potential
    use orb_c, only: NAK => orbital_kappa
    use pote_c, only: YP, XP, XQ
    use scf_c, only: NYCOF, NDCOF, NDA, DA
    use wave_c, only: MF => orbital_last_grid_index, PF => radial_large_components, QF => radial_small_components
    use ypot_I
    use dacon_I
    use rinti_I
#ifdef NPOL_MPI
    use mpi
    use mpi_c, only: myid => mpi_rank, nprocs => mpi_process_count
    use npol_mpi, only: npol_init_mpi
#endif
    implicit none (type, external)
    real(real64) :: coul(NNNP), base_xp(NNNP), base_xq(NNNP), base_yp(NNNP), reduced(NNNP)
    real(real64) :: initial(2,2), delta, expected, actual, term, q
    integer :: i, j, mode, status, ranks, rank, ierr
    character(NPOL_MESSAGE_LEN) :: message
    ranks = 1
    rank = 0
#ifdef NPOL_MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    ranks = nprocs
    rank = myid
#endif
    N = 4
    R = 0
    R(:N) = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]
    RP = 1
    RPOR = 0.5_real64
    H = 1
    C = 2
    LDBPR = .false.
    NYCOF = 0
    NDCOF = 1
    allocate(NDA(1), DA(1), PF(NNNP,2), QF(NNNP,2))
    NDA = 2
    q = 0.25_real64
    DA = q
    NAK = 1
    MF(:2) = N
    PF = 0
    QF = 0
    ! Opposite-sign small components cancel the diagonal external pieces in
    ! the second fixture below while leaving an off-diagonal overlap.
    PF(2:N,1) = [1.0_real64, 2.0_real64, 1.0_real64]
    QF(2:N,1) = [0.5_real64, 1.0_real64, 0.5_real64]
    PF(2:N,2) = [2.0_real64, 1.0_real64, 3.0_real64]
    QF(2:N,2) = [1.0_real64, 0.5_real64, 1.5_real64]
    coul = 7
    ZZ = coul
    call snapshot()
    ! Absent configuration goes through the real initialization seam.
#ifdef NPOL_MPI
    call npol_init_mpi("missing.inp", 2.0_real64, N, R, ZZ, MPI_COMM_WORLD, status, message)
#else
    call npol_init("missing.inp", 2.0_real64, N, R, ZZ, status, message)
#endif
    call require(status == 0 .and. all(ZZ == coul), "disabled preserves full ZZ exactly")
    call YPOT(1)
    call require(all(YP(:N) == base_yp(:N)), "disabled YPOT unchanged")
    XP = 0
    XQ = 0
    call DACON()
    call require(all(XP == base_xp) .and. all(XQ == base_xq), "disabled DACON unchanged")
    do j = 1, 2
        do mode = 0, 1
            actual = RINTI(1,j,mode)
            call require(actual == initial(j,mode+1), "disabled RINTI unchanged")
        end do
    end do
    if (rank == 0) call write_fixture("npol.inp")
#ifdef NPOL_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call npol_init_mpi("npol.inp", 2.0_real64, N, R, ZZ, MPI_COMM_WORLD, status, message)
#else
    call npol_init("npol.inp", 2.0_real64, N, R, ZZ, status, message)
#endif
    call require(status == 0, message)
    call YPOT(1)
    call require(all(YP(:N) == ZZ(:N)/ranks), "production YPOT contribution")
#ifdef NPOL_MPI
    call MPI_Allreduce(YP, reduced, N, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    reduced(:N) = YP(:N)
#endif
    do i = 1, N
        call near(reduced(i), ZZ(i), "reduced YPOT contains one full external potential")
    end do
    XP = 0
    XQ = 0
    call DACON()
    do i = 2, N
        delta = ZZ(i)-coul(i)
        call near(XP(i)-base_xp(i), -q*delta*QF(i,2)/C, "production DACON XP sign and scale")
        call near(XQ(i)-base_xq(i), q*delta*PF(i,2)/C, "production DACON XQ sign and scale")
    end do
    do j = 1, 2
        expected = 0
        do i = 2, N
            term = real(i,real64)*RPOR(i)*(ZZ(i)-coul(i))*(PF(i,1)*PF(i,j)+QF(i,1)*QF(i,j))
            expected = expected-term
        end do
        actual = RINTI(1,j,0)
        call near(actual-initial(j,1), expected, "production RINTI MODE=0 diagonal/offdiagonal")
        actual = RINTI(1,j,1)
        call require(actual == initial(j,2), "production RINTI MODE=1 unchanged")
    end do
    ! Isolate only a non-diagonal matrix element using a signed external
    ! perturbation: diagonal weighted overlaps cancel exactly.
    PF = 0
    QF = 0
    PF(2:3,1) = 1
    PF(2:3,2) = [1.0_real64, -1.0_real64]
    ZZ = coul
    call snapshot()
    ZZ(2) = coul(2)+3
    ZZ(3) = coul(3)-2
    actual = RINTI(1,1,0)
    call near(actual-initial(1,1), 0.0_real64, "isolated fixture diagonal cancels")
    actual = RINTI(1,2,0)
    call near(actual-initial(2,1), -6.0_real64, "isolated offdiagonal external integral")
#ifdef NPOL_MPI
    call MPI_Finalize(ierr)
#endif
    print "(A)", "PASS production consumers"
contains
    subroutine snapshot()
        integer :: jj, mm
        call YPOT(1)
        base_yp = YP
        XP = 0
        XQ = 0
        call DACON()
        base_xp = XP
        base_xq = XQ
        do jj = 1, 2
            do mm = 0, 1
                initial(jj,mm+1) = RINTI(1,jj,mm)
            end do
        end do
    end subroutine snapshot
end program npol_consumers

! Deterministic test-only derivative/quadrature seams. Never linked into apps.
