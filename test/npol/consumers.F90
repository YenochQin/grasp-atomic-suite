program npol_consumers
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters, only: max_radial_points
    use npol_test_support, only: near, require, write_fixture
    use nuclear_polarization_state, only: npol_message_len
    use nuclear_polarization_initialization, only: npol_init
    use debug_output_settings, only: rmcdhf_debug_flags
    use atomic_calculation_state, only: speed_of_light_au
    use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_derivative, radial_derivative_over_radius, &
        transformed_grid_step
    use nuclear_potential_state, only: negative_radius_times_potential
    use orbital_csf_basis, only: orbital_kappa
    use scf_radial_potentials, only: orbital_y_potential, large_exchange_source, small_exchange_source
    use scf_iteration_state, only: y_potential_coefficient_count, direct_coefficient_count, direct_orbital_indices, &
        direct_coefficients
    use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
    use build_direct_radial_potential_interface, only: build_direct_radial_potential
    use apply_off_diagonal_integrals_interface, only: apply_off_diagonal_integrals
    use one_electron_radial_integral, only: integrate_one_electron_radial
#ifdef NPOL_MPI
    use mpi_f08, only: mpi_allreduce, mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_double_precision, &
       & mpi_finalize, mpi_init, mpi_sum
    use mpi_state, only: mpi_rank, mpi_process_count
    use mpi_nuclear_polarization, only: npol_init_mpi
#endif
    implicit none (type, external)
    real(real64) :: coul(max_radial_points), base_xp(max_radial_points), base_xq(max_radial_points), &
        base_yp(max_radial_points), reduced(max_radial_points)
    real(real64) :: initial(2,2), delta, expected, actual, term, q
    integer :: i, j, mode, status, ranks, rank, ierr
    character(NPOL_MESSAGE_LEN) :: message
    ranks = 1
    rank = 0
#ifdef NPOL_MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_process_count, ierr)
    ranks = mpi_process_count
    rank = mpi_rank
#endif
    radial_point_count = 4
    radial_grid = 0
    radial_grid(:radial_point_count) = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]
    radial_grid_derivative = 1
    radial_derivative_over_radius = 0.5_real64
    transformed_grid_step = 1
    speed_of_light_au = 2
    rmcdhf_debug_flags = .false.
    y_potential_coefficient_count = 0
    direct_coefficient_count = 1
    allocate(direct_orbital_indices(1), direct_coefficients(1), radial_large_components(max_radial_points,2), &
        radial_small_components(max_radial_points,2))
    direct_orbital_indices = 2
    q = 0.25_real64
    direct_coefficients = q
    orbital_kappa = 1
    orbital_last_grid_index(:2) = radial_point_count
    radial_large_components = 0
    radial_small_components = 0
    ! Opposite-sign small components cancel the diagonal external pieces in
    ! the second fixture below while leaving an off-diagonal overlap.
    radial_large_components(2:radial_point_count,1) = [1.0_real64, 2.0_real64, 1.0_real64]
    radial_small_components(2:radial_point_count,1) = [0.5_real64, 1.0_real64, 0.5_real64]
    radial_large_components(2:radial_point_count,2) = [2.0_real64, 1.0_real64, 3.0_real64]
    radial_small_components(2:radial_point_count,2) = [1.0_real64, 0.5_real64, 1.5_real64]
    coul = 7
    negative_radius_times_potential = coul
    call snapshot()
    ! Absent configuration goes through the real initialization seam.
#ifdef NPOL_MPI
    call npol_init_mpi("missing.inp", 2.0_real64, radial_point_count, radial_grid, negative_radius_times_potential, &
        MPI_COMM_WORLD, status, message)
#else
    call npol_init("missing.inp", 2.0_real64, radial_point_count, radial_grid, negative_radius_times_potential, status, message)
#endif
    call require(status == 0 .and. all(negative_radius_times_potential == coul), "disabled preserves full ZZ exactly")
    call build_direct_radial_potential(1)
    call require(all(orbital_y_potential(:radial_point_count) == base_yp(:radial_point_count)), "disabled YPOT unchanged")
    large_exchange_source = 0
    small_exchange_source = 0
    call apply_off_diagonal_integrals()
    call require(all(large_exchange_source == base_xp) .and. all(small_exchange_source == base_xq), "disabled DACON unchanged")
    do j = 1, 2
        do mode = 0, 1
            actual = integrate_one_electron_radial(1,j,mode)
            call require(actual == initial(j,mode+1), "disabled RINTI unchanged")
        end do
    end do
    if (rank == 0) call write_fixture("npol.inp")
#ifdef NPOL_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call npol_init_mpi("npol.inp", 2.0_real64, radial_point_count, radial_grid, negative_radius_times_potential, &
        MPI_COMM_WORLD, status, message)
#else
    call npol_init("npol.inp", 2.0_real64, radial_point_count, radial_grid, negative_radius_times_potential, status, message)
#endif
    call require(status == 0, message)
    call build_direct_radial_potential(1)
    call require(all(orbital_y_potential(:radial_point_count) == negative_radius_times_potential(:radial_point_count)/ranks), &
        "production YPOT contribution")
#ifdef NPOL_MPI
    call MPI_Allreduce(orbital_y_potential, reduced, radial_point_count, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    reduced(:radial_point_count) = orbital_y_potential(:radial_point_count)
#endif
    do i = 1, radial_point_count
        call near(reduced(i), negative_radius_times_potential(i), "reduced YPOT contains one full external potential")
    end do
    large_exchange_source = 0
    small_exchange_source = 0
    call apply_off_diagonal_integrals()
    do i = 2, radial_point_count
        delta = negative_radius_times_potential(i)-coul(i)
        call near(large_exchange_source(i)-base_xp(i), -q*delta*radial_small_components(i,2)/speed_of_light_au, &
            "production DACON XP sign and scale")
        call near(small_exchange_source(i)-base_xq(i), q*delta*radial_large_components(i,2)/speed_of_light_au, &
            "production DACON XQ sign and scale")
    end do
    do j = 1, 2
        expected = 0
        do i = 2, radial_point_count
            term = real(i, &
                real64)*radial_derivative_over_radius(i)*(negative_radius_times_potential(i)-coul(i))* &
                (radial_large_components(i,1)*radial_large_components(i,j)+ &
                radial_small_components(i,1)*radial_small_components(i,j))
            expected = expected-term
        end do
        actual = integrate_one_electron_radial(1,j,0)
        call near(actual-initial(j,1), expected, "production RINTI MODE=0 diagonal/offdiagonal")
        actual = integrate_one_electron_radial(1,j,1)
        call require(actual == initial(j,2), "production RINTI MODE=1 unchanged")
    end do
    ! Isolate only a non-diagonal matrix element using a signed external
    ! perturbation: diagonal weighted overlaps cancel exactly.
    radial_large_components = 0
    radial_small_components = 0
    radial_large_components(2:3,1) = 1
    radial_large_components(2:3,2) = [1.0_real64, -1.0_real64]
    negative_radius_times_potential = coul
    call snapshot()
    negative_radius_times_potential(2) = coul(2)+3
    negative_radius_times_potential(3) = coul(3)-2
    actual = integrate_one_electron_radial(1,1,0)
    call near(actual-initial(1,1), 0.0_real64, "isolated fixture diagonal cancels")
    actual = integrate_one_electron_radial(1,2,0)
    call near(actual-initial(2,1), -6.0_real64, "isolated offdiagonal external integral")
#ifdef NPOL_MPI
    call MPI_Finalize(ierr)
#endif
    print "(A)", "PASS production consumers"
contains
    subroutine snapshot()
        integer :: jj, mm
        call build_direct_radial_potential(1)
        base_yp = orbital_y_potential
        large_exchange_source = 0
        small_exchange_source = 0
        call apply_off_diagonal_integrals()
        base_xp = large_exchange_source
        base_xq = small_exchange_source
        do jj = 1, 2
            do mm = 0, 1
                initial(jj,mm+1) = integrate_one_electron_radial(1,jj,mm)
            end do
        end do
    end subroutine snapshot
end program npol_consumers

! Deterministic test-only derivative/quadrature seams. Never linked into apps.
