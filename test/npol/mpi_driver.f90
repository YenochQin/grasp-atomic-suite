program npol_mpi_driver
    use mpi
    use npol_C
    use npol_mpi
    use npol_test_support
    use parameter_def, only: NNNP
    implicit none
    real(double) :: r(NNNP), zz(NNNP), coul(NNNP), reference(3), z
    integer :: rank, ierr, n, status, ios
    character(4096) :: root_path, directory
    character(32) :: mode
    character(NPOL_MESSAGE_LEN) :: message
    type(npol_config) :: expected
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call get_command_argument(1, mode)
    call get_command_argument(2, directory)
    ! Every worker has its own directory with no input or output file. Passing
    ! an existing invalid input to workers also detects accidental worker I/O.
    call chdir(trim(directory)//'/rank'//achar(iachar('0')+rank), ios)
    call require(ios == 0, 'rank cwd')
    root_path = 'worker-must-not-read.inp'
    if (rank == 0) root_path = trim(directory)//'/input/npol.inp'
    n = 3
    r = 0
    r(:3) = [0.0_double, 2.0_double, 4.0_double]
    coul = 7
    zz = coul
    z = 2
    if (mode == 'length' .and. rank /= 0) n = NNNP+1
    if (mode == 'mismatch' .and. rank /= 0) n = 2
    if (mode == 'charge' .and. rank /= 0) z = 3
    if (mode == 'compose' .and. rank /= 0) then
        r(2) = huge(1.0_double)
        coul(2) = huge(1.0_double)
        zz = coul
    end if
    call npol_init_mpi(trim(root_path), z, n, r, zz, MPI_COMM_WORLD, status, message)
    if (status /= 0) then
        call require(all(zz == coul), 'collective failure leaves all Coulomb arrays intact')
        call require(.not.npol_ready .and. .not.npol_settings%enabled, 'collective failure stays disabled')
        ! Reaching this collective proves that no worker was left in a Bcast.
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        if (rank == 0) print '(A)', 'EXPECTED_COLLECTIVE_FAILURE: '//trim(message)
        call npol_abort_mpi(MPI_COMM_WORLD, message)
    end if
    call require(mode == 'ok' .or. mode == 'disabled' .or. mode == 'missing' .or. mode == 'zero', &
        'failure scenario unexpectedly succeeded')
    call require(all(zz(4:) == 7), 'MPI untouched tail')
    if (mode == 'ok' .or. mode == 'zero') then
        expected = fixture()
        call require(npol_settings%schema_version == expected%schema_version, 'schema broadcast')
        call require(npol_settings%enabled .and. npol_settings%use_e1 .and. npol_settings%use_e2, 'switch broadcast')
        call require(npol_settings%isotope_z == 2 .and. npol_settings%isotope_a == 4, 'isotope broadcast')
        call require(npol_settings%units == expected%units, 'units broadcast')
        call require(npol_settings%model_id == expected%model_id, 'model broadcast')
        call require(npol_settings%parameter_source == expected%parameter_source, 'source broadcast')
        call require(npol_settings%b1 == 2 .and. npol_settings%b2 == 2, 'cutoff broadcast')
        if (mode == 'ok') then
            call require(npol_settings%k1 == 16 .and. npol_settings%k2 == 128, 'strength broadcast')
            call near(vnp_e1(2), -0.5_double, 'MPI E1 not divided by ranks')
            call near(vnp_e2(2), -1.0_double, 'MPI E2 not divided by ranks')
            call near(zz(2), 10.0_double, 'MPI full potential')
        else
            call require(npol_settings%k1 == 0 .and. npol_settings%k2 == 0, 'zero strength broadcast')
            call require(all(zz == coul), 'MPI zero exact copy')
        end if
    else
        call require(all(zz == coul) .and. .not.npol_settings%enabled, 'MPI disabled exact copy')
    end if
    reference = zz(:3)
    call MPI_Bcast(reference, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call require(all(reference == zz(:3)), 'rank identity')
    call MPI_Finalize(ierr)
    print '(A)', 'PASS MPI initialization'
end program
