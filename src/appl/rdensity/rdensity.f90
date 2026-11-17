! Shared entry point for rdensity and rdensity_mpi.
program rdensity
    use default_c, only: ndef
    use debug_c, only: cutoff
    use density_parallel, only: broadcast_options, density_fail, parallel_finish, parallel_start, rank, workers
    use density_calculation, only: calculate_density
    use getyn_i, only: getyn
    use setdbg_i, only: setdbg
    use setmc_i, only: setmc
    use setcon_i, only: setcon
    use setcsla_i, only: setcsla
    use getsmd_i, only: getsmd
    use getmixblock_i, only: getmixblock
    use factt_i, only: factt
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    character(24) :: name
    character(256) :: arg
    integer :: nci, ordering, status, i, ios, length, ncore
    logical :: yes, exists, interactive

    call parallel_start
    name = ""
    nci = 1
    ndef = 0
    ordering = 2
    status = 0
    cutoff = 1.0e-10_real64
    interactive = command_argument_count() == 0
    if (rank == 0) then
        if (interactive) then
            write(*,*) "Default settings?"
            yes = getyn()
            if (.not.yes) ndef = 1
            if (workers > 1 .and. ndef /= 0) then
                write(*,*) "MPI runs support default settings only, as in RHFS/RIS."
                status = 2
            else
                write(*,*) "Name of state"
                read(*,"(a)",iostat=ios) arg
                if (ios /= 0) status = 2
                call accept_name(arg)
                write(*,*) "Mixing coefficients from a CI calc.?"
                yes = getyn()
                if (yes) nci = 0
            end if
        else
            do i = 1,command_argument_count()
                call get_command_argument(i, arg, length=length, status=ios)
                if (ios /= 0 .or. length > len(arg)) then
                    status = 2
                    exit
                end if
                select case(trim(arg))
                case("--ci")
                    nci = 0
                case("--nonci")
                    nci = 1
                case("--order=1")
                    ordering = 1
                case("--order=2")
                    ordering = 2
                case("--help","-h")
                    status = 1
                    exit
                case default
                    if (i == 1 .and. arg(1:1) /= "-") then
                        call accept_name(arg)
                    else
                        write(*,*) "Unknown argument: ", trim(arg)
                        status = 2
                    end if
                end select
            end do
        end if
        if (status == 0 .and. len_trim(name) == 0) status = 2
        if (status /= 0) write(*,*) "Usage: rdensity[_mpi] name [--ci|--nonci] [--order=1|--order=2]"
    end if
    call broadcast_options(name, nci, ndef, ordering, status)
    if (status == 1) then
        call parallel_finish
        stop
    end if
    if (status /= 0) call density_fail("Invalid command line or interactive input")

    call require_file("isodata")
    call require_file(trim(name)//".c")
    call require_file(trim(name)//".w")
    if (nci == 0) then
        call require_file(trim(name)//".cm")
        arg = trim(name)//".cd"
    else
        call require_file(trim(name)//".m")
        arg = trim(name)//".d"
    end if
    if (rank == 0) then
        inquire(file=trim(arg), exist=exists)
        if (exists) call density_fail("Density output already exists: "//trim(arg))
    end if
    call setdbg
    call setmc
    call setcon
    ! Read-only inputs are shared; angular and orbital state is rank-local.
    call setcsla(name, ncore)
    call getsmd(name)
    call getmixblock(name, nci)
    call factt
    if (rank == 0 .and. interactive) then
        write(*,*) "How do you want to order your eigenvectors?"
        write(*,*) "1) Dominant component; 2) Decreasing occupation (preserving pure orbitals)"
        read(*,*,iostat=ios) ordering
        if (ios /= 0) status = 2
        if (ordering /= 1 .and. ordering /= 2) status = 2
    end if
    call broadcast_options(name, nci, ndef, ordering, status)
    if (status /= 0) call density_fail("Orbital ordering must be 1 or 2")
    call calculate_density(name, nci, ordering)
    call parallel_finish
contains
    subroutine accept_name(text)
        character(*), intent(in) :: text
        ! Legacy readers find the trailing blank in CHARACTER(24).
        if (len_trim(text) < 1 .or. len_trim(text) > 23) then
            status = 2
        else if (index(trim(text)," ") /= 0) then
            status = 2
        else
            name = trim(text)
        end if
    end subroutine accept_name

    subroutine require_file(filename)
        character(*), intent(in) :: filename
        logical :: found
        inquire(file=filename, exist=found)
        if (.not.found) call density_fail("Missing input file: "//filename)
    end subroutine require_file
end program rdensity
